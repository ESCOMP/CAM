module gmean_mod
   !-----------------------------------------------------------------------
   !
   ! Purpose:
   ! Perform mixed layer global calculations for energy conservation checks.
   !
   ! Method:
   ! Reproducible (semi-scalable):
   !    Gather to a master processor which sums all columns in order, or,
   !    gather to a few processors, each of which sums its block in order
   !    then sends its sum to the master processor for a final sum.
   !    Chunking is controlled by the namelist variable, <global_max_fields>.
   !    Sums should be reproducible as long as the value of this variable
   !    does not change.
   !
   !-----------------------------------------------------------------------
   use shr_kind_mod,  only: r8 => shr_kind_r8
   use perf_mod,      only: t_startf, t_stopf, t_adj_detailf
   use spmd_utils,    only: MPI_COMM_NULL, column_redist_t
   use cam_logfile,   only: iulog

   implicit none
   private
   save

   public :: gmean_init ! Initialize reproducible sum algorithm (if enabled)
   public :: gmean ! compute global mean of 2D fields on physics decomposition
   public :: gmean_finalize

   interface gmean
      module procedure gmean_arr
      module procedure gmean_scl
   end interface gmean

   ! Private data
   logical                       :: initialized = .false.
   logical                       :: use_repro_sum = .false.
   integer                       :: max_nflds = -1 ! # fields for disp. decomp
   type(column_redist_t)         :: column_reorder

   private :: test_gmean
   private :: bits_to_string
!!XXgoldyXX: v remove once new algorithm is in place
private :: gmean_float_repro
private :: gmean_fixed_repro
!!XXgoldyXX: ^ remove once new algorithm is in place

CONTAINS

   !
   !========================================================================
   !

   subroutine gmean_init(max_nflds_in)
      use phys_control,   only: phys_getopts
      use cam_abortutils, only: endrun, handle_allocate_error
      use spmd_utils,     only: masterprocid, masterproc, mpicom, npes, iam
      use spmd_utils,     only: MPI_INTEGER, MPI_REAL8, MPI_SUM, MPI_MAX
      use ppgrid,         only: pcols, begchunk, endchunk
      use phys_grid,      only: global_max_fields => phys_global_max_fields
      use phys_grid,      only: num_global_phys_cols, columns_on_task
      use phys_grid,      only: init_col_assem_p
      !-----------------------------------------------------------------------
      !
      ! Purpose:
      ! Pre-compute layout and communication information for efficient
      !  run-time compuatation of global sums
      !
      !-----------------------------------------------------------------------
      !
      ! Dummy arguments
      !
      integer,          intent(in)  :: max_nflds_in ! maximum number of fields
      ! Local variables
      !
      integer                       :: total_size    ! For memory diagnostic
      integer                       :: num_sum_tasks ! # PEs for one partial sum
      integer                       :: num_sum_blks  ! # partial sum blocks
      integer                       :: max_blck_size ! Largest partial sum block
      integer                       :: num_dest_tasks
      integer                       :: max_nflds_use
      integer                       :: dest_spacing
      integer                       :: next_field    ! next field to distribute
      integer                       :: next_start    ! next first col on task
      integer                       :: next_dtask    ! next destination task
      integer                       :: round_num
      integer                       :: findex
      integer                       :: color
      integer                       :: ierr
      integer                       :: max_fields    ! Global count
      integer                       :: block_sum_comm
      integer,          allocatable :: send_cnts(:)
      integer,          allocatable :: col_info(:)
      character(len=256)            :: errmsg
      integer,          parameter   :: max_recv = 16 ! max for gather
      integer,          parameter   :: max_pcnt = 3  ! max ratio for gather
      character(len=*), parameter   :: subname = 'gmean_init: '

      !
      !-----------------------------------------------------------------------
      !
      call phys_getopts(cam_repro_sum_out=use_repro_sum)
      if (.not. use_repro_sum) then
         call gmean_finalize()
         return
      end if
      if (initialized) then
         if (max_nflds /= max_nflds_in) then
            call gmean_finalize()
         else
            return ! Should not need to do anything
         end if
      end if
      ! Perform reproducible sum. Note, this depends on the value
      ! of <global_max_fields> not changing between runs
      max_nflds_use = max_nflds_in
      if (max_nflds_use < 1) then
         if (masterproc) then
            write(iulog, '(2a,i0,a)') subname, "WARNING, max_nflds_use = ", &
                 max_nflds_use, ", setting to 1 (minimum value)"
         end if
      end if
      total_size = 7 * storage_size(1)
      ! First try
      if (global_max_fields <= 0) then
         ! max global size not configured, go for one 2-D field so that the
         ! algorithm will work on 1PE
         max_fields = num_global_phys_cols
      else
         max_fields = global_max_fields
      end if
      if (max_fields < num_global_phys_cols) then
         ! We do not have enough space to fit an entire field on one PE
         ! so we have to do 'partial sums' on a number of PEs.
         num_sum_tasks = num_global_phys_cols / max_fields
         if ((num_sum_tasks * max_fields) < num_global_phys_cols) then
            num_sum_tasks = num_sum_tasks + 1
         end if
         num_sum_blks = 1
         ! Check to see if there are enough PEs for this choice of max_fields
         if (num_sum_tasks > npes) then
            write(errmsg, '(a,2(i0,a))') 'phys_global_max_fields = ',         &
                 max_fields, ' requires at least ', num_sum_tasks,            &
                 ' PEs for this resolution in cam_repro_sum mode'
            call endrun(subname//trim(errmsg))
         end if
         max_blck_size = num_global_phys_cols / num_sum_tasks
         if ((max_blck_size * num_sum_tasks) < num_global_phys_cols) then
            max_blck_size = max_blck_size + 1
         end if
         ! Since we are skimping on space, can we use more PEs?
         num_sum_blks = (npes / num_sum_tasks)
         column_reorder%max_nflds = num_sum_blks
      else
         ! We have room for at least one whole field on a task -- how many?
         num_sum_tasks = 1
         findex = max_fields / num_global_phys_cols
         max_blck_size = num_global_phys_cols * findex
         num_sum_blks = max_nflds_use / findex
         if ((num_sum_blks * findex) < max_nflds_use) then
            num_sum_blks = num_sum_blks + 1
         end if
         if (num_sum_blks > npes) then
            num_sum_blks = npes
            ! The MIN here should be redundant but I'm paranoid
            column_reorder%max_nflds = MIN(npes * findex, max_nflds_use)
         else
            column_reorder%max_nflds = max_nflds_use
         end if
      end if
      num_dest_tasks = num_sum_tasks * num_sum_blks
      if (num_dest_tasks > npes) then
         ! Something has gone off the rails, report out
         write(errmsg, '(2(a,i0))') "INTERNAL ERROR: num_dest_tasks = ",      &
              num_dest_tasks, ", npes = ", npes
         call endrun(subname//trim(errmsg))
      end if
      max_nflds_use = column_reorder%max_nflds
      ! If we get here, we should not have any allocated arrays, check.
      if (associated(column_reorder%dest_tasks)) then
         call endrun(subname//'dest_tasks allocated (should not be)')
      end if
      if (associated(column_reorder%col_starts)) then
         call endrun(subname//'col_starts allocated (should not be)')
      end if
      if (associated(column_reorder%num_rflds)) then
         call endrun(subname//'num_rflds allocated (should not be)')
      end if
      if (associated(column_reorder%task_sizes)) then
         call endrun(subname//'task_sizes allocated (should not be)')
      end if
      if (associated(column_reorder%recv_cnts)) then
         call endrun(subname//'recv_cnts allocated (should not be)')
      end if
      if (associated(column_reorder%recv_disps)) then
         call endrun(subname//'recv_disps allocated (should not be)')
      end if
      if (associated(column_reorder%recv_reorder)) then
         call endrun(subname//'recv_reorder allocated (should not be)')
      end if
      ! Allocate data about destination tasks (available on all pes)
      allocate(column_reorder%dest_tasks(0:num_dest_tasks-1), stat=ierr)
      call handle_allocate_error(ierr, subname, 'dest_tasks')
      total_size = total_size + (num_dest_tasks * storage_size(1))
      allocate(column_reorder%col_starts(0:num_dest_tasks-1), stat=ierr)
      call handle_allocate_error(ierr, subname, 'col_starts')
      total_size = total_size + (num_dest_tasks * storage_size(1))
      ! dest_spacing spaces out the collection tasks across the communicator
      dest_spacing = MAX(npes / num_dest_tasks, 1)
      do findex = 0, num_dest_tasks - 1
         column_reorder%dest_tasks(findex) = dest_spacing * findex ! < npes
      end do
      ! How many fields are gathered for each partial sum round?
      allocate(column_reorder%num_rflds(num_sum_blks), stat=ierr)
      call handle_allocate_error(ierr, subname, 'num_rflds')
      total_size = total_size + (num_sum_blks * storage_size(1))
      ! Distribute the partial summation columns over the receiving tasks
      next_field = 1
      next_dtask = 0
      round_num = 1
      column_reorder%num_rflds = 0
      column_reorder%my_nflds = 0
      do
         if (next_field <= column_reorder%max_nflds) then
            next_start = 1
            do findex = 0, num_sum_tasks - 1
               ! Sanity check
               if (next_field > column_reorder%max_nflds) then
                  write(errmsg, '(4(i0,a),i0)') iam,                          &
                       ') next_field out of range, ', next_field, ' > ',      &
                       column_reorder%max_nflds, ', index = ',                &
                       findex + 1, ' out of ', num_sum_tasks
                  call endrun(subname//trim(errmsg))
               end if
               column_reorder%col_starts(next_dtask) = next_start
               if (iam == column_reorder%dest_tasks(next_dtask)) then
                  column_reorder%strt_nfld = next_field
                  if (next_start == 1) then
                     column_reorder%my_nflds = column_reorder%my_nflds + 1
                  end if
               end if
               next_dtask = next_dtask + 1
               next_start = next_start + max_blck_size
               if (next_start > num_global_phys_cols) then
                  ! We are done with this field, start a new one
                  ierr = MAX(num_global_phys_cols / max_blck_size, 1)
                  if (round_num > SIZE(column_reorder%num_rflds, 1)) then
                     write(errmsg, '(2(i0,a),i0)') iam, ') round_num = ',     &
                          round_num, ', should be <= ',                       &
                          SIZE(column_reorder%num_rflds, 1)
                     call endrun(subname//trim(errmsg))
                  end if
                  if (ierr + next_field > column_reorder%max_nflds) then
                     column_reorder%num_rflds(round_num) =                     &
                          column_reorder%max_nflds - next_field + 1
                  else
                     column_reorder%num_rflds(round_num) = ierr
                  end if
                  round_num = round_num + 1
                  next_field = next_field + ierr
                  next_start = 1
               end if
            end do
         else
            exit
         end if
      end do
      column_reorder%num_rounds = COUNT(column_reorder%col_starts == 1)
      ! Sanity checks
      if (column_reorder%num_rounds /= num_sum_blks) then
         write(errmsg, '(2(i0,a),i0)') iam, ') num_rounds = ',                &
              column_reorder%num_rounds, ' should be ', num_sum_blks
         call endrun(subname//trim(errmsg))
      end if
      if (ANY(column_reorder%col_starts < 0)) then
         if (masterproc) then
            do findex = 0, num_dest_tasks - 1
               if (column_reorder%col_starts(findex) < 0)  then
                  write(iulog, '(a,i0,a)') 'column_reorder%col_starts(',      &
                       findex, ') < 0'
               end if
            end do
         end if
         call MPI_barrier(mpicom, ierr)
         write(errmsg, '(a,i0,a)') ' column_reorder%col_starts < 0 in ',      &
              COUNT(column_reorder%col_starts < 0), ' tasks'
         call endrun(subname//trim(errmsg))
      end if
      if (SUM(column_reorder%num_rflds) /= column_reorder%max_nflds) then
         write(errmsg, '(2(i0,a),i0)') iam, ') num_flds counting error, ',    &
              SUM(column_reorder%num_rflds), ', should be ',                  &
              column_reorder%max_nflds
         call endrun(subname//trim(errmsg))
      end if
      ! Get redistribution info from physics grid
      call init_col_assem_p(column_reorder, max_blck_size)
      ! Create communicator for receiving tasks
      block_sum_comm = mpicom
      if (ANY(column_reorder%dest_tasks == iam)) then
         color = 1
      else
         color = 2
      end if
      if (size(column_reorder%dest_tasks, 1) < npes) then
         call MPI_Comm_split(mpicom, color, iam, block_sum_comm, ierr)
      end if
      if (color == 1) then
         call mpi_comm_rank(block_sum_comm, column_reorder%recv_iam, ierr)
         if (masterproc) then
            column_reorder%recv_master_id = column_reorder%recv_iam
         end if
      else
         call mpi_comm_free(block_sum_comm, ierr)
         block_sum_comm = MPI_COMM_NULL
         if (masterproc) then
            ! Ack, we need masterproc to be a receiving task for broadcast
            errmsg = ': Internal Error, masterproc not a receiving task'
            call endrun(subname//trim(errmsg))
         end if
      end if
      column_reorder%mpi_comm = block_sum_comm
      ! Now that things are okay, broadcast out the correct receive master
      call MPI_bcast(column_reorder%recv_master_id, 1, MPI_INTEGER, &
           masterprocid, mpicom, ierr)
      if (column_reorder%mpi_comm == MPI_COMM_NULL) then
         column_reorder%recv_master_id = -1
      end if
      ! Retrieve number of columns (column_reorder%task_sizes) for each block
      allocate(column_reorder%recv_cnts(0:npes-1), stat=ierr)
      call handle_allocate_error(ierr, subname, '%recv_cnts')
      total_size = total_size + (npes * storage_size(1))
      column_reorder%recv_cnts = 0
      allocate(send_cnts(0:npes-1), stat=ierr)
      call handle_allocate_error(ierr, subname, 'send_cnts')
      total_size = total_size + (npes * storage_size(1))
      ! Initialize send_cnts to the send count for each destination task
      send_cnts = 0
      send_cnts(column_reorder%dest_tasks(:)) = column_reorder%task_sizes(:)
      if ( (num_dest_tasks >= max_recv) .or.                                  &
           ((npes / num_dest_tasks) < max_pcnt)) then
         call MPI_Alltoall(send_cnts, 1, MPI_INTEGER,                         &
              column_reorder%recv_cnts, 1, MPI_INTEGER, mpicom, ierr)
      else
         do findex = 0, num_dest_tasks - 1
            call MPI_gather(column_reorder%task_sizes(findex), 1,             &
                 MPI_INTEGER, column_reorder%recv_cnts, 1, MPI_INTEGER,       &
                 column_reorder%dest_tasks(findex), mpicom, ierr)
         end do
      end if
      ! Compute displacements from counts
      allocate(column_reorder%send_disps(0:npes - 1), stat=ierr)
      call handle_allocate_error(ierr, subname, '%send_disps')
      total_size = total_size + (npes * storage_size(1))
      column_reorder%send_disps(0) = 0
      allocate(column_reorder%recv_disps(0:npes - 1), stat=ierr)
      call handle_allocate_error(ierr, subname, '%recv_disps')
      total_size = total_size + (npes * storage_size(1))
      column_reorder%recv_disps(0) = 0
      do findex = 1, npes - 1
         column_reorder%send_disps(findex) =                                  &
              column_reorder%send_disps(findex - 1) + send_cnts(findex - 1)
         column_reorder%recv_disps(findex) =                                  &
              column_reorder%recv_disps(findex - 1) +                         &
              column_reorder%recv_cnts(findex - 1)
      end do
      ! Compute the summation order for each task with a block of data
      allocate(col_info(SUM(column_reorder%recv_cnts)), stat=ierr)
      call handle_allocate_error(ierr, subname, '%recv_cnts')
      allocate(column_reorder%recv_reorder(size(col_info, 1)), stat=ierr)
      call handle_allocate_error(ierr, subname, '%recv_reorder')
      col_info = 0
      total_size = total_size + (size(col_info, 1) * storage_size(1))
      call MPI_alltoallv(column_reorder%task_indices, send_cnts,              &
           column_reorder%send_disps, MPI_INTEGER, col_info,                  &
           column_reorder%recv_cnts, column_reorder%recv_disps,               &
           MPI_INTEGER, mpicom, ierr)
      ! col_info now has the global column index for each received quantity
      ierr = 0
      do findex = 1, size(col_info)
         ! color is the location where this column belongs in summation order
         color = col_info(findex) -                                           &
              column_reorder%col_starts(column_reorder%recv_iam) + 1
         if (color < 1) then
            write(iulog, '(i0,3a,i0,a,i0)') iam, ') ', subname,               &
                 'ERROR: recv_order underflow on column ', findex, ', ', color
            ierr = ierr + 1
         else if (color > size(column_reorder%recv_reorder)) then
            write(iulog, '(i0,2a,3(a,i0))') iam, ') ', subname,               &
                 'ERROR: recv_order overflow on column, ', findex, ', ',      &
                 color, ' should be between 1 and ',                          &
                 size(column_reorder%recv_reorder)
            ierr = ierr + 1
         else
            column_reorder%recv_reorder(color) = findex
         end if
      end do
      if (ierr > 0) then
         write(errmsg, '(a,i0,a)') ': ', ierr, ' recv_order errors'
         call endrun(subname//trim(errmsg))
      end if
      ! Done !
      max_nflds = max_nflds_use ! Bookkeeping for possible re-initialization
      if (max_nflds /= column_reorder%max_nflds) then
         write(errmsg, '(2(a,i0))') 'Internal Error, max_nflds = ',           &
              max_nflds, ' /= ', column_reorder%max_nflds
         call endrun(subname//trim(errmsg))
      end if
      initialized = .true.
      ! Cleanup
      deallocate(send_cnts)
      deallocate(col_info)
      ! Log some basic information on this repro-sum configuration
      if (masterproc) then
         write(iulog, *) '*****************************************'
         write(iulog, *) 'CAM Reproducible Global Sum Configuration'
         write(iulog, '(a,i0)') "  Max summation fields = ", max_nflds
         ! At this point, total_size is the size to store the algorithm
         write(iulog, '(2a)') "  Algorithm storage is ",                      &
              trim(bits_to_string(total_size))
         num_sum_tasks = total_size
         ! Runtime size
         total_size = (4 * npes) * storage_size(1) ! For MPI routines
         total_size = total_size + (storage_size(1.0_r8) *                    &
              ! Partial sum array
              ((column_reorder%max_nflds * SUM(column_reorder%recv_cnts)) +   &
              ! sizes for MPI_gatherv for final sum
              (size(column_reorder%dest_tasks, 1) * 2) +                      &
              ! Weighted field send array
              (MAXVAL(column_reorder%num_rflds) * columns_on_task *           &
              column_reorder%num_rounds)))
         write(iulog, '(2a)') "  Runtime max Memory per task  = ",            &
              trim(bits_to_string(total_size))
         write(iulog, '(2a)') "  Max Memory per task  = ",                    &
              trim(bits_to_string(total_size + num_sum_tasks))
         call MPI_comm_size(column_reorder%mpi_comm, total_size, ierr)
         write(iulog, '(a,i0,a)') "  Intermediate sum on ", &
              total_size, ' tasks'
         if (column_reorder%max_nflds < max_nflds_in) then
            total_size = max_nflds_in / column_reorder%max_nflds
            if (total_size * column_reorder%max_nflds < max_nflds_in) then
               total_size = total_size + 1
            end if
            write(iulog, '(2a,i0,a)') "  Sum is memory limited and may ", &
                 "require up to ", total_size, " summation rounds"
         else
            write(iulog, '(a,i0)') "  Sum is not memory limited"
         end if
         write(iulog, *) '*****************************************'
      end if
   end subroutine gmean_init

   subroutine gmean_arr(arr, arr_gmean, nflds)
      use shr_reprosum_mod, only: shr_reprosum_tolExceeded
      use shr_reprosum_mod, only: shr_reprosum_reldiffmax
      use shr_reprosum_mod, only: shr_reprosum_recompute
      use spmd_utils,     only: masterproc
      use cam_abortutils, only: endrun, handle_allocate_error
      use physconst,      only: pi
      use spmd_utils,     only: masterprocid, mpicom, npes, iam
      use spmd_utils,     only: MPI_REAL8, MPI_SUM, MPI_INTEGER
      use ppgrid,         only: pcols, begchunk, endchunk
      use phys_grid,      only: num_global_phys_cols, columns_on_task
      use phys_grid,      only: weighted_sum_p, weighted_field_p
      !-----------------------------------------------------------------------
      !
      ! Purpose:
      ! Compute the global mean of each field in "arr" in the physics
      ! chunked decomposition
      !
      !-----------------------------------------------------------------------
      !
      ! Dummy arguments
      !
      integer,           intent(in)  :: nflds            ! number of fields
      real(r8),          intent(in)  :: arr(pcols, begchunk:endchunk, nflds)
      real(r8),          intent(out) :: arr_gmean(nflds) ! global means
      !
      ! Local variables
      !
      integer                       :: fld_ind
      integer                       :: ierr
      ! relative diffs between 'fast' reproducible and nonreproducible means
      real(r8) :: rel_diff(2,nflds)
      real(r8)                      :: wsum(nflds)   ! Weighted sums
      real(r8),         parameter   :: inv4pi = 1.0_r8 / (4.0_r8 * pi)
      character(len=*), parameter   :: subname = 'gmean_arr: '
      !
      !-----------------------------------------------------------------------
      !
      call t_adj_detailf(+1)
      call t_startf('gmean_arr')
      arr_gmean(:) = 0.0_r8
      if (use_repro_sum) then
         call t_startf('gmean_fixed_repro')
         call gmean_fixed_repro(arr, arr_gmean, rel_diff, nflds)
         call t_stopf('gmean_fixed_repro')

         ! check that "fast" reproducible sum is accurate enough. If not,
         ! calculate using old method
         if (shr_reprosum_tolExceeded('gmean', nflds, masterproc,             &
              iulog, rel_diff)) then
            if (shr_reprosum_recompute) then
               do fld_ind = 1, nflds
                  if (rel_diff(1, fld_ind) > shr_reprosum_reldiffmax) then
                     call t_startf('gmean_float_repro')
                     call gmean_float_repro(arr(:,:,fld_ind),                 &
                          arr_gmean(fld_ind), 1)
                     call t_stopf('gmean_float_repro')
                  end if
               end do
            end if
         end if
         ! Normalize
         do fld_ind = 1, nflds
            arr_gmean(fld_ind) = arr_gmean(fld_ind) / (4.0 * pi)
         end do
      else
         ! Non reproducible algorithm
         do fld_ind = 1, nflds
            ! Get weighted sum from physics
            call weighted_sum_p(arr(:, :, fld_ind), wsum(fld_ind))
         end do
         call MPI_allreduce(wsum, arr_gmean, nflds, MPI_REAL8, MPI_SUM,       &
              mpicom, ierr)
         ! Normalize
         do fld_ind = 1, nflds
            arr_gmean(fld_ind) = arr_gmean(fld_ind) * inv4pi
         end do
      end if
      call t_stopf('gmean_arr')
      call t_adj_detailf(-1)
   end subroutine gmean_arr

   !
   !========================================================================
   !

   subroutine gmean_scl (arr, gmean)
      use ppgrid,        only: pcols, begchunk, endchunk
      !-----------------------------------------------------------------------
      !
      ! Purpose:
      ! Compute the global mean of each field in "arr" in the physics
      ! chunked decomposition
      !
      !-----------------------------------------------------------------------
      !
      ! Arguments
      !
      real(r8), intent(in)  :: arr(pcols,begchunk:endchunk)
      ! Input array, chunked
      real(r8), intent(out) :: gmean ! global means
      !
      ! Local workspace
      !
      integer, parameter    :: nflds = 1
      real(r8)              :: gmean_array(nflds)

      call gmean_arr(reshape(arr, (/ pcols, endchunk-begchunk+1, nflds /)),   &
           gmean_array, nflds)
      gmean = gmean_array(1)

   end subroutine gmean_scl

   !
   !========================================================================
   !

   subroutine gmean_finalize()
      !-----------------------------------------------------------------------
      !
      ! Purpose: Deallocate arrays and reset to uninitialized state
      !
      !-----------------------------------------------------------------------

      ! Local variable
      integer :: ierr

      initialized = .false.
      use_repro_sum = .false.
      column_reorder%recv_iam = -1
      column_reorder%recv_master_id = -1
      column_reorder%max_nflds = 0
      column_reorder%num_rounds = 0
      column_reorder%strt_nfld = -1
      column_reorder%my_nflds = 0
      if (column_reorder%mpi_comm /= MPI_COMM_NULL) then
         call MPI_comm_free(column_reorder%mpi_comm, ierr)
         column_reorder%mpi_comm = MPI_COMM_NULL
      end if
      if (associated(column_reorder%dest_tasks)) then
         deallocate(column_reorder%dest_tasks)
         nullify(column_reorder%dest_tasks)
      end if
      if (associated(column_reorder%col_starts)) then
         deallocate(column_reorder%col_starts)
         nullify(column_reorder%col_starts)
      end if
      if (associated(column_reorder%num_rflds)) then
         deallocate(column_reorder%num_rflds)
         nullify(column_reorder%num_rflds)
      end if
      if (associated(column_reorder%recv_cnts)) then
         deallocate(column_reorder%recv_cnts)
         nullify(column_reorder%recv_cnts)
      end if
      if (associated(column_reorder%recv_disps)) then
         deallocate(column_reorder%recv_disps)
         nullify(column_reorder%recv_disps)
      end if
      if (associated(column_reorder%recv_reorder)) then
         deallocate(column_reorder%recv_reorder)
         nullify(column_reorder%recv_reorder)
      end if
      if (associated(column_reorder%task_sizes)) then
         deallocate(column_reorder%task_sizes)
         nullify(column_reorder%task_sizes)
      end if
      if (associated(column_reorder%task_indices)) then
         deallocate(column_reorder%task_indices)
         nullify(column_reorder%task_indices)
      end if
      if (associated(column_reorder%send_disps)) then
         deallocate(column_reorder%send_disps)
         nullify(column_reorder%send_disps)
      end if
      if (associated(column_reorder%send_reorder)) then
         deallocate(column_reorder%send_reorder)
         nullify(column_reorder%send_reorder)
      end if
   end subroutine gmean_finalize

   !
   !========================================================================
   !

   subroutine test_gmean(max_diff)
      use physconst,      only: pi
      use spmd_utils,     only: iam, masterproc
      use cam_abortutils, only: endrun
      use cam_logfile,    only: iulog
      use phys_grid,      only: get_chunk_info_p, get_gcol_p, get_area_p
      use phys_grid,      only: pgcols => num_global_phys_cols
      use phys_grid,      only: columns_on_task
      use ppgrid,         only: pcols, begchunk, endchunk

      ! Dummy argument
      real(r8), optional, intent(in) :: max_diff
      ! Local variables
      integer                     :: lchnk, icol, gcol
      integer                     :: findex, col_ind
      real(r8)                    :: test_arr(pcols, begchunk:endchunk, 3)
      real(r8)                    :: test_mean(3)
      real(r8)                    :: expect(3)
      real(r8)                    :: diff
      real(r8)                    :: max_diff_use
      real(r8),         parameter :: fact2 = 0.01_r8
      real(r8),         parameter :: pi4 = 4.0_r8 * pi
      real(r8),         parameter :: max_diff_def = 1.0e-14_r8
      character(len=256)          :: errmsg
      character(len=*), parameter :: subname = 'test_gmean: '

      if (present(max_diff)) then
         max_diff_use = max_diff
      else
         max_diff_use = max_diff_def
      end if
      do col_ind = 1, columns_on_task
         call get_chunk_info_p(col_ind, lchnk, icol)
         gcol = get_gcol_p(lchnk, icol)
         test_arr(icol,lchnk,1) = 1.0_r8
         test_arr(icol,lchnk,2) = real(gcol, r8) * pi4 / get_area_p(lchnk, icol)
         test_arr(icol,lchnk,3) = test_arr(icol,lchnk,2) * fact2
      end do
      test_mean(:) = (/ -1.0_r8, -2.71828_r8, -3.14159_r8 /)
      expect(1) = 1.0_r8
      expect(2) = real((pgcols + 1) * pgcols / 2, r8)
      expect(3) = expect(2) * fact2
      call gmean(test_arr, test_mean, 3)
      do findex = 1, 3
         diff = abs(test_mean(findex) - expect(findex)) / expect(findex)
         if (diff > max_diff_use) then
            write(errmsg, '(i0,a,i0,2(a,e20.13e2))') iam, ': test_mean(',     &
                 findex, ') FAIL: ', test_mean(findex), ' /= ', expect(findex)
            call endrun(subname//trim(errmsg))
         end if
      end do
      if (masterproc) then
         write(iulog, *) subname, test_mean(:)
      end if
   end subroutine test_gmean

   !
   !========================================================================
   !
   function bits_to_string(num_bits)
      integer, intent(in) :: num_bits
      character(len=32)   :: bits_to_string
      ! Local variables
      real(r8) :: value

      ! Start with bytes
      value = real(num_bits, r8) / 8.0_r8
      ! Maybe convert to higher units
      if (value > (1024.0_r8 * 1024.0_r8)) then
         value = value / (1024.0_r8 * 1024.0_r8)
         write(bits_to_string, '(F8.2," Mbytes")') value
      else if (value > 1024.0_r8) then
         value = value / 1024.0_r8
         write(bits_to_string, '(F8.2," Kbytes")') value
      else
         write(bits_to_string, '(F8.2," bytes")') value
      end if
   end function bits_to_string

!!XXgoldyXX: v remove once new algorithm is in place
   subroutine gmean_float_repro(arr, arr_gmean, nflds)
!-----------------------------------------------------------------------
!
! Purpose:
! Compute the global mean of each field in "arr" in the physics
! chunked decomposition - all work is done on the masterproc to avoid
! order of operations differences and assure bfb reproducibility.
!
!-----------------------------------------------------------------------

      use cam_abortutils, only: endrun
      use spmd_utils,     only: iam, masterproc, npes, masterprocid, mpicom
      use spmd_utils,     only: MPI_INTEGER, MPI_REAL8
      use ppgrid,         only: pcols, begchunk, endchunk
      use phys_grid,      only: get_chunk_info_p, get_gcol_p
      use phys_grid,      only: pgcols => num_global_phys_cols
      use phys_grid,      only: columns_on_task, weighted_field_p

!
! Arguments
!
      integer,  intent(in)  :: nflds
      real(r8), intent(in)  :: arr(pcols,begchunk:endchunk,nflds)
      real(r8), intent(out) :: arr_gmean(nflds)
!
! Local workspace
!
      integer                       :: itemp
      integer                       :: lchnk, icol
      integer                       :: ierr           ! MPI error return
      integer                       :: col_ind, findex
      integer,          allocatable :: gcols(:)
      integer,          allocatable :: ind_order(:)   ! global index ordering
      integer,          allocatable :: ccounts(:)     ! per-task Column counts
      integer,          allocatable :: disps(:)       ! task displacements
      integer,          allocatable :: lgcols(:)      ! global col # on task
      real(r8),         allocatable :: wght_fld(:,:)  ! weighted <arr>
      real(r8),         allocatable :: loc_field(:,:) ! Local field info
      real(r8),         allocatable :: arr_field(:,:) ! Global field info
      character(len=128)            :: errmsg
      character(len=*), parameter   :: subname = 'gmean_float_repro: '
!
!-----------------------------------------------------------------------
!
      if (masterproc) then
         allocate(arr_field(pgcols, nflds))
         allocate(gcols(pgcols))
         allocate(ind_order(pgcols))
         allocate(disps(npes))
         allocate(ccounts(npes))
      end if
      allocate(lgcols(columns_on_task))
      allocate(loc_field(nflds, columns_on_task))

      arr_field(:,:) = 0.0_r8
      arr_gmean = 0.0_r8

      ! We need global column number for each task
      do col_ind = 1, columns_on_task
         call get_chunk_info_p(col_ind, lchnk, icol)
         lgcols(col_ind) = get_gcol_p(lchnk, icol)
      end do
      ! Gather displacements from each tasks
      ind_order = 0
      call MPI_gather(columns_on_task, 1, MPI_INTEGER, ccounts, 1,            &
           MPI_INTEGER, masterprocid, mpicom, ierr)
      if (masterproc) then
         ! Check that disps adds up.
         if (SUM(ccounts) /= pgcols) then
            write(errmsg, '(a,2(i0,a))') 'ccounts (', SUM(ccounts),           &
                 ') /= pgcols (', pgcols, ')'
            call endrun(subname//trim(errmsg))
         end if
         ! Create a displacement array
         disps(1) = 0
         do col_ind = 2, npes
            disps(col_ind) = disps(col_ind - 1) + ccounts(col_ind - 1)
         end do
      end if
      ! Now, gather up all global columns
      call MPI_gatherv(lgcols, columns_on_task, MPI_INTEGER,                  &
           gcols, ccounts, disps, MPI_INTEGER, masterprocid, mpicom, ierr)
      ! Create an indirect addressing array for global addition order
      if (masterproc) then
         ind_order = -1
         do col_ind = 1, pgcols
            ind_order(gcols(col_ind)) = col_ind
         end do
      end if
      ! Compute the weighted field on each task
      call weighted_field_p(arr, wght_fld)
      ! Gather data onto masterproc
      if (masterproc) then
         do col_ind = 1, npes
            disps(col_ind) = disps(col_ind) * nflds
            ccounts(col_ind) = ccounts(col_ind) * nflds
         end do
      end if
      call MPI_gatherv(wght_fld, columns_on_task, MPI_REAL8,                  &
           loc_field, ccounts, disps, MPI_REAL8, masterprocid, mpicom, ierr)
      ! Add everything up in global_order
      do col_ind = 1, pgcols
         itemp = ind_order(col_ind)
         do findex = 1, nflds
            arr_gmean(findex) = arr_gmean(findex) + loc_field(findex, itemp)
         end do
      end do

      ! Finally, broadcast the sums to all tasks
      call MPI_bcast(arr_gmean, nflds, MPI_REAL8, masterprocid, mpicom, ierr)

      ! Cleanup
      if (allocated(arr_field)) then
         deallocate(arr_field)
      end if
      if (allocated(gcols)) then
         deallocate(gcols)
      end if
      if (allocated(ccounts)) then
         deallocate(ccounts)
      end if
      if (allocated(ind_order)) then
         deallocate(ind_order)
      end if
      if (allocated(disps)) then
         deallocate(disps)
      end if
      if (allocated(lgcols)) then
         deallocate(lgcols)
      end if
      if (allocated(loc_field)) then
         deallocate(loc_field)
      end if
      if (allocated(wght_fld)) then
         deallocate(wght_fld)
      end if

   end subroutine gmean_float_repro

!
!========================================================================
!
   subroutine gmean_fixed_repro (arr, arr_gmean, rel_diff, nflds)
      use shr_reprosum_mod, only: shr_reprosum_calc
      use spmd_utils,       only: mpicom
      use ppgrid,           only: pcols, begchunk, endchunk
      use phys_grid,        only: pgcols => num_global_phys_cols
      use phys_grid,        only: nlcols => columns_on_task
      use phys_grid,        only: get_ncols_p, get_wght_all_p
!-----------------------------------------------------------------------
!
! Purpose:
! Compute the global mean of each field in "arr" in the physics
! chunked decomposition with a reproducible yet scalable implementation
! based on a fixed-point algorithm.
!
!-----------------------------------------------------------------------
      !
      ! Dummy Arguments
      !
      integer,  intent(in)  :: nflds                              ! num fields
      real(r8), intent(in)  :: arr(pcols,begchunk:endchunk,nflds) ! Input array
      real(r8), intent(out) :: arr_gmean(nflds)                   ! global means
      ! rel_diff = relative and absolute differences between
      !            reproducible and nonreproducible means
      real(r8), intent(out) :: rel_diff(2,nflds)
      !
      ! Local workspace
      !
      integer               :: lchnk, icol ! chunk, column indices
      integer               :: ifld        ! field index
      integer               :: ncols       ! num columns in current chunk
      integer               :: count       ! summand count
      integer               :: ierr        ! MPI error return

      real(r8)              :: wght(pcols) ! column for integration weights
      real(r8), allocatable :: xfld(:,:)   ! weighted summands
      !
      !-----------------------------------------------------------------------
      !
      allocate(xfld(nlcols, nflds))

      ! pre-weight summands
      do ifld = 1, nflds
         count = 0
         do lchnk = begchunk, endchunk
            ncols = get_ncols_p(lchnk)
            call get_wght_all_p(lchnk, ncols, wght)
            do icol = 1, ncols
               count = count + 1
               xfld(count, ifld) = arr(icol, lchnk, ifld) * wght(icol)
            end do
         end do
      end do

      ! call fixed-point algorithm
      call shr_reprosum_calc (xfld, arr_gmean, count, nlcols, nflds,          &
           gbl_count=pgcols, commid=mpicom, rel_diff=rel_diff)

      deallocate(xfld)

   end subroutine gmean_fixed_repro
!!XXgoldyXX: ^ remove once new algorithm is in place

end module gmean_mod
