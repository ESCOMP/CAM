!-------------------------------------------------------------------------------
! This localizes ESMF regridding operations to allow for multiple instances of
! CAM.
!-------------------------------------------------------------------------------
module edyn_grid_comp
   use shr_kind_mod,   only: r8 => shr_kind_r8, cs=>shr_kind_cs, cl=>shr_kind_cl
   use ESMF,           only: ESMF_KIND_I4, ESMF_Mesh, ESMF_DistGrid
   use ESMF,           only: ESMF_State, ESMF_Clock, ESMF_GridComp
   use ppgrid,         only: pcols
   use cam_logfile,    only: iulog
   use shr_sys_mod,    only: shr_sys_flush
   use cam_abortutils, only: endrun

   implicit none

   private

   public :: edyn_grid_comp_init
   public :: edyn_grid_comp_run1
   public :: edyn_grid_comp_run2
   public :: edyn_grid_comp_final

   ! Private data and interfaces
   ! phys_mesh: Local copy of physics grid
   type(ESMF_Mesh)     :: phys_mesh
   ! edyn_comp: ESMF gridded component for the ionosphere models
   type(ESMF_GridComp) :: phys_comp
   ! Local copy of ionosphere epotential model
   character(len=16)   :: ionos_epotential_model = 'none'
   ! Total number of columns on this task
   integer             :: total_cols = 0
   integer             :: col_start = 1
   integer             :: col_end = -1
   integer             :: nlev = 0
   ! dist_grid_2d: DistGrid for 2D fields
   type(ESMF_DistGrid) :: dist_grid_2d
   ! Which run?
   integer             :: do_run
   ! Pointers for run1 output
   real(r8), pointer   :: prescr_efx_phys(:) => NULL()
   real(r8), pointer   :: prescr_kev_phys(:) => NULL()
   logical             :: ionos_epotential_amie
   logical             :: ionos_epotential_ltr
   ! Pointers for run2
   real(r8), pointer   :: omega_blck(:,:) => NULL()
   real(r8), pointer   :: pmid_blck(:,:) => NULL()
   real(r8), pointer   :: zi_blck(:,:) => NULL()
   real(r8), pointer   :: hi_blck(:,:) => NULL()
   real(r8), pointer   :: u_blck(:,:) => NULL()
   real(r8), pointer   :: v_blck(:,:) => NULL()
   real(r8), pointer   :: tn_blck(:,:) => NULL()
   real(r8), pointer   :: sigma_ped_blck(:,:) => NULL()
   real(r8), pointer   :: sigma_hall_blck(:,:) => NULL()
   real(r8), pointer   :: te_blck(:,:) => NULL()
   real(r8), pointer   :: ti_blck(:,:) => NULL()
   real(r8), pointer   :: mbar_blck(:,:) => NULL()
   real(r8), pointer   :: n2mmr_blck(:,:) => NULL()
   real(r8), pointer   :: o2mmr_blck(:,:) => NULL()
   real(r8), pointer   :: o1mmr_blck(:,:) => NULL()
   real(r8), pointer   :: o2pmmr_blck(:,:) => NULL()
   real(r8), pointer   :: nopmmr_blck(:,:) => NULL()
   real(r8), pointer   :: n2pmmr_blck(:,:) => NULL()
   real(r8), pointer   :: opmmr_blck(:,:) => NULL()
   real(r8), pointer   :: opmmrtm1_blck(:,:) => NULL()
   real(r8), pointer   :: ui_blck(:,:) => NULL()
   real(r8), pointer   :: vi_blck(:,:) => NULL()
   real(r8), pointer   :: wi_blck(:,:) => NULL()
   real(r8)            :: rmassO2p
   real(r8)            :: rmassNOp
   real(r8)            :: rmassN2p
   real(r8)            :: rmassOp

CONTAINS

   subroutine edyn_gcomp_init(comp, importState, exportState, clock, rc)
      use ESMF,         only: ESMF_DistGridCreate, ESMF_MeshCreate
      use ESMF,         only: ESMF_FILEFORMAT_ESMFMESH,  ESMF_MeshGet
      use cam_instance, only: inst_name
      use phys_control, only: phys_getopts
      use phys_grid,    only: get_ncols_p, get_gcol_p, get_rlon_all_p, get_rlat_all_p
      use ppgrid,       only: begchunk, endchunk
      use edyn_esmf,    only: edyn_esmf_chkerr, edyn_esmf_update_phys_mesh
      use shr_const_mod,only: shr_const_pi

      ! Dummy arguments
      type(ESMF_GridComp)  :: comp
      type(ESMF_State)     :: importState
      type(ESMF_State)     :: exportState
      type(ESMF_Clock)     :: clock
      integer, intent(out) :: rc

      ! Local variables
      integer                               :: ncols
      integer                               :: chnk, col, dindex
      integer,                allocatable   :: decomp(:)
      character(len=cl)                     :: grid_file
      character(len=*),           parameter :: subname = 'edyn_gcomp_init'
      real(r8)        ,           parameter :: radtodeg = 180.0_r8/shr_const_pi
      integer                 :: spatialDim
      integer                 :: numOwnedElements
      real(r8), pointer       :: ownedElemCoords(:)
      real(r8), pointer       :: lat(:), latMesh(:)
      real(r8), pointer       :: lon(:), lonMesh(:)
      real(r8)                :: lats(pcols)                       ! array of chunk latitudes
      real(r8)                :: lons(pcols)                       ! array of chunk longitude
      integer :: i, c, n
      character(len=cs)       :: tempc1,tempc2
      character(len=300)      :: errstr

      real(r8), parameter :: abstol =  1.e-6_r8

      ! Find the physics grid file
      call phys_getopts(physics_grid_out=grid_file)
      ! Compute the local decomp
      total_cols = 0
      do chnk = begchunk, endchunk
         total_cols = total_cols + get_ncols_p(chnk)
      end do
      allocate(decomp(total_cols))
      dindex = 0
      do chnk = begchunk, endchunk
         ncols = get_ncols_p(chnk)
         do col = 1, ncols
            dindex = dindex + 1
            decomp(dindex) = get_gcol_p(chnk, col)
         end do
      end do
      ! Create a DistGrid based on the physics decomp
      dist_grid_2d = ESMF_DistGridCreate(arbSeqIndexList=decomp, rc=rc)
      call edyn_esmf_chkerr(subname, 'ESMF_DistGridCreate phys decomp', rc)
      ! Create an ESMF_mesh for the physics decomposition
      phys_mesh = ESMF_MeshCreate(trim(grid_file), ESMF_FILEFORMAT_ESMFMESH,  &
           elementDistgrid=dist_grid_2d, rc=rc)
      call edyn_esmf_chkerr(subname, 'ESMF_MeshCreateFromFile', rc)
      call edyn_esmf_update_phys_mesh(phys_mesh)
      do_run = 1


      ! Check that the mesh coordinates are consistent with the model physics column coordinates

      ! obtain mesh lats and lons
      call ESMF_MeshGet(phys_mesh, spatialDim=spatialDim, numOwnedElements=numOwnedElements, rc=rc)
      call edyn_esmf_chkerr(subname, 'ESMF_MeshGet', rc)

      if (numOwnedElements /= total_cols) then
         write(tempc1,'(i10)') numOwnedElements
         write(tempc2,'(i10)') total_cols
         call endrun(trim(subname)//": ERROR numOwnedElements "// &
                     trim(tempc1) //" not equal to local size "// trim(tempc2))
      end if

      allocate(ownedElemCoords(spatialDim*numOwnedElements))
      allocate(lonMesh(total_cols), latMesh(total_cols))
      call ESMF_MeshGet(phys_mesh, ownedElemCoords=ownedElemCoords)

      do n = 1,total_cols
         lonMesh(n) = ownedElemCoords(2*n-1)
         latMesh(n) = ownedElemCoords(2*n)
      end do

      ! obtain internally generated cam lats and lons
      allocate(lon(total_cols)); lon(:) = 0._r8
      allocate(lat(total_cols)); lat(:) = 0._r8
      n=0
      do c = begchunk, endchunk
         ncols = get_ncols_p(c)
         ! latitudes and longitudes returned in radians
         call get_rlat_all_p(c, ncols, lats)
         call get_rlon_all_p(c, ncols, lons)
         do i=1,ncols
            n = n+1
            lat(n) = lats(i)*radtodeg
            lon(n) = lons(i)*radtodeg
         end do
      end do

      errstr = ''
      ! error check differences between internally generated lons and those read in
      do n = 1,total_cols
         if (abs(lonMesh(n) - lon(n)) > abstol) then
            if ( (abs(lonMesh(n)-lon(n)) > 360._r8+abstol) .or. (abs(lonMesh(n)-lon(n)) < 360._r8-abstol) ) then
               write(errstr,100) n,lon(n),lonMesh(n), abs(lonMesh(n)-lon(n))
               write(iulog,*) trim(errstr)
            endif
         end if
         if (abs(latMesh(n) - lat(n)) > abstol) then
            ! poles in the 4x5 SCRIP file seem to be off by 1 degree
            if (.not.( (abs(lat(n))>88.0_r8) .and. (abs(latMesh(n))>88.0_r8) )) then
               write(errstr,101) n,lat(n),latMesh(n), abs(latMesh(n)-lat(n))
               write(iulog,*) trim(errstr)
            endif
         end if
      end do

      if ( len_trim(errstr) > 0 ) then
         call endrun(subname//': physics mesh coords do not match model coords')
      end if

      ! deallocate memory
      deallocate(ownedElemCoords)
      deallocate(lon, lonMesh)
      deallocate(lat, latMesh)
      deallocate(decomp)

100   format('edyn_gcomp_init: coord mismatch... n, lon(n), lonmesh(n), diff_lon = ',i6,2(f21.13,3x),d21.5)
101   format('edyn_gcomp_init: coord mismatch... n, lat(n), latmesh(n), diff_lat = ',i6,2(f21.13,3x),d21.5)

   end subroutine edyn_gcomp_init

   !-----------------------------------------------------------------------
   subroutine edyn_gcomp_run(comp, importState, exportState, clock, rc)
      use ESMF,              only: ESMF_SUCCESS, ESMF_Array, ESMF_ArrayGet
      use ESMF,              only: ESMF_StateGet
      use epotential_params, only: epot_crit_colats
      use edyn_esmf,         only: edyn_esmf_chkerr
      use dpie_coupling,     only: d_pie_epotent
      use dpie_coupling,     only: d_pie_coupling

      ! Dummy arguments
      type(ESMF_GridComp)  :: comp
      type(ESMF_State)     :: importState
      type(ESMF_State)     :: exportState
      type(ESMF_Clock)     :: clock
      integer, intent(out) :: rc
      ! Local variables
      type(ESMF_Array)                    :: run_type
      integer                             :: cols, cole, blksize
      character(len=cs)                   :: errmsg
      character(len=*), parameter         :: subname = 'edyn_gcomp_run'

      if (do_run == 1) then
         if ( ionos_epotential_amie .or. ionos_epotential_ltr) then
            call d_pie_epotent(ionos_epotential_model, epot_crit_colats,      &
                 cols=col_start, cole=col_end,                                &
                 efx_phys=prescr_efx_phys, kev_phys=prescr_kev_phys,          &
                 amie_in=ionos_epotential_amie, ltr_in=ionos_epotential_ltr )
         else
            call d_pie_epotent(ionos_epotential_model, epot_crit_colats)
         end if
      else if (do_run == 2) then
         call d_pie_coupling(omega_blck, pmid_blck, zi_blck, hi_blck,         &
              u_blck, v_blck, tn_blck, sigma_ped_blck, sigma_hall_blck,       &
              te_blck, ti_blck, mbar_blck, n2mmr_blck, o2mmr_blck,            &
              o1mmr_blck, o2pmmr_blck, nopmmr_blck, n2pmmr_blck,              &
              opmmr_blck, opmmrtm1_blck, ui_blck, vi_blck, wi_blck,           &
              rmassO2p, rmassNOp, rmassN2p, rmassOp, col_start, col_end, nlev)
      else
         write(errmsg, '(2a,i0)') subname, ': Unknown run number, ', do_run
         call endrun(trim(errmsg))
      end if

      rc = ESMF_SUCCESS

   end subroutine edyn_gcomp_run
   !-----------------------------------------------------------------------
   subroutine edyn_gcomp_final(comp, importState, exportState, clock, rc)
      use ESMF,      only: ESMF_MeshDestroy
      use ESMF,      only: ESMF_SUCCESS
      use edyn_esmf, only: edyn_esmf_chkerr

      ! Dummy arguments
      type(ESMF_GridComp)  :: comp
      type(ESMF_State)     :: importState
      type(ESMF_State)     :: exportState
      type(ESMF_Clock)     :: clock
      integer, intent(out) :: rc
      ! Local variables
      character(len=*), parameter :: subname = 'edyn_gcomp_final'

      call ESMF_MeshDestroy(phys_mesh, rc=rc)
      rc = ESMF_SUCCESS

   end subroutine edyn_gcomp_final

   !-----------------------------------------------------------------------
   subroutine edyn_gcomp_SetServices(comp, rc)
      use ESMF,         only: ESMF_GridCompSetEntryPoint
      use ESMF,         only: ESMF_METHOD_INITIALIZE, ESMF_METHOD_RUN
      use ESMF,         only: ESMF_METHOD_FINALIZE, ESMF_SUCCESS
      use edyn_esmf,    only: edyn_esmf_chkerr

      type(ESMF_GridComp)  :: comp
      integer, intent(out) :: rc
      character(len=*), parameter :: subname = 'edyn_gcomp_SetServices'

      ! Set the entry points for standard ESMF Component methods
      call ESMF_GridCompSetEntryPoint(comp, ESMF_METHOD_INITIALIZE, &
           userRoutine=edyn_gcomp_Init, rc=rc)
      call edyn_esmf_chkerr(subname, 'ESMF_GridCompSetEntryPoint init', rc)
      call ESMF_GridCompSetEntryPoint(comp, ESMF_METHOD_RUN, &
           userRoutine=edyn_gcomp_Run, rc=rc)
      call edyn_esmf_chkerr(subname, 'ESMF_GridCompSetEntryPoint run', rc)
      call ESMF_GridCompSetEntryPoint(comp, ESMF_METHOD_FINALIZE, &
           userRoutine=edyn_gcomp_Final, rc=rc)
      call edyn_esmf_chkerr(subname, 'ESMF_GridCompSetEntryPoint final', rc)

   end subroutine edyn_gcomp_SetServices

   subroutine edyn_grid_comp_init(mpi_comm)
      use mpi,          only: MPI_INTEGER
      use ESMF,         only: ESMF_StateCreate, ESMF_GridCompInitialize
      use ESMF,         only: ESMF_GridCompCreate, ESMF_GridCompSetServices
      use ESMF,         only: ESMF_VM, ESMF_VMGetCurrent, ESMF_VMGet
      use cam_instance, only: inst_index, inst_name
      use edyn_esmf,    only: edyn_esmf_chkerr

      ! Dummy argument
      integer, intent(in)           :: mpi_comm
      ! Local variables
      integer,          allocatable :: petlist(:)
      integer                       :: iam
      integer                       :: npes
      integer                       :: localPet
      integer                       :: petCount
      integer                       :: rc
      type(ESMF_VM)                 :: vm_init
      character(len=*), parameter   :: subname = 'edyn_grid_comp_init'

      !! Gather PE information for this instance
      call ESMF_VMGetCurrent(vm_init, rc=rc)
      call edyn_esmf_chkerr(subname, 'ESMF_VMGetCurrent', rc)
      call ESMF_VMGet(vm_init, localPet=localPet, petCount=petCount)
      call edyn_esmf_chkerr(subname, 'ESMF_VMGet', rc)
      call mpi_comm_size(mpi_comm, npes, rc)
      call mpi_comm_rank(mpi_comm, iam, rc)
      ! Collect all the PETS for each instance for phys grid
      allocate(petlist(npes))
      call mpi_allgather(localPet, 1, MPI_INTEGER, petlist, 1, MPI_INTEGER, mpi_comm, rc)
      ! Now, we should be able to create a gridded component
      phys_comp = ESMF_GridCompCreate(name=trim(inst_name), petList=petlist, rc=rc)
      call edyn_esmf_chkerr(subname, 'ESMF_GridCompCreate '//trim(inst_name), rc)
      call ESMF_GridCompSetServices(phys_comp, edyn_gcomp_SetServices, rc=rc)
      call edyn_esmf_chkerr(subname, 'ESMF_GridCompSetServices '//trim(inst_name), rc)
      ! Initialize the required component arguments
      call ESMF_GridCompInitialize(phys_comp, rc=rc)
      call edyn_esmf_chkerr(subname, 'ESMF_GridCompInitialize', rc)

   end subroutine edyn_grid_comp_init

   subroutine edyn_grid_comp_run1(ionos_epotential_model_in,  &
              cols, cole, efx_phys, kev_phys, amie_in, ltr_in)

      use ESMF,         only: ESMF_GridCompRun
      use edyn_esmf,    only: edyn_esmf_chkerr

      ! Dummy arguments
      character(len=*),           intent(in) :: ionos_epotential_model_in
      integer,          optional, intent(in) :: cols
      integer,          optional, intent(in) :: cole
      real(r8),         optional, target, intent(out) :: efx_phys(:)
      real(r8),         optional, target, intent(out) :: kev_phys(:)
      logical,          optional, intent(in) :: amie_in
      logical,          optional, intent(in) :: ltr_in

      ! Local variables
      integer                                :: rc
      character(len=*), parameter            :: subname = 'edyn_grid_comp_run1'
      logical :: args_present(6)

      do_run = 1
      args_present(:) = (/ present(cols), present(cole),  present(efx_phys), present(kev_phys), &
                           present(amie_in), present(ltr_in) /)

      if ( any( args_present ) ) then
         if (.not. all( args_present ) ) then
            call endrun(subname//': all optional arguments must be present for AMIE/LTR')
         endif

         ionos_epotential_amie = amie_in
         ionos_epotential_ltr = ltr_in
         prescr_efx_phys => efx_phys
         prescr_kev_phys => kev_phys
         col_start = cols
         col_end = cole
      else
         ! No else check assume no optional arguments are passed
         nullify(prescr_efx_phys)
         nullify(prescr_kev_phys)
      end if
      ionos_epotential_model = ionos_epotential_model_in
      call ESMF_GridCompRun(phys_comp, rc=rc)
      call edyn_esmf_chkerr(subname, 'ESMF_GridCompRun', rc)

   end subroutine edyn_grid_comp_run1

   subroutine edyn_grid_comp_run2(omega_blck_in, pmid_blck_in, zi_blck_in,    &
        hi_blck_in, u_blck_in, v_blck_in, tn_blck_in, sigma_ped_blck_in,      &
        sigma_hall_blck_in, te_blck_in, ti_blck_in, mbar_blck_in,             &
        n2mmr_blck_in, o2mmr_blck_in, o1mmr_blck_in, o2pmmr_blck_in,          &
        nopmmr_blck_in, n2pmmr_blck_in, opmmr_blck_in, opmmrtm1_blck_in,      &
        ui_blck_in, vi_blck_in, wi_blck_in, rmassO2p_in, rmassNOp_in,         &
        rmassN2p_in, rmassOp_in, cols, cole, pver)
      use ESMF,         only: ESMF_GridCompRun
      use edyn_esmf,    only: edyn_esmf_chkerr

      ! Dummy arguments
      real(r8),         pointer   :: omega_blck_in(:,:)
      real(r8),         pointer   :: pmid_blck_in(:,:)
      real(r8),         pointer   :: zi_blck_in(:,:)
      real(r8),         pointer   :: hi_blck_in(:,:)
      real(r8),         pointer   :: u_blck_in(:,:)
      real(r8),         pointer   :: v_blck_in(:,:)
      real(r8),         pointer   :: tn_blck_in(:,:)
      real(r8),         pointer   :: sigma_ped_blck_in(:,:)
      real(r8),         pointer   :: sigma_hall_blck_in(:,:)
      real(r8),         pointer   :: te_blck_in(:,:)
      real(r8),         pointer   :: ti_blck_in(:,:)
      real(r8),         pointer   :: mbar_blck_in(:,:)
      real(r8),         pointer   :: n2mmr_blck_in(:,:)
      real(r8),         pointer   :: o2mmr_blck_in(:,:)
      real(r8),         pointer   :: o1mmr_blck_in(:,:)
      real(r8),         pointer   :: o2pmmr_blck_in(:,:)
      real(r8),         pointer   :: nopmmr_blck_in(:,:)
      real(r8),         pointer   :: n2pmmr_blck_in(:,:)
      real(r8),         pointer   :: opmmr_blck_in(:,:)
      real(r8),         pointer   :: opmmrtm1_blck_in(:,:)
      real(r8),         pointer   :: ui_blck_in(:,:)
      real(r8),         pointer   :: vi_blck_in(:,:)
      real(r8),         pointer   :: wi_blck_in(:,:)
      real(r8)                    :: rmassO2p_in
      real(r8)                    :: rmassNOp_in
      real(r8)                    :: rmassN2p_in
      real(r8)                    :: rmassOp_in
      integer, intent(in)         :: cols
      integer, intent(in)         :: cole
      integer, intent(in)         :: pver
 ! Local variables
      integer                     :: rc
      character(len=*), parameter :: subname = 'edyn_grid_comp_run2'

      do_run = 2
      omega_blck => omega_blck_in
      pmid_blck => pmid_blck_in
      zi_blck => zi_blck_in
      hi_blck => hi_blck_in
      u_blck => u_blck_in
      v_blck => v_blck_in
      tn_blck => tn_blck_in
      sigma_ped_blck => sigma_ped_blck_in
      sigma_hall_blck => sigma_hall_blck_in
      te_blck => te_blck_in
      ti_blck => ti_blck_in
      mbar_blck => mbar_blck_in
      n2mmr_blck => n2mmr_blck_in
      o2mmr_blck => o2mmr_blck_in
      o1mmr_blck => o1mmr_blck_in
      o2pmmr_blck => o2pmmr_blck_in
      nopmmr_blck => nopmmr_blck_in
      n2pmmr_blck => n2pmmr_blck_in
      opmmr_blck => opmmr_blck_in
      opmmrtm1_blck => opmmrtm1_blck_in
      ui_blck => ui_blck_in
      vi_blck => vi_blck_in
      wi_blck => wi_blck_in
      rmassO2p = rmassO2p_in
      rmassNOp = rmassNOp_in
      rmassN2p = rmassN2p_in
      rmassOp = rmassOp_in
      col_start = cols
      col_end = cole
      nlev = pver
      call ESMF_GridCompRun(phys_comp, rc=rc)
      call edyn_esmf_chkerr(subname, 'ESMF_GridCompRun', rc)

   end subroutine edyn_grid_comp_run2

   subroutine edyn_grid_comp_final()
      use ESMF,         only: ESMF_GridCompFinalize
      use edyn_esmf,    only: edyn_esmf_chkerr

      ! Local variables
      integer                     :: rc
      character(len=*), parameter :: subname = 'edyn_grid_comp_final'

      call ESMF_GridCompFinalize(phys_comp, rc=rc)
      call edyn_esmf_chkerr(subname, 'ESMF_GridCompInitialize', rc)

   end subroutine edyn_grid_comp_final


end module edyn_grid_comp
