module gmean_mod
   !-----------------------------------------------------------------------
   !
   ! Purpose:
   ! Perform global mean calculations for energy conservation and other checks.
   !
   ! Method:
   ! Reproducible (scalable):
   !    Convert to fixed point (integer representation) to enable
   !    reproducibility when using MPI collectives.
   ! If error checking is on (via setting reprosum_diffmax > 0 and
   !    reprosum_recompute = .true. in user_nl_cpl), shr_reprosum_calc will
   !    check the accuracy of its computation with a fast but
   !    non-reproducible algorithm. If any error is reported, report
   !    the difference and the expected sum and abort run (call endrun)
   !
   !
   !-----------------------------------------------------------------------
   use shr_kind_mod,     only: r8 => shr_kind_r8
   use ppgrid,           only: pcols, begchunk, endchunk
   use shr_reprosum_mod, only: shr_reprosum_calc, shr_reprosum_tolExceeded
   use shr_reprosum_mod, only: shr_reprosum_reldiffmax, shr_reprosum_recompute
   use perf_mod,         only: t_startf, t_stopf
   use cam_logfile,      only: iulog

   implicit none
   private

   public :: gmean ! compute global mean of 2D fields on physics decomposition
   public :: gmean_init  ! Initialize gmean (maybe run tests)
   public :: test_gmean  ! test accuracy of gmean

   interface gmean
      module procedure gmean_arr
      module procedure gmean_scl
   end interface gmean

   private :: gmean_fixed_repro
   private :: gmean_float_norepro

   ! Set do_gmean_tests to .true. to run a gmean challenge test
   logical, private    :: do_gmean_tests = .false.

CONTAINS

  !
  !========================================================================
  !

   subroutine gmean_init(do_test)
      !-----------------------------------------------------------------------
      !
      ! Purpose: Possibly run a test
      !
      !-----------------------------------------------------------------------
      !
      logical, optional, intent(in) :: do_test

      logical                       :: do_test_use

      if (present(do_test)) then
         do_test_use = do_test
      else
         do_test_use = do_gmean_tests
      end if

      if (do_test_use) then
         call test_gmean()
      end if

   end subroutine gmean_init

   !
   !========================================================================
   !

   subroutine gmean_arr (arr, arr_gmean, nflds)
      use shr_strconvert_mod, only: toString
      use spmd_utils,         only: masterproc
      use cam_abortutils,     only: endrun
      !-----------------------------------------------------------------------
      !
      ! Purpose:
      ! Compute the global mean of each field in "arr" in the physics
      !    chunked decomposition
      !
      ! Method is to call shr_reprosum_calc (called from gmean_fixed_repro)
      !-----------------------------------------------------------------------
      !
      ! Arguments
      !
      integer,  intent(in)  :: nflds            ! number of fields
      real(r8), intent(in)  :: arr(pcols, begchunk:endchunk, nflds)
      real(r8), intent(out) :: arr_gmean(nflds) ! global means
      !
      ! Local workspace
      !
      real(r8)                   :: rel_diff(2, nflds)
      integer                    :: ifld ! field index
      integer                    :: num_err
      logical                    :: write_warning
      !
      !-----------------------------------------------------------------------
      !
      call t_startf('gmean_arr')
      call t_startf ('gmean_fixed_repro')
      call gmean_fixed_repro(arr, arr_gmean, rel_diff, nflds)
      call t_stopf ('gmean_fixed_repro')

      ! check that "fast" reproducible sum is accurate enough. If not, calculate
      ! using old method
      write_warning = masterproc
      num_err = 0
      if (shr_reprosum_tolExceeded('gmean', nflds, write_warning,             &
           iulog, rel_diff)) then
         if (shr_reprosum_recompute) then
            do ifld = 1, nflds
               if (rel_diff(1, ifld) > shr_reprosum_reldiffmax) then
                  call gmean_float_norepro(arr(:,:,ifld), arr_gmean(ifld), ifld)
                  num_err = num_err + 1
               end if
            end do
         end if
      end if
      call t_stopf('gmean_arr')
      if (num_err > 0) then
         call endrun('gmean: '//toString(num_err)//' reprosum errors found')
      end if

   end subroutine gmean_arr

   !
   !========================================================================
   !

   subroutine gmean_scl (arr, gmean)
      use phys_grid, only: get_ncols_p

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
      real(r8), intent(in) :: arr(pcols, begchunk:endchunk)
      ! Input array, chunked
      real(r8), intent(out):: gmean      ! global means
      !
      ! Local workspace
      !
      integer, parameter :: nflds = 1
      real(r8)           :: gmean_array(nflds)
      real(r8)           :: array(pcols, begchunk:endchunk, nflds)
      integer            :: ncols, lchnk

      do lchnk = begchunk, endchunk
         ncols = get_ncols_p(lchnk)
         array(:ncols, lchnk, 1) = arr(:ncols, lchnk)
      end do
      call gmean_arr(array, gmean_array, nflds)
      gmean = gmean_array(1)

   end subroutine gmean_scl

   !
   !========================================================================
   !

   subroutine gmean_float_norepro(arr, repro_sum, index)
      !-----------------------------------------------------------------------
      !
      ! Purpose:
      ! Compute the global mean of <arr> in the physics chunked
      !    decomposition using a fast but non-reproducible algorithm.
      !    Log that value along with the value computed by
      !    shr_reprosum_calc (<repro_sum>)
      !
      !-----------------------------------------------------------------------

      use physconst,  only: pi
      use spmd_utils, only: masterproc, masterprocid, MPI_REAL8, MPI_SUM, mpicom
      use phys_grid,  only: get_ncols_p, get_wght_p
      !
      ! Arguments
      !
      real(r8), intent(in) :: arr(pcols, begchunk:endchunk)
      real(r8), intent(in) :: repro_sum ! Value computed by reprosum
      integer,  intent(in) :: index     ! Index of field in original call
      !
      ! Local workspace
      !
      integer             :: lchnk, ncols, icol
      integer             :: ierr
      real(r8)            :: wght
      real(r8)            :: check
      real(r8)            :: check_sum
      real(r8), parameter :: pi4 = 4.0_r8 * pi

      !
      !-----------------------------------------------------------------------
      !
      ! Calculate and print out non-reproducible value
      check = 0.0_r8
      do lchnk = begchunk, endchunk
         ncols = get_ncols_p(lchnk)
         do icol = 1, ncols
            wght = get_wght_p(lchnk, icol)
            check = check + arr(icol, lchnk) * wght
         end do
      end do
      call MPI_reduce(check, check_sum, 1, MPI_REAL8, check_sum, MPI_SUM,     &
                       masterprocid, mpicom, ierr)
      if (masterproc) then
         write(iulog, '(a,i0,2(a,e20.13e2))') 'gmean(', index, ') = ',        &
              check_sum / pi4, ', reprosum reported ', repro_sum
      end if

   end subroutine gmean_float_norepro

   !
   !========================================================================
   !
   subroutine gmean_fixed_repro (arr, arr_gmean, rel_diff, nflds)
      !-----------------------------------------------------------------------
      !
      ! Purpose:
      ! Compute the global mean of each field in "arr" in the physics
      ! chunked decomposition with a reproducible yet scalable implementation
      ! based on a fixed-point algorithm.
      !
      !-----------------------------------------------------------------------
      use spmd_utils, only: mpicom
      use phys_grid,  only: get_ncols_p, get_wght_all_p, get_nlcols_p
      use phys_grid,  only: ngcols_p => num_global_phys_cols
      use physconst,  only: pi
      !
      ! Arguments
      !
      integer,  intent(in)  :: nflds ! number of fields
      real(r8), intent(in)  :: arr(pcols,begchunk:endchunk,nflds)
      ! arr_gmean: output global sums
      real(r8), intent(out) :: arr_gmean(nflds)
      ! rel_diff: relative and absolute differences from shr_reprosum_calc
      real(r8), intent(out) :: rel_diff(2, nflds)
      !
      ! Local workspace
      !
      integer               :: lchnk, icol, ifld ! chunk, column, field indices
      integer               :: ncols             ! # columns in current chunk
      integer               :: count             ! summand count

      real(r8)              :: wght(pcols)       ! integration weights
      real(r8), allocatable :: xfld(:,:)         ! weighted summands
      integer               :: nlcols
      !
      !-----------------------------------------------------------------------
      !
      nlcols = get_nlcols_p()
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
           gbl_count=ngcols_p, commid=mpicom, rel_diff=rel_diff)

      deallocate(xfld)
      ! final normalization
      arr_gmean(:) = arr_gmean(:) / (4.0_r8 * pi)

   end subroutine gmean_fixed_repro

   subroutine test_gmean(max_diff)
      ! Test gmean on some different field patterns
      ! Test 1: Just 1, easy peasy
      ! Test 2: Positive definite, moderate dynamic range
      ! Test 3: Positive definite, large dynamic range (pattern 1)
      ! Test 4: Positive definite, large dynamic range (pattern 2)
      ! Test 5: Large dynamic range (pattern 1)
      ! Test 6: Large dynamic range (pattern 2)
      use shr_kind_mod,    only: SHR_KIND_CL, INT64 => SHR_KIND_I8
      use physconst,       only: pi
      use spmd_utils,      only: iam, masterproc
      use cam_abortutils,  only: endrun
      use cam_logfile,     only: iulog
      use phys_grid,       only: get_ncols_p, get_gcol_p, get_wght_p
      use phys_grid,       only: ngcols_p => num_global_phys_cols

      ! Dummy argument
      real(r8), optional, intent(in) :: max_diff
      ! Local variables
      integer,          parameter :: num_tests = 6
      integer                     :: lchnk, ncols, icol, gcol, findex
      integer(INT64)              :: test_val
      real(r8)                    :: test_arr(pcols,begchunk:endchunk,num_tests)
      real(r8)                    :: test_mean(num_tests)
      real(r8)                    :: expect(num_tests)
      real(r8)                    :: diff, wght
      real(r8)                    :: max_diff_use
      real(r8),         parameter :: fact2 = 1.0e-8_r8
      real(r8),         parameter :: ifact = 1.0e6_r8
      real(r8),         parameter :: pi4 = 4.0_r8 * pi
      real(r8),         parameter :: max_diff_def = 1.0e-14_r8
      character(len=SHR_KIND_CL)  :: errmsg(num_tests)
      character(len=*), parameter :: subname = 'test_gmean: '

      if (present(max_diff)) then
         max_diff_use = max_diff
      else
         max_diff_use = max_diff_def
      end if
      test_arr = 0.0_r8
      do lchnk = begchunk, endchunk
         ncols = get_ncols_p(lchnk)
         do icol = 1, ncols
            gcol = get_gcol_p(lchnk, icol)
            test_arr(icol, lchnk, 1) = 1.0_r8
            wght = get_wght_p(lchnk, icol)
            test_arr(icol, lchnk, 2) = real(gcol, r8) * pi4 / wght
            test_arr(icol, lchnk, 3) = test_arr(icol, lchnk, 2) * fact2
            if (mod(gcol, 2) == 1) then
               test_arr(icol, lchnk, 4) = test_arr(icol, lchnk, 3) + ifact
               test_arr(icol, lchnk, 6) = test_arr(icol, lchnk, 3) + ifact
            else
               test_arr(icol, lchnk, 4) = test_arr(icol, lchnk, 3)
               test_arr(icol, lchnk, 6) = test_arr(icol, lchnk, 3) - ifact
            end if
            if (gcol > (ngcols_p / 2)) then
               test_arr(icol, lchnk, 5) = test_arr(icol, lchnk, 3) + ifact
               test_arr(icol, lchnk, 3) = test_arr(icol, lchnk, 3) + ifact
            else
               ! test_arr 3 already has correct value
               test_arr(icol, lchnk, 5) = test_arr(icol, lchnk, 3) - ifact
            end if
         end do
      end do
      test_mean(:) = -2.71828_r8 * pi
      expect(1) = 1.0_r8
      test_val = int(ngcols_p, INT64)
      test_val = (test_val + 1) * test_val / 2_INT64
      expect(2) = real(test_val, r8)
      expect(3) = (expect(2) * fact2) + (ifact / 2.0_r8)
      expect(4) = expect(3)
      expect(5) = expect(2) * fact2
      expect(6) = expect(5)
      call gmean(test_arr, test_mean, num_tests)
      errmsg = ''
      do findex = 1, num_tests
         diff = abs(test_mean(findex) - expect(findex)) / expect(findex)
         if (diff > max_diff_use) then
            write(errmsg(findex), '(i0,a,i0,3(a,e20.13e2))') iam,             &
                 ': test_mean(', findex, ') FAIL: ', test_mean(findex),       &
                 ' /= ', expect(findex), ', diff = ', diff
         end if
      end do
      if (ANY(len_trim(errmsg) > 0)) then
         call endrun(subname//trim(errmsg(1))//'\n'//trim(errmsg(2))//'\n'//  &
              trim(errmsg(3))//'\n'//trim(errmsg(4))//'\n'//                  &
              trim(errmsg(5))//'\n'//trim(errmsg(6)))
      end if
      if (masterproc) then
         do findex = 1, num_tests
            write(iulog, '(2a,i0,a,e20.13e2)') subname, 'test_mean(', findex, &
                 ') = ', test_mean(findex)
         end do
      end if
   end subroutine test_gmean

   !
   !========================================================================
   !

end module gmean_mod
