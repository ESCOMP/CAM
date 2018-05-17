module phys_gmean
!-----------------------------------------------------------------------
!
! Purpose:
! Perform mixed layer global calculations for energy conservation checks.
!
! Methods:
! Reproducible (nonscalable):
!    Gather to a master processor who does all the work.
! Reproducible (scalable):
!    Convert to fixed point (integer representation) to enable
!    reproducibility when using MPI collectives. Results compared with
!    a nonreproducible (but scalable) algorithm using floating point
!    and MPI_Allreduce to verify the results are good enough.
!
! Author: Byron Boville from SOM code by Jim Rosinski/Bruce Briegleb
! Modified: P. Worley to aggregate calculations (4/04)
! Modified: J. White/P. Worley to introduce scalable algorithms;
!           B. Eaton to remove dycore-specific dependencies and to
!           introduce gmean_mass (10/07)
! Modified: P. Worley to replace in-place implementation with call
!           to repro_sum.
!
!-----------------------------------------------------------------------
   use shr_kind_mod,  only: r8 => shr_kind_r8
   use physconst,     only: pi
   use spmd_utils,    only: masterproc, MPI_REAL8, MPI_MAX, MPI_MIN, mpicom
   use gmean_mod,     only: gmean
   use ppgrid,        only: pcols, begchunk, endchunk
   use shr_reprosum_mod, only: shr_reprosum_calc, shr_reprosum_tolExceeded, &
                            shr_reprosum_reldiffmax, shr_reprosum_recompute
   use perf_mod
   use cam_logfile,   only: iulog

   implicit none
   private
   save

   public :: gmean_mass ! compute global mean mass of constituent fields on physics decomposition

   CONTAINS

!
!========================================================================
!

   subroutine gmean_mass(title, state)
!-----------------------------------------------------------------------
!
! Purpose:
! Computes global mean mass, max and min mmr, of constituents on the
! physics decomposition. Prints diagnostics to log file.
!
! Author: B. Eaton (based on gavglook)
!
!-----------------------------------------------------------------------
      use ppgrid,         only: pver
      use physconst,      only: gravit
      use phys_grid,      only: get_ncols_p
      use physics_types,  only: physics_state
      use constituents,   only: pcnst, cnst_name
!
! Arguments
!
      character(len=*),    intent(in) :: title    ! location of this call
      type(physics_state), intent(in) :: state(begchunk:endchunk)
!
! Local workspace
!
      character(len=*), parameter :: sub_name='gmean_mass: '

      integer :: c, i, k, m
      integer :: ierr
      integer :: ncols

      real(r8), pointer :: mass_wet(:,:,:) ! constituent masses assuming moist mmr
      real(r8), pointer :: mass_dry(:,:,:) ! constituent masses assuming dry mmr
      real(r8) :: mass_wet_mean(pcnst)     ! global mean constituent masses assuming moist mmr
      real(r8) :: mass_dry_mean(pcnst)     ! global mean constituent masses assuming dry mmr
      real(r8) :: mmr_max(pcnst)           ! maximum constituent mmr in this process
      real(r8) :: mmr_min(pcnst)           ! minimum constituent mmr in this process
      real(r8) :: mmr_max_glob(pcnst)      ! global maximum constituent mmr
      real(r8) :: mmr_min_glob(pcnst)      ! global minimum constituent mmr
!
!-----------------------------------------------------------------------
!
      allocate(mass_wet(pcols,begchunk:endchunk,pcnst), stat=ierr)
      if (ierr /= 0) write(iulog,*) sub_name // 'FAIL to allocate mass_wet'

      allocate(mass_dry(pcols,begchunk:endchunk,pcnst), stat=ierr)
      if (ierr /= 0) write(iulog,*) sub_name // 'FAIL to allocate mass_wet'

      mmr_max(:) = -1.e36_r8
      mmr_min(:) =  1.e36_r8
      do m = 1, pcnst
         do c = begchunk, endchunk
            ncols = get_ncols_p(c)
            do i = 1, ncols

               ! Compute column masses assuming both dry and wet mixing ratios

               mass_wet(i,c,m) = 0.0_r8
               do k = 1, pver
                  mass_wet(i,c,m) = mass_wet(i,c,m) + &
                                    state(c)%pdel(i,k)*state(c)%q(i,k,m)
                  mmr_max(m) = max(mmr_max(m), state(c)%q(i,k,m))
                  mmr_min(m) = min(mmr_min(m), state(c)%q(i,k,m))
               end do
               mass_wet(i,c,m) = mass_wet(i,c,m)/gravit

               mass_dry(i,c,m) = 0.0_r8
               do k = 1, pver
                  mass_dry(i,c,m) = mass_dry(i,c,m) + &
                                    state(c)%pdeldry(i,k)*state(c)%q(i,k,m)
               end do
               mass_dry(i,c,m) = mass_dry(i,c,m)/gravit

            end do
         end do
      end do

      ! compute global mean mass
      call gmean(mass_wet, mass_wet_mean, pcnst)
      call gmean(mass_dry, mass_dry_mean, pcnst)

      ! global min/max mmr
      call mpi_reduce(mmr_max, mmr_max_glob, pcnst, MPI_REAL8, MPI_MAX, 0, mpicom, ierr)
      call mpi_reduce(mmr_min, mmr_min_glob, pcnst, MPI_REAL8, MPI_MIN, 0, mpicom, ierr)

      ! report to log file
      if (masterproc) then

         do m = 1, pcnst
               write (6,66) trim(title)//' m=',m, &
                  'name='//trim(cnst_name(m))//' gavg dry, wet, min, max ', &
                  mass_dry_mean(m), mass_wet_mean(m), mmr_min_glob(m), mmr_max_glob(m)
66             format (a24,i2,a36,1p,4e25.13)
         end do

      endif

      deallocate(mass_wet)
      deallocate(mass_dry)

   end subroutine gmean_mass

end module phys_gmean
