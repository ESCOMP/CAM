module phys_gmean
!-----------------------------------------------------------------------
!
! Purpose:
! Computes global mean mass, max and min mmr, of constituents on the
! physics decomposition. Prints diagnostics to log file.
!
! Author: B. Eaton (based on gavglook)
!
!-----------------------------------------------------------------------
   implicit none
   private

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
      use shr_kind_mod,   only: r8 => shr_kind_r8
      use spmd_utils,     only: MPI_REAL8, MPI_MAX, MPI_MIN
      use spmd_utils,     only: masterproc, masterprocid, mpicom
      use cam_logfile,    only: iulog
      use cam_abortutils, only: endrun
      use ppgrid,         only: pver, pcols, begchunk, endchunk
      use physconst,      only: gravit
      use gmean_mod,      only: gmean
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
      integer                     :: ick, icol, lind, im
      integer                     :: ierr
      integer                     :: ncols

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
      if (ierr /= 0) then
         call endrun(sub_name //'FAIL to allocate mass_wet')
      end if

      allocate(mass_dry(pcols,begchunk:endchunk,pcnst), stat=ierr)
      if (ierr /= 0) then
         call endrun(sub_name//'FAIL to allocate mass_wet')
      end if

      mmr_max(:) = -1.e36_r8
      mmr_min(:) =  1.e36_r8
      do im = 1, pcnst
         do ick = begchunk, endchunk
            ncols = get_ncols_p(ick)
            do icol = 1, ncols

               ! Compute column masses assuming both dry and wet mixing ratios

               mass_wet(icol, ick, im) = 0.0_r8
               do lind = 1, pver
                  mass_wet(icol,ick,im) = mass_wet(icol,ick,im) +     &
                       state(ick)%pdel(icol,lind)*state(ick)%q(icol,lind,im)
                  mmr_max(im) = max(mmr_max(im), state(ick)%q(icol,lind,im))
                  mmr_min(im) = min(mmr_min(im), state(ick)%q(icol,lind,im))
               end do
               mass_wet(icol,ick,im) = mass_wet(icol,ick,im)/gravit

               mass_dry(icol,ick,im) = 0.0_r8
               do lind = 1, pver
                  mass_dry(icol,ick,im) = mass_dry(icol,ick,im) +             &
                       state(ick)%pdeldry(icol,lind)*state(ick)%q(icol,lind,im)
               end do
               mass_dry(icol,ick,im) = mass_dry(icol,ick,im)/gravit

            end do
         end do
      end do

      ! compute global mean mass
      call gmean(mass_wet, mass_wet_mean, pcnst)
      call gmean(mass_dry, mass_dry_mean, pcnst)

      ! global min/max mmr
      call mpi_reduce(mmr_max, mmr_max_glob, pcnst, MPI_REAL8, MPI_MAX,       &
           masterprocid, mpicom, ierr)
      call mpi_reduce(mmr_min, mmr_min_glob, pcnst, MPI_REAL8, MPI_MIN,       &
           masterprocid, mpicom, ierr)

      ! report to log file
      if (masterproc) then
         write(iulog, *) 'vvvvv ', sub_name, trim(title), ' vvvvv'
         write(iulog, *) 'm                  name               ',            &
              '       gavg dry       ', '       gavg wet       ',             &
              '       gavg min       ', '       gavg max       '
         do im = 1, pcnst
            write (iulog, '(i2,a36,4("  ",e20.13e2))') im,                    &
                 trim(cnst_name(im)), mass_dry_mean(im), mass_wet_mean(im),   &
                 mmr_min_glob(im), mmr_max_glob(im)
         end do
         write(iulog, *) '^^^^^ ', sub_name, trim(title), ' ^^^^^'
      end if

      deallocate(mass_wet)
      deallocate(mass_dry)

   end subroutine gmean_mass

end module phys_gmean
