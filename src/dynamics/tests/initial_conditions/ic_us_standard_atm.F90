module ic_us_standard_atmosphere

!-------------------------------------------------------------------------------
!
! Set analytic initial conditions to be static (u=v=0) with temperature profile
! from the US standard atmosphere.
!
!-------------------------------------------------------------------------------

use shr_kind_mod,        only: r8 => shr_kind_r8
use spmd_utils,          only: masterproc

use hycoef,              only: ps0, hyam, hybm
use physconst,           only: gravit
use constituents,        only: cnst_name
use const_init,          only: cnst_init_default

use std_atm_profile,     only: std_atm_pres, std_atm_height, std_atm_temp

use cam_logfile,         only: iulog
use cam_abortutils,      only: endrun

implicit none
private
save

public :: us_std_atm_set_ic

!=========================================================================================
CONTAINS
!=========================================================================================

subroutine us_std_atm_set_ic(latvals, lonvals, U, V, T, PS, PHIS_IN, &
       PHIS_OUT, Q, m_cnst, mask, verbose)
    
   !----------------------------------------------------------------------------
   !
   ! Set initial values for static atmosphere with vertical profile from US
   ! Standard Atmosphere.  
   !
   !----------------------------------------------------------------------------

   ! Arguments
   real(r8),           intent(in)    :: latvals(:) ! lat in degrees (ncol)
   real(r8),           intent(in)    :: lonvals(:) ! lon in degrees (ncol)
   real(r8), optional, intent(inout) :: U(:,:)     ! zonal velocity
   real(r8), optional, intent(inout) :: V(:,:)     ! meridional velocity
   real(r8), optional, intent(inout) :: T(:,:)     ! temperature
   real(r8), optional, intent(inout) :: PS(:)      ! surface pressure
   real(r8), optional, intent(in)    :: PHIS_IN(:) ! surface geopotential
   real(r8), optional, intent(out)   :: PHIS_OUT(:)! surface geopotential
   real(r8), optional, intent(inout) :: Q(:,:,:)   ! tracer (ncol, lev, m)
   integer,  optional, intent(in)    :: m_cnst(:)  ! tracer indices (reqd. if Q)
   logical,  optional, intent(in)    :: mask(:)    ! Only init where .true.
   logical,  optional, intent(in)    :: verbose    ! For internal use

   ! Local variables
   logical, allocatable              :: mask_use(:)
   logical                           :: verbose_use
   integer                           :: i, k, m
   integer                           :: ncol
   integer                           :: nlev
   integer                           :: ncnst
   character(len=*), parameter       :: subname = 'us_std_atm_set_ic'
   real(r8)                          :: psurf(1)
   real(r8), allocatable             :: pmid(:), zmid(:)
   !----------------------------------------------------------------------------

   ncol = size(latvals, 1)
   allocate(mask_use(ncol))
   if (present(mask)) then
      if (size(mask_use) /= size(mask)) then
         call endrun(subname//': input, mask, is wrong size')
      end if
      mask_use = mask
   else
      mask_use = .true.
   end if

   if (present(verbose)) then
      verbose_use = verbose
   else
      verbose_use = .true.
   end if

   if (present(PHIS_OUT)) then
      PHIS_OUT = 0.0_r8
      if (masterproc .and. verbose_use) then
         write(iulog,*) '          PHIS initialized by '//subname
      end if
   end if

   nlev = -1
   if (present(U)) then
      nlev = size(U, 2)
      do k = 1, nlev
         where(mask_use)
            U(:,k) = 0.0_r8
         end where
      end do
      if(masterproc .and. verbose_use) then
         write(iulog,*) '          U initialized by '//subname
      end if
   end if

   if (present(V)) then
      nlev = size(V, 2)
      do k = 1, nlev
         where(mask_use)
            V(:,k) = 0.0_r8
         end where
      end do
      if(masterproc .and. verbose_use) then
         write(iulog,*) '          V initialized by '//subname
      end if
   end if

   if (present(T)) then
      if (.not.present(PHIS_IN)) then
         call endrun(subname//': PHIS_IN must be specified to initiallize T')
      end if
      nlev = size(T, 2)
      allocate(pmid(nlev), zmid(nlev))
      do i = 1, ncol
         if (mask_use(i)) then
            ! get surface pressure
            call std_atm_pres(PHIS_IN(i:i)/gravit, psurf)
            ! get pressure levels
            do k = 1, nlev
               pmid(k) = hyam(k)*ps0 + hybm(k)*psurf(1)
            end do
            ! get height of pressure level            
            call std_atm_height(pmid, zmid)
            ! given height get temperature
            call std_atm_temp(zmid, T(i,:))
         end if
      end do
      deallocate(pmid, zmid)

      if(masterproc .and. verbose_use) then
         write(iulog,*) '          T initialized by "',subname,'"'
      end if
   end if

   if (present(PS)) then
      if (.not.present(PHIS_IN)) then
         call endrun(subname//': PHIS_IN must be specified to initiallize PS')
      end if
      
      do i = 1, ncol
         if (mask_use(i)) then
            call std_atm_pres(PHIS_IN(i:i)/gravit, PS(i:i))
         end if
      end do
      if(masterproc .and. verbose_use) then
         write(iulog,*) '          PS initialized by "',subname,'"'
      end if
   end if

   if (present(Q)) then
      nlev = size(Q, 2)
      ncnst = size(m_cnst, 1)
      do m = 1, ncnst
         if (m_cnst(m) == 1) then
            ! No water vapor in profile
            do k = 1, nlev
               where(mask_use)
                  Q(:,k,m_cnst(m)) = 0.0_r8
               end where
            end do
            if(masterproc .and. verbose_use) then
               write(iulog,*) '          ', trim(cnst_name(m_cnst(m))), ' initialized by '//subname
            end if
         else
            call cnst_init_default(m_cnst(m), latvals, lonvals, Q(:,:,m_cnst(m)),&
               mask=mask_use, verbose=verbose_use, notfound=.false.)
         end if
      end do
   end if

   deallocate(mask_use)

end subroutine us_std_atm_set_ic

!=========================================================================================

end module ic_us_standard_atmosphere
