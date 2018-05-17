
module tracers_suite

!---------------------------------------------------------------
!
! Implements artificial suite of passive tracers
!    1) low tracer with unit mixing ratio at about 800 mb
!    2) med tracer with unit mixing ratio at about 500 mb
!    3) high tracer with unit mixing ratio at about 200 mb
!    4) reverse med tracer with unit mixing ratio everywhere except about 500 mb
!    5) unit tracer with unit mixing ratio everywhere
!
!  D Bundy June 2003
!  modified Feb 2004 to include TT_UN and smoothing
!
!  A Mirin and B Eaton, August 2007
!  Modified to create up to 1000 distinct copies of the 5 basic tracers
!  by appending up to a 3 digit number to the base tracer name.
!  RESTRICTION - trac_ncnstmx cannot exceed 5000 unless the algorithm for
!                constructing new tracer names is extended.
!
!---------------------------------------------------------------

use shr_kind_mod,   only: r8 => shr_kind_r8
use ppgrid,         only: pcols, pver
use ref_pres,       only: pref_mid
use cam_abortutils, only: endrun
use cam_logfile,    only: iulog

implicit none
private
save

public get_tracer_name  ! generate names of tracers
public init_cnst_tr     ! initialize tracer fields

integer, parameter :: trac_ncnstmx=5000 ! Max no. of tracers based on current algorithm for
                                        ! constructing tracer names.  This could easily be extended.
integer, parameter :: trac_names=5      ! No. of base tracers

logical, parameter :: smooth = .false.
  
!======================================================================
contains
!======================================================================

function get_tracer_name(n)

   ! The tracer names are only defined in this module. This function is for
   ! outside programs to grab the name for each tracer number. 

   integer, intent(in) :: n
   character(len=8)    :: get_tracer_name

   ! Local variables
   character(len=5), dimension(trac_names), parameter :: & ! constituent names
      tracer_names  =  (/ 'TT_LW', 'TT_MD', 'TT_HI', 'TTRMD' , 'TT_UN'/)

   integer :: nbase  ! Corresponding base tracer index
   integer :: ncopy  ! No. of copies of base tracers
   character(len=1) :: c1
   character(len=2) :: c2
   character(len=3) :: c3
   !-----------------------------------------------------------------------

   if ( n > trac_ncnstmx ) then
      write(iulog,*) 'tracers_suite:get_tracer_name()','requested tracer',n
      write(iulog,*) 'only ',trac_ncnstmx,' tracers available'
      call endrun('tracers_suite: ERROR in get_tracer_name(); n too large')
   else
      nbase = mod(n-1, trac_names) + 1
      ncopy = (n-1)/trac_names
      if ( ncopy == 0 ) then
         get_tracer_name = tracer_names(nbase)
      else if ( ncopy >= 1  .and.  ncopy <= 9 ) then
         write (c1,'(i1)') ncopy
         get_tracer_name = tracer_names(nbase) // c1
      else if ( ncopy >= 10  .and.  ncopy <= 99 ) then
         write (c2,'(i2)') ncopy
         get_tracer_name = tracer_names(nbase) // c2
      else if ( ncopy >= 100  .and.  ncopy <= 999 ) then
         write (c3,'(i3)') ncopy
         get_tracer_name = tracer_names(nbase) // c3
      end if
   endif

end function get_tracer_name

!======================================================================

subroutine init_cnst_tr(m, latvals, lonvals, mask, q)

   ! calls initialization routine for tracer m, returns mixing ratio in q

   integer,  intent(in)  :: m          ! index of tracer
   real(r8), intent(in)  :: latvals(:) ! lat in radians (ncol)
   real(r8), intent(in)  :: lonvals(:) ! lon in radians (ncol)
   logical,  intent(in)  :: mask(:)    ! Only initialize where .true.
   real(r8), intent(out) :: q(:,:)     ! kg tracer/kg dry air (gcol,plev)

   integer nbase ! Corresponding base tracer index

   if ( m > trac_ncnstmx ) then
      write(iulog,*) 'tracers_suite:init_cnst_tr()'
      write(iulog,*) ' asked to initialize tracer number ',m
      write(iulog,*) ' but there are only trac_ncnstmx = ',trac_ncnstmx,' tracers'
      call endrun('tracers_suite: ERROR in init_cnst_tr(); m too large')
   endif

   nbase = mod(m-1,trac_names)+1

   if ( nbase == 1 ) then
      call init_cnst_lw(latvals, lonvals, mask, q)
   else if ( nbase == 2 ) then
      call init_cnst_md(latvals, lonvals, mask, q)
   else if ( nbase == 3 ) then
      call init_cnst_hi(latvals, lonvals, mask, q)
   else if ( nbase == 4 ) then
      call init_cnst_md(latvals, lonvals, mask, q, rev_in=1)
   else if ( nbase == 5 ) then
      call init_cnst_un(latvals, lonvals, mask, q)
   else
      write(iulog,*) 'tracers_suite:init_cnst_tr()'
      write(iulog,*) 'no initialization routine specified for tracer',nbase
      call endrun('tracers_suite: ERROR in init_cnst_tr(); no init routine available')
   endif
      
end subroutine init_cnst_tr

!======================================================================

subroutine init_cnst_lw(latvals, lonvals, mask, q)

   ! Initialize test tracer TT_LW
   ! Initialize low tracer to zero except at 800 level

   ! Arguments
   real(r8), intent(in)  :: latvals(:) ! lat in radians (ncol)
   real(r8), intent(in)  :: lonvals(:) ! lon in radians (ncol)
   logical,  intent(in)  :: mask(:)    ! Only initialize where .true.
   real(r8), intent(out) :: q(:,:)    ! kg tracer/kg dry air (gcol,plev)

   ! Local
   integer :: indx, k
   !-----------------------------------------------------------------------

   indx = setpindxtr(800._r8)

   if ( smooth ) then 
      call setsmoothtr(indx,q,.876_r8, mask)
   else 
      do k = 1, size(q, 2)
         if (k == indx) then
            where(mask)
               q(:,k) = 1.0_r8
            end where
         else
            where(mask)
               q(:,k) = 0.0_r8
            end where
         end if
      end do
   end if

end subroutine init_cnst_lw

!======================================================================

subroutine init_cnst_md(latvals, lonvals, mask, q, rev_in)

   ! Initialize test tracer TT_MD
   ! Initialize med tracer to zero except at 500 level

   ! Arguments
   real(r8), intent(in)  :: latvals(:) ! lat in radians (ncol)
   real(r8), intent(in)  :: lonvals(:) ! lon in radians (ncol)
   logical,  intent(in)  :: mask(:)    ! Only initialize where .true.
   real(r8), intent(out) :: q(:,:)    ! kg tracer/kg dry air
   integer,  intent(in), optional :: rev_in         ! reverse the mixing ratio

   ! Local
   integer :: indx, k
   integer :: rev
   !-----------------------------------------------------------------------

   rev = 0
   if (present(rev_in)) then
      if (rev_in == 1) then
         rev = 1
      endif
   endif

   indx = setpindxtr(500._r8)

   if ( smooth ) then 
      call setsmoothtr(indx,q,.876_r8,mask,rev_in=rev)
   else
      do k = 1, size(q, 2)
         if (rev == 1 ) then
            if (k == indx) then
               where(mask)
                  q(:,indx) = 0.0_r8
               end where
            else
               where(mask)
                  q(:,k) = 1.0_r8
               end where
            end if
         else
            if (k == indx) then
               where(mask)
                  q(:,indx) = 1.0_r8
               end where
            else
               where(mask)
                  q(:,k) = 0.0_r8
               end where
            end if
         end if
      end do
   end if
   
end subroutine init_cnst_md

!======================================================================

subroutine init_cnst_hi(latvals, lonvals, mask, q)

   ! Initialize test tracer TT_HI
   ! Initialize high tracer to zero except at 200 level

   ! Arguments
   real(r8), intent(in)  :: latvals(:) ! lat in radians (ncol)
   real(r8), intent(in)  :: lonvals(:) ! lon in radians (ncol)
   logical,  intent(in)  :: mask(:)    ! Only initialize where .true.
   real(r8), intent(out) :: q(:,:)    ! kg tracer/kg dry air

   ! Local
   integer :: indx, k
   !-----------------------------------------------------------------------

   indx = setpindxtr(200._r8)

   if ( smooth ) then 
      call setsmoothtr(indx,q,.3_r8,mask)
   else
      do k = 1, size(q, 2)
         if (k == indx) then
            where(mask)
               q(:,k) = 1.0_r8
            end where
         else
            where(mask)
               q(:,k) = 0.0_r8
            end where
         end if
      end do
   end if

end subroutine init_cnst_hi

!======================================================================

subroutine init_cnst_un(latvals, lonvals, mask, q)

   ! Initialize test unit tracer TT_UN

   real(r8), intent(in)  :: latvals(:) ! lat in radians (ncol)
   real(r8), intent(in)  :: lonvals(:) ! lon in radians (ncol)
   logical,  intent(in)  :: mask(:)    ! Only initialize where .true.
   real(r8), intent(out) :: q(:,:)    ! kg tracer/kg dry air
   !-----------------------------------------------------------------------
   integer :: k

   do k = 1, size(q, 2)
      where(mask)
         q(:,k) = 1.0_r8
      end where
   end do

end subroutine init_cnst_un

!======================================================================

subroutine setsmoothtr(indx,q,width,mask,rev_in)

   ! Arguments
   integer, intent(in)     :: indx               ! k index of pressure level
   real(r8), intent(inout) :: q(:,:)  ! kg tracer/kg dry air
   real(r8), intent(in)    :: width              ! eta difference from unit level where q = 0.1
   logical,  intent(in)    :: mask(:) ! Only set q where mask is .true.
   integer,  intent(in), optional :: rev_in      ! reverse the mixing ratio

   ! Local variables
   integer  :: k
   real(r8) :: alpha ! guassian width, determined by width, T
   real(r8) :: pdist ! pressure distance (eta.e4) from k=indx
   real(r8) :: T     ! desired m.r. in level specified by pdiff from k=indx
   integer  :: rev  ! = 1 then reverse (q = 1, q(k=indx) = 0 )

   rev = 0
   if (present(rev_in)) then
      if (rev_in == 1) then
         rev = 1
      endif
   endif

   T = 0.1_r8
   alpha = -log(T)/(width*1.e4_r8)**2  ! s.t. in level width from indx, mr = T

   !  alpha = 3.e-8  ! m.r. ~ 0.1 in adjacent levels, where change eta ~ 0.08

   do k=1,pver
      pdist = pref_mid(k) - pref_mid(indx)

      if ( rev == 1 ) then
         where(mask)
            q(:,k) = 1.0_r8 - exp(-alpha*(pdist**2))
         end where
      else
         where(mask)
            q(:,k) =       exp(-alpha*(pdist**2))
         end where
      endif
   end do
  
end subroutine setsmoothtr

!======================================================================

integer function setpindxtr(pmb)

   ! find the index of layer nearest pmb

   real(r8), intent(in) :: pmb

   integer  :: indx, k
   real(r8) :: pmin, pdist
  
   indx = 0
   pmin = 1.e36_r8
   pdist = 1.e36_r8
   do k=1,pver
      pdist = abs(pref_mid(k) - pmb*100._r8)
      if (pdist < pmin) then
         indx = k
         pmin = pdist
      end if
   end do

   setpindxtr = indx

end function setpindxtr

!======================================================================

end module tracers_suite
