 
      module exbdrift
!---------------------------------------------------------------------- 
! description: calculates ExB drift velocities UI,VI,WI
!     uses the electric field which is calculated in module efield
!     on a regular magnetic grid (MLT deg/ mag. latitude)
!
!     0. initilize called before time-loop (exbdrift_init)
!     every timestep and every processor
!     1. map from magn. grid to geographic grid WACCM (map_mag2geo)
!     2. rotate e-field (rot_efield)
!     3. calculate ExB drift velocities  ui,vi,wi [m/s] (iondrift)
!
! input ed1,ed2,    &    ! zonal/meridional elect. field [V/m]
!      nmlon,nmlat, &    ! dimension of mag. grid 
!      dlatm,dlonm, &    ! grid spacing of mag. grid 
!      ylatm,ylonm       ! magnetic MLT deg./latitudes (deg)
!
! ExB electromagnetic drift velocity [m/s] (east,north,upward)
! These are output to the physics buffer.
!      real(r8) ::  ui,vi,wi
!
! notes: 
! - assume regular magnetic grid for the e-field (0:360 MLT/0:90) -> mapping
!
! Author: A. Maute Dec 2003
!         B. Foster adding physics buffer Feb, 2004.
!---------------------------------------------------------------------- 

      use shr_kind_mod,  only: r8 => shr_kind_r8
      use physconst,     only: pi
      
      use ppgrid      ,  only: pcols, pver
      use cam_logfile,   only: iulog

      use efield, only : &  ! inputs from efield module
        ed1, ed2,     &     ! global zonal/meridional elect. field [V/m]
        potent,       &     ! electric potential [V]
        nmlon, nmlat, &     ! dimension of mag. grid 
        dlatm, dlonm, &     ! grid spacing of mag. grid 
        ylatm, ylonm        ! magnetic longitudes,latitudes (deg) (0:nmlat),(0:nmlon)
      use apex, only : apex_subsol, apex_magloctm
      use cam_history,  only: outfld, addfld, add_default, horiz_only ! for history saves
      use shr_assert_mod, only: shr_assert_in_domain
      use phys_control, only: phys_getopts

      implicit none

      private

      save

!---------------------------------------------------------------------- 
! Public interfaces:
!---------------------------------------------------------------------- 
      public :: exbdrift_init
      public :: exbdrift_register       ! register drift velocities with pbuf 
      public :: exbdrift_ion_vels       ! updates empirical ion drift velocities (in pbuf)

!---------------------------------------------------------------------- 
! Indices to drift velocities in physics buffer:
!---------------------------------------------------------------------- 
      integer  :: ndx_ui, ndx_vi, ndx_wi
      real(r8), parameter :: rtd = 180._r8/pi ! radians to degrees
      real(r8), parameter :: hr2d = 360._r8/24._r8

      logical :: state_debug_checks = .false.

      contains

      subroutine exbdrift_init( empirical_ion_vels )
!-----------------------------------------------------------------------
! Purpose: Prepare fields for histories.
!
! Method: 
!
! Author: A. Maute Dec 2003  am 12/30/03    
!-----------------------------------------------------------------------

      use phys_control, only: phys_getopts

      logical, intent(in) :: empirical_ion_vels

      logical :: history_waccm


      call phys_getopts(history_waccm_out=history_waccm, state_debug_checks_out=state_debug_checks)

!-----------------------------------------------------------------------
! Add mag field output to master field list:
!-----------------------------------------------------------------------
      call addfld('EF_EAST', horiz_only,'I','V/m', 'eastward electric field')
      call addfld('EF_NORTH', horiz_only,'I','V/m', 'northward electric field')
      call addfld('EF_UP',   horiz_only,'I','V/m', 'upward electric field')
      call addfld('EF1_MAP', horiz_only,'I','V/m', 'map. mag. eastward ef')
      call addfld('EF2_MAP', horiz_only,'I','V/m', 'map. mag. northward ef')
      call addfld('EPOTEN',  horiz_only,'I','V', 'Electric Potential')
!-----------------------------------------------------------------------
!  Write these fields to WACCM history by default:
!-----------------------------------------------------------------------
      if (history_waccm .and. empirical_ion_vels) then
         call add_default ('EPOTEN ' , 1, ' ')
      end if

      end subroutine exbdrift_init

      subroutine exbdrift_register

      use physics_buffer, only : pbuf_add_field, dtype_r8

!-----------------------------------------------------------------------
! Register drift velocity outputs with physics buffer:
!
! Ion velocities are 2d fields in pbuf (no vertical dimension),
!   so fdim = mdim = ldim = 1.
! Indices are saved as module data above.
!-----------------------------------------------------------------------

      call pbuf_add_field("UI", "global", dtype_r8, (/pcols,pver/), ndx_ui)
      call pbuf_add_field("VI", "global", dtype_r8, (/pcols,pver/), ndx_vi)
      call pbuf_add_field("WI", "global", dtype_r8, (/pcols,pver/), ndx_wi)

      end subroutine exbdrift_register

      subroutine exbdrift_ion_vels( lchnk, ncol, pbuf)
      use physics_buffer, only : physics_buffer_desc
!-----------------------------------------------------------------------
! Purpose: calculate ion drift velocities  [m/s]
!
! Method: v = E x B/B^2
!
! Author: A. Maute Dec 2003  am 12/30/03    
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!	... dummy arguments
!-----------------------------------------------------------------------
      integer, intent(in) :: ncol    ! np. of atmospheric columns
      integer, intent(in) :: lchnk   ! chunk identifier
      
      type(physics_buffer_desc), pointer :: pbuf(:)

!-----------------------------------------------------------------------
!	... local variables
!-----------------------------------------------------------------------
      real(r8) :: ed1_geo(ncol), & ! electric field on geographic grid [V/m]
                  ed2_geo(ncol), &
                  epot_geo(ncol)   ! electric potential on geographic grid
      real(r8) :: elfld(3,ncol)    ! electric field in geog. direction on geographic grid [V/m]
      real(r8) :: mlt(ncol)        ! mag.local time of WACCM geo. grid point

!-----------------------------------------------------------------------
! calculate the magn. local time of WACCM geogr. grid points
!-----------------------------------------------------------------------
      call cal_mlt( mlt, lchnk, ncol )
!-----------------------------------------------------------------------
! map the electric field from regular mag.grid to WACCM grid 
!-----------------------------------------------------------------------
      call map_mag2geo( mlt, lchnk, ncol, ed1_geo, ed2_geo, epot_geo )
!-----------------------------------------------------------------------
! rotate the electric field in geographic direction 
!-----------------------------------------------------------------------
      call rot_efield( lchnk, ncol, ed1_geo, ed2_geo, elfld )
!-----------------------------------------------------------------------
! calculate ExB drift velocities  ui,vi,wi [m/s]
!-----------------------------------------------------------------------
      call iondrift( lchnk, ncol, elfld, pbuf)

      end subroutine exbdrift_ion_vels

      subroutine map_mag2geo( mlt, lchnk, pcol, ed1_geo, ed2_geo, epot_geo )
!-----------------------------------------------------------------------
! Purpose: map electric field from regular magnetic grid 
!    to WACCM physics geographic grid
!
! Method:  bilinear interpolation
!    assumptions: magnetic grid regular (MLT[deg]=0:360 deg/lat=0:180 deg)
!    ed1_geo(:,:) = Ed1_geo <- mapping of - 1/[R cos lam_m] d PHI/d phi_m
!    ed2_geo(:,:) = Ed2_geo <- mapping of 1/R d PHI/d lam_m/ sinIm
!
! Author: A.Maute Dec 2003 am 12/17/03 
!-----------------------------------------------------------------------

   use mo_apex,     only:  alatm                ! apex mag latitude at each geographic grid point (radians)

!-----------------------------------------------------------------------
! dummy arguments
!-----------------------------------------------------------------------
   integer, intent(in)   :: pcol      ! np. of atmospheric columns
   integer, intent(in)   :: lchnk     ! chunk identifier
   real(r8), intent(in)  :: mlt(pcol) ! mag.local time of WACCM geo. grid point
   real(r8), intent(out) :: &
     ed1_geo(pcol), &                 ! electric field on geog. grid
     ed2_geo(pcol), &                 ! electric field on geog. grid
     epot_geo(pcol)                   ! electric potential on geog grid

!-----------------------------------------------------------------------
! local variables
!-----------------------------------------------------------------------
     integer  :: i, iphi1, iphi2, ilam1, ilam2
     real(r8) :: t, u, collat, collon

     if (state_debug_checks) then
        call shr_assert_in_domain(ed1, is_nan=.false., varname="ed1", msg="NaN found in exbdrift::map_mag2geo")
        call shr_assert_in_domain(ed2, is_nan=.false., varname="ed2", msg="NaN found in exbdrift::map_mag2geo")
        call shr_assert_in_domain(potent, is_nan=.false., varname="potent", msg="NaN found in exbdrift::map_mag2geo")
     end if

     do i = 1,pcol
       collat = alatm(i,lchnk)*rtd      ! mag lats (deg) [-90:90]
       collon = mlt(i)*hr2d             ! mlt (deg) [0:360]
       iphi1 = int(collon/dlonm)	! indices of the 4 surrounding points (regular mag. grid
       iphi2 = iphi1 + 1	        ! with dlonm/dlatm deg grid spacing)
       ilam1 = int( (collat + 90._r8)/dlatm )
       ilam2 = ilam1 + 1
       if(iphi1 == nmlon ) then  ! boundaries 
           iphi1 = nmlon -1
           iphi2 = nmlon
       end if
       if(ilam1 == nmlat ) then 
           ilam1 = nmlat-1
           ilam2 = nmlat
       end if
       if(collon == 360._r8) then
          collon= 0._r8
          iphi1 = 0
          iphi2 = 1
       end if

       t = (collon - ylonm(iphi1))/(ylonm(iphi2) - ylonm(iphi1))
       u = (collat + 90._r8 - ylatm(ilam1))/(ylatm(ilam2) - ylatm(ilam1))
       ed1_geo(i) = (1._r8 - t)*(1._r8 - u)*ed1(iphi1,ilam1) + &
                    t*(1._r8 - u)*ed1(iphi2,ilam1) + &
                    t*u*          ed1(iphi2,ilam2) + &
                    (1._r8 - t)*u*ed1(iphi1,ilam2) 
       ed2_geo(i) = (1._r8 - t)*(1._r8 - u)*ed2(iphi1,ilam1) + &
                    t*(1._r8 - u)*ed2(iphi2,ilam1) + &
                    t*u*          ed2(iphi2,ilam2) + &
                    (1._r8 - t)*u*ed2(iphi1,ilam2) 
       epot_geo(i)= (1._r8 - t)*(1._r8 - u)*potent(iphi1,ilam1) + &
                    t*(1._r8 - u)*potent(iphi2,ilam1) + &
                    t*u*          potent(iphi2,ilam2) + &
                    (1._r8 - t)*u*potent(iphi1,ilam2)
       
     end do ! i = 1,pcol

     if (state_debug_checks) then
        call shr_assert_in_domain(ed1_geo, is_nan=.false., varname="ed1_geo", msg="NaN found in exbdrift::map_mag2geo")
        call shr_assert_in_domain(ed2_geo, is_nan=.false., varname="ed2_geo", msg="NaN found in exbdrift::map_mag2geo")
        call shr_assert_in_domain(epot_geo, is_nan=.false., varname="epot__geo", msg="NaN found in exbdrift::map_mag2geo")
     endif

     call outfld( 'EF1_MAP', ed1_geo, pcol, lchnk)
     call outfld( 'EF2_MAP', ed2_geo, pcol, lchnk)
     call outfld( 'EPOTEN', epot_geo, pcol, lchnk)

   end subroutine map_mag2geo

   subroutine rot_efield( lchnk, pcol, ed1_geo, ed2_geo, elfld )
!-----------------------------------------------------------------------
! Purpose: rotate the electric field to get geographic east and westward direction
!
! Method: Richmond: J Geomag. Geoelectr. [1995] eqn. (4.5)
!      rotation of the electric field to get geog. eastward and downward/equatorward
!      E(k) = d_1(k)*ed1_geo + d_2(k)*ed2_geo  for k = 1,3
!
! Author: A. Maute Dec 2003  am 12/17/03 
!-----------------------------------------------------------------------

   use mo_apex,     only: d1vec, d2vec ! base vectors (3,pcols,begchunk:endchunk)      

!-----------------------------------------------------------------------
! dummy arguments
!-----------------------------------------------------------------------
   integer, intent(in)   :: pcol                         ! np. of atmospheric columns
   integer, intent(in)   :: lchnk                        ! chunk identifier
   real(r8), intent(in ) :: ed1_geo(pcol),ed2_geo(pcol)  ! electric field on geog. grid
   real(r8), intent(out) :: elfld(3,pcol)                ! electric field on geog. grid geog. direction

!-----------------------------------------------------------------------
! local
!-----------------------------------------------------------------------
   integer :: i, k

   do i = 1,pcol
     do k = 1,3
       elfld(k,i) = ed1_geo(i)*d1vec(k,i,lchnk) + & 
                    ed2_geo(i)*d2vec(k,i,lchnk)
     end do
   end do

   call outfld( 'EF_EAST', elfld(1,:), pcol, lchnk)
   call outfld( 'EF_NORTH', elfld(2,:), pcol, lchnk)
   call outfld( 'EF_UP', elfld(3,:), pcol, lchnk)

   end subroutine rot_efield

   subroutine iondrift( lchnk, pcol, elfld, pbuf)
!-----------------------------------------------------------------------
! Purpose: calculate ion drift velocity
!
! Method: v_i = ExB/B^2
!      for high altitudes where collisionfrequency mue_in << omega_i and
!      mue_en,vert << omega_e
!      v_i,vert = v_e,vert = v_E
!      B magnetic field component from apex code
!         bnorth northward gauss
!         beast  eastward  gauss
!         bdown  downward  gauss -> for upward component -bdown
!         bmag   magnitude gauss
!
! Author: A. Maute Dec 2003  am 12/17/03 
!-----------------------------------------------------------------------

   use mo_apex,        only: beast, bnorth, bdown, bmag ! component of B-field, |B| [Gauss]
   use physics_buffer, only: physics_buffer_desc, pbuf_get_field


!-----------------------------------------------------------------------
! dummy arguments
!-----------------------------------------------------------------------
   integer, intent(in)  :: pcol           ! np. of atmospheric columns
   integer, intent(in)  :: lchnk          ! chunk id
   real(r8), intent(in) :: elfld(3,pcol)  ! electric field on geog. grid geog. direction
!-----------------------------------------------------------------------
! Ion velocities are saved in physics buffer:
!-----------------------------------------------------------------------
   
   type(physics_buffer_desc), pointer :: pbuf(:)

!-----------------------------------------------------------------------
! local variables
!-----------------------------------------------------------------------
   integer  :: i
   real(r8) :: fac
   real(r8), pointer :: ui(:,:), vi(:,:), wi(:,:)

!-----------------------------------------------------------------------
! Set local pointers to respective fields in pbuf:
!-----------------------------------------------------------------------
   call pbuf_get_field(pbuf, ndx_ui, ui )
   call pbuf_get_field(pbuf, ndx_vi, vi )
   call pbuf_get_field(pbuf, ndx_wi, wi )

   do i = 1,pcol    ! number of columns in each chunk
!-----------------------------------------------------------------------
! Magnitude of magnetic field bmag [nT] is from apex module. 
!-----------------------------------------------------------------------
      fac = 1.e9_r8/bmag(i,lchnk)**2           ! nT to T (nT in denominator)
      ui(i,:) =  -(elfld(2,i)*bdown(i,lchnk) + elfld(3,i)*bnorth(i,lchnk))
      vi(i,:) =   elfld(3,i)*beast(i,lchnk) + elfld(1,i)*bdown(i,lchnk)
      wi(i,:) =   elfld(1,i)*bnorth(i,lchnk) - elfld(2,i)*beast(i,lchnk)
      ui(i,:) = ui(i,:)*fac
      vi(i,:) = vi(i,:)*fac
      wi(i,:) = wi(i,:)*fac
   end do ! i = 1,pcol	

   if (state_debug_checks) then
      call shr_assert_in_domain(ui(:pcol,:), is_nan=.false., varname="ui", msg="NaN found in exbdrift::iondrift")
      call shr_assert_in_domain(vi(:pcol,:), is_nan=.false., varname="vi", msg="NaN found in exbdrift::iondrift")
      call shr_assert_in_domain(wi(:pcol,:), is_nan=.false., varname="wi", msg="NaN found in exbdrift::iondrift")
   endif

#ifdef SW_DEBUG
   if( lchnk == 25 ) then
      write(iulog,*) ' '
      write(iulog,*) '---------------------------------------'
      write(iulog,*) 'iondrift: elfld @ lchnk,i = ',lchnk,' 14'
      write(iulog,'(1p,3g15.7)') elfld(:,14)
      write(iulog,*) 'iondrift: bdown,bnorth,beast,bmag @ lchnk,i = ',lchnk,' 14'
      write(iulog,'(1p,4g15.7)') bdown(14,lchnk), bnorth(14,lchnk), beast(14,lchnk), bmag(14,lchnk)
      write(iulog,*) '---------------------------------------'
      write(iulog,*) ' '
   end if
#endif

   end subroutine iondrift

   subroutine cal_mlt( mlt, lchnk, ncol )

!-------------------------------------------------------------------------------
! Purpose: calculate the magnetic local time of WACCM geog. point
!
! Method: using the location of the geomagnetic dipole north pole,
!    the subsolar point location and the apex longitude of the 
!    geographic WACCM point the magn. local time can be calculated
!    subroutines from Roy Barnes HAO Feb. 2004
!
! Author: A. Maute Feb 2004  
!-------------------------------------------------------------------------------

   use mo_apex, only:  &
     alonm,   &  ! apex mag longitude at each geographic grid point (radians) 
     colatp,  &  ! geocentric colatitude of geomagnetic dipole north pole (deg)
     elonp	 ! East longitude of geomagnetic dipole north pole (deg)
   use time_manager, only : get_curr_calday, get_curr_date
   use mo_apex, only : geomag_year
!-------------------------------------------------------------------------------
!	... dummy arguments
!-------------------------------------------------------------------------------
     integer, intent(in)   :: ncol      ! np. of atmospheric columns
     integer, intent(in)   :: lchnk     ! chunk identifier
     real(r8), intent(out) :: mlt(ncol) ! mag.local time of WACCM geo. grid point

!-------------------------------------------------------------------------------
! local arguments
!-------------------------------------------------------------------------------
     integer  :: i
     integer  :: iyear, iday, ihr, imn, imo, iday_m, tod 	! time of day [s] 
     real(r8) :: collonm, sbsllat, sbsllon, sec

!-------------------------------------------------------------------------------
! get current calendar day of year & date components 
! valid at end of current timestep
!-------------------------------------------------------------------------------
     call get_curr_date (iyear,imo,iday_m,tod)  ! year, time of day [sec]
     iyear = int(geomag_year)
     iday  = get_curr_calday()                  ! day of year

     ihr  = tod/3600
     imn  = mod( tod,3600 )/60
     sec  = tod - 60*(ihr*60 + imn)

!-------------------------------------------------------------------------------
!  find subsolar geographic latitude and longitude 
!-------------------------------------------------------------------------------
     call apex_subsol( iyear, iday, ihr, imn, sec, sbsllat, sbsllon )

!-------------------------------------------------------------------------------
!  computes magnetic local time   
!-------------------------------------------------------------------------------
     do i = 1,ncol
       collonm = alonm(i,lchnk)*rtd   ! mag lons (deg)
       call apex_magloctm( collonm, sbsllat, sbsllon, colatp, elonp, mlt(i) )
     end do

     end subroutine cal_mlt

   end module exbdrift
