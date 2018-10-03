
      module mo_aurora
!-----------------------------------------------------------------------
!
! Auroral oval parameterization. See reference:
! R.G. Roble, E.C. Ridley
! An auroral model for the NCAR thermospheric general circulation model (TGCM)
! Annales Geophysicae,5A, (6), 369-382, 1987. 
!
! The aurora oval is a circle in auroral circle coordinates.  Auroral circle
!  coordinates are offset from magnetic coordinates by offa degrees (radians)
!  towards 0 MLT and by dskofa degrees (radians) towards dusk (18 MLT).
! The aurora assumes a Maxwellian in energy, so that the characteristic
!  energy is half of the mean energy (or mean energy = 2*alfa, where alfa
!  is the characteristic energy).  The Maxwellian is approximated in the
!  aion subroutine.
! The aurora oval is assumed to be a Gaussian in auroral latitude, with
!  peak values on the day (=1) and night (=2) sides that change from one to
!  the other using cosines of the auroral longitude coordinate.
! There is provision for a low energy (~75 eV) aurora at the location of the
!  regular (~1-6 keV) aurora in order to simulate the energy flux found
!  at higher altitudes that is non-Maxwellian, but the flux is usually
!  set to zero (1.e-80).
! There is provision for a proton (MeV) aurora, but the flux is usually
!  set to zero (1.e-20).
! The drizzle is a constant low energy electron flux over the polar cap,
!  which goes to 1/e over twice the half-width of the aurora at the
!  radius of the aurora.
! The cusp is a low energy electron flux centered over the dayside convection
!  entrance at phid at the convection reversal boundary theta0.  The cusp
!  falls off over 5 degrees in latitude and over 20 degrees in longitude
!  to 1/e values of the peak at the center.
! 1.e-20 and 1.e-80 are used to give a near zero answer.
!
! The polar drizzle and cusp electron energies are low, and soft particles
!  have great influence on the over-all thermospheric and ionospheric
!  structure, especially on the electron density profiles at mid-latitudes
!  and in winter since low energy electrons produce ionization at high
!  altitudes where loss rates are very low.  (Comment by Wenbin Wang.)
! The original energies for drizzle and cusp were alfad=0.75, alfac=0.5 keV.
! The original guess at energy fluxes were: ed=0.1+2.0*power/100.,ec=0.1+0.9*power/100.
! The next guess at energy fluxes were: ed=0.01+0.2*power/100., ec=0.01+0.09*power/100.
! The values below reflect higher estimates for the electron energy (lower alt)
!
! Calling sequence (all subs in mo_aurora, mo_aurora.F):
!   1) sub aurora_cons called once per time step from advance.
!   2) sub aurora called from dynamics, inside parallel latitude scan.
!   3) subs aurora_cusp and aurora_heat called from sub aurora.
!   4) sub aurora_ions called from sub aurora. 
!
!-----------------------------------------------------------------------

      use shr_kind_mod,  only: r8 => shr_kind_r8
      use mo_constants,  only: pi, gask => rgas_cgs
      use cam_logfile,   only: iulog
      use spmd_utils,    only: masterproc
      use aurora_params, only: power=>hpower, plevel, aurora_params_set
      use aurora_params, only: ctpoten, theta0, dskofa, offa, phid, rrad
      use aurora_params, only: amie_period

      implicit none

      interface aurora
         module procedure aurora_prod
         module procedure aurora_hrate
      end interface

      save


      private
      public :: aurora_inti, aurora_timestep_init, aurora
      public :: aurora_register
      
      integer, parameter  :: isouth = 1
      integer, parameter  :: inorth = 2

      ! g = 8.7 m/s^2? Because this is 400 km up?
      real(r8), parameter :: grav   = 870._r8          ! (cm/s^2)

      integer  :: lev1 = 1
      real(r8) :: rmass_o1
      real(r8) :: rmass_o2
      real(r8) :: rmass_n2
      real(r8) :: rmassinv_o1
      real(r8) :: rmassinv_o2
      real(r8) :: rmassinv_n2
      real(r8), parameter :: twopi = 2._r8*pi
      real(r8), parameter :: d2r = pi/180._r8
      real(r8), parameter :: r2d = 180._r8/pi

!-----------------------------------------------------------------------
! 	... polar drizzle parameters:
!   alfad: Characteristic Maxwellian energy of drizzle electrons (keV)
!   ed   : Column energy input of drizzle electrons (ergs/cm**2/s)
!   fd   : Electron particle flux of drizzle electrons (particles/cm**2/s)
!-----------------------------------------------------------------------
      real(r8), parameter :: alfad = 0.5_r8
      real(r8) :: ed
      real(r8) :: fd                     ! set in sub aurora_ions

!-----------------------------------------------------------------------
! 	... polar cusp parameters:
!   alfac: Characteristic Maxwellian energy of polar cusp electons (keV)
!   ec   : Column energy input of polar cusp electrons (ergs/cm**2/s)
!   fc   : Electron particle flux of polar cusp electrons (particles/cm**2/s)
!-----------------------------------------------------------------------
      real(r8), parameter :: alfac = 0.1_r8
      real(r8) :: ec
      real(r8) :: fc                     ! set in sub aurora_ions

!-----------------------------------------------------------------------
! e1: Peak energy flux in noon sector of the aurora (ergs/cm**2/s)
! e2: Peak energy flux in midnight sector of the aurora (ergs/cm**2/s)
! h1: Gaussian half-width of the noon auroral oval in degrees
! h2: Gaussian half-width of the midnight auroral oval in degrees
!-----------------------------------------------------------------------
      real(r8) :: &
        e1, e2, &                        ! set in sub aurora_cons (function of hem power)
        h1, h2                           ! set in sub aurora_cons (function of hem power)

!-----------------------------------------------------------------------
! 	... additional auroral parameters
!-----------------------------------------------------------------------
      real(r8) :: &
        alfa0, &        ! average of noon and midnight characteristic Maxw energies
        ralfa,ralfa2, & ! difference ratios of characteristic energies
        rrote, &        ! clockwise rotation from noon of peak dayside energy flux (e1)
        rroth, &        ! clockwise rotation from noon of dayside h1 Gaussian half-width
        h0, &           ! average of noon and midnight Gaussian half-widths
        rh, &           ! difference ratio of half-widths (rh=(h2-h1)/(h2+h1))
        e0,e20, &       ! e0 = average of noon and midnight electrons
        ree,re2, &      ! difference ratios of peak energy fluxes (ree=(e2-e1)/(e2+e1))
        alfa20          ! average of noon and midnight char energies for high alt aurora

      logical :: aurora_active = .false.
      integer :: indxAIPRS    = -1
      integer :: indxQTe = -1
      integer :: indxAMIEefxg = -1        ! am_amie_201712
      integer :: indxAMIEkevg = -1        ! am_amie_201712

      real(r8), parameter :: h2deg = 15._r8   ! hour to degree

      contains

        
      !----------------------------------------------------------------------
      !----------------------------------------------------------------------
      subroutine aurora_register
        use ppgrid,       only : pver,pcols
        use physics_buffer, only : pbuf_add_field, dtype_r8

        ! add ionization rates to phys buffer for waccmx ionosphere module

        call pbuf_add_field('AurIPRateSum', 'physpkg', dtype_r8, (/pcols,pver/), indxAIPRS)     ! Sum of ion auroral production rates for O2
        call pbuf_add_field('QTeAur', 'physpkg', dtype_r8, (/pcols/), indxQTe)     ! for electron temperature

      endsubroutine aurora_register

      subroutine aurora_inti(pbuf2d)
!-----------------------------------------------------------------------
! 	... initialize aurora module
!-----------------------------------------------------------------------

      use ppgrid,       only : pver
      use constituents, only : cnst_get_ind, cnst_mw
      use ref_pres,     only : pref_mid
      use mo_chem_utls, only : get_spc_ndx
      use cam_history,  only : addfld, horiz_only
      use physics_buffer,only: pbuf_get_index
      use infnan,       only : nan, assignment(=)
      use physics_buffer, only: physics_buffer_desc, pbuf_set_field

      implicit none

      type(physics_buffer_desc), pointer :: pbuf2d(:,:)
!-----------------------------------------------------------------------
! 	... local variables
!-----------------------------------------------------------------------
      integer             :: k, m
      real(r8), parameter :: e = 1.e-10_r8

      real(r8) :: plb
      real(r8) :: alfa_1, alfa_2, alfa21, alfa22
      real(r8) :: e21, e22
      integer  :: op_ndx,o2p_ndx,np_ndx,n2p_ndx,e_ndx
      integer  :: ierr
      real(r8) :: x_nan

      indxAMIEefxg = pbuf_get_index('AMIE_efxg', errcode=ierr)
      indxAMIEkevg = pbuf_get_index('AMIE_kevg', errcode=ierr)

      if (indxAMIEefxg>0 .and. indxAMIEkevg>0) then
         x_nan = nan
         call pbuf_set_field(pbuf2d, indxAMIEefxg, x_nan)
         call pbuf_set_field(pbuf2d, indxAMIEkevg, x_nan)
      endif

      theta0(:) = nan
      offa(:) = nan
      dskofa(:) = nan
      phid(:) = nan
      rrad(:) = nan
      ctpoten = nan
      power   = nan

      op_ndx   = get_spc_ndx( 'Op' )
      o2p_ndx  = get_spc_ndx( 'O2p' )
      np_ndx   = get_spc_ndx( 'Np' )
      n2p_ndx  = get_spc_ndx( 'N2p' )
      e_ndx    = get_spc_ndx( 'e' )

      aurora_active = op_ndx > 0 .and. o2p_ndx > 0 .and. np_ndx > 0 .and. n2p_ndx > 0 .and. e_ndx > 0 &
                     .and. pref_mid(1) < 0.1_r8 ! need high-top

      if (.not. aurora_active) return

!-----------------------------------------------------------------------
!	... initialize module variables
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!	... set molecular weights
!-----------------------------------------------------------------------
      call cnst_get_ind( 'O2', m )
      rmass_o2    = cnst_mw(m)
      rmassinv_o2 = 1._r8/rmass_o2
      call cnst_get_ind( 'O', m )
      rmass_o1    = cnst_mw(m)
      rmassinv_o1 = 1._r8/rmass_o1
      call cnst_get_ind( 'N', m )
      rmass_n2    = 2._r8*cnst_mw(m)
      rmassinv_n2 = 1._r8/rmass_n2

      offa(isouth)   = 1.0_r8*d2r
      offa(inorth)   = 1.0_r8*d2r

      alfa_1         = 1.5_r8
      alfa_2         = 2._r8

!-----------------------------------------------------------------------
! Values from 10/05/94 HPI estimates (50% or more higher than old estimates):
!     alfa_1 = amin1(1.5,1.25+0.05*plevel)
!     alfa_2 = 1.2 + 0.095*plevel
!-----------------------------------------------------------------------
      alfa0  = 0.5_r8*(alfa_1 + alfa_2)
      ralfa  = (alfa_2 - alfa_1) / (alfa_1 + alfa_2 + e)
      alfa21 = 0.075_r8
      alfa22 = 0.075_r8
      alfa20 = 0.5_r8 * (alfa21 + alfa22)
      ralfa2 = (alfa22 - alfa21) / (alfa21 + alfa22 + e)
      e21    = 1.e-80_r8
      e22    = 1.e-80_r8
      e20    = 0.5_r8 * (e21 + e22)
      re2    = (e22 - e21) / (e21 + e22)

!-----------------------------------------------------------------------
! 	... set auroral lower bndy index
!-----------------------------------------------------------------------
      plb = 5.e-4_r8*exp( 7._r8 ) * .1_r8             ! Pa
      do k = 1,pver
	 if( pref_mid(k) >= plb ) then
	    lev1 = k-1
	    exit
	 end if
      end do

      if (masterproc) write(iulog,*) ' '
      if (masterproc) write(iulog,*) 'aurora_inti: aurora will go down to lev,p = ',lev1,pref_mid(lev1)
      if (masterproc) write(iulog,*) ' '

!-----------------------------------------------------------------------
! Report to stdout:
!-----------------------------------------------------------------------
#ifdef AURORA_DIAGS
        write(iulog,"(/,'aurora_cons:')")
!        write(iulog,"('  cusp:    alfac=',f8.3,' ec=',f8.3,' fc=',e10.4)") &
!          alfac,ec,fc
!        write(iulog,"('  drizzle: alfad=',f8.3,' ed=',f8.3,' fd=',e10.4)") &
!          alfad,ed,fd
        write(iulog,"('  half-widths = h1,h2=',2f10.3)") h1,h2
        write(iulog,"('  energy flux = e1,e2=',2f10.3)") e1,e2
        write(iulog,"('  add_sproton = ',l1)") add_sproton
        write(iulog,"(' ')")
#endif
        call addfld('ALATM', horiz_only,  'I','degrees', &
             'Magnetic latitude at each geographic coordinate')
        call addfld('ALONM', horiz_only,  'I','degrees', &
             'Magnetic longitude at each geographic coordinate')
        call addfld( 'QSUM', (/ 'lev' /), 'I','/s',      &
             'total ion production' )

      end subroutine aurora_inti

      subroutine aurora_timestep_init( )
!-----------------------------------------------------------------------
! 	... per timestep initialization
!-----------------------------------------------------------------------

      use heelis_mod,    only : heelis_update

!-----------------------------------------------------------------------
!	... local variables
!-----------------------------------------------------------------------
      real(r8) :: roth, rote
      real(r8), parameter :: convert = 3.1211e8_r8

      if (.not. aurora_active) return

      if (.not.aurora_params_set) then
         call heelis_update() ! use heelis if aurora parameters are not already set
      endif

      if( power >= 1.0_r8 ) then
         plevel = 2.09_r8*log( power )
      else
         plevel = 0._r8
      end if

#ifdef AURORA_DIAGS
      if( masterproc ) then
         write(iulog,*) '----------------------------------------'
         write(iulog,*) 'aurora_timestep_init: power,ctpoten = ',power,ctpoten
         write(iulog,*) '----------------------------------------'
      end if
#endif

!-----------------------------------------------------------------------
! h1 = Gaussian half-width of the noon auroral oval in degrees
! h2 = Gaussian half-width of the midnight auroral oval in degrees
!-----------------------------------------------------------------------
! produce realistic oval compared to NOAA empirical auroral oval and TIMED/GUVI
! h1 formula given by Wenbin base on POLARVIS image;
! h2 formula based on Emery et al original auroral parameterization report
      h1 = min(2.35_r8, 0.83_r8 + 0.33_r8*plevel)
      h2 = 2.5_r8+0.025_r8*max(power,55._r8)+0.01_r8*min(0._r8,power-55._r8)

!-----------------------------------------------------------------------
! Values from corrections to Emery et al Parameterization report:
!     h1 = amin1(2.35, 0.83 + 0.33*plevel)
!     h2 = 2.87 + 0.15*plevel
!-----------------------------------------------------------------------

      rh = (h2 - h1) / (h1 + h2)
      h0     = 0.5_r8 * (h1 + h2) * d2r   
      
      
! roth = MLT of max width of aurora in hours
! rote = MLT of max energy flux of aurora in hours

      roth = 0.81_r8 - 0.06_r8 * plevel
      rote = 0.17_r8 - 0.04_r8 * plevel

!  Convert MLT from hours to degrees to radians

      rroth = roth * h2deg * d2r
      rrote = rote * h2deg * d2r

!-----------------------------------------------------------------------
! e1 = energy flux in the noon sector of the aurora (ergs/cm**2/s)
! e2 = energy flux in the midnight sector of the aurora (ergs/cm**2/s)
!-----------------------------------------------------------------------
! produce realistic oval compared to NOAA empirical auroral oval and TIMED/GUVI
! e1 formula given by Wenbin base on POLARVIS image;
! e2 formula based on Emery et al original auroral parameterization report
!-----------------------------------------------------------------------
      e1 = max(0.50_r8, -2.15_r8 + 0.62_r8*plevel)
      e2=1._r8+0.11_r8*power

!-----------------------------------------------------------------------
!   ed   : Column energy input of drizzle electrons (ergs/cm**2/s)
!   ec   : Column energy input of polar cusp electrons (ergs/cm**2/s)
!-----------------------------------------------------------------------
      ed = .0012_r8+.0006_r8*power
      ec = (0.24_r8+0.0067_r8*power)/5._r8

!-----------------------------------------------------------------------
! Set cusp and drizzle parameters:
! (conversion between particle number density and characteristic
!  energy and column energy input)
!-----------------------------------------------------------------------
      fc = convert * ec / alfac
      fd = convert * ed / alfad

!-----------------------------------------------------------------------
! Values from corrections to Emery et al Parameterization report:
!-----------------------------------------------------------------------
      e0  = 0.5_r8 * (e1 + e2)
      ree = (e2 - e1) / (e1 + e2)

      end subroutine aurora_timestep_init

      subroutine aurora_prod( tn, o2, o1, mbar, rlats, &
                              qo2p, qop, qn2p, qnp, pmid, &
                              lchnk, calday,  ncol, rlons, pbuf )
!-----------------------------------------------------------------------
! 	... auroral parameterization driver
!-----------------------------------------------------------------------

      use mo_apex,     only : alatm, alonm                      ! magnetic latitude,longitude grid (radians)
      use mo_apex,     only : maglon0
      use ppgrid,      only : pcols, pver
      use cam_history, only : outfld
      use physics_buffer,only: physics_buffer_desc,pbuf_get_field

      implicit none

!-----------------------------------------------------------------------
! 	... dummy arguments
!-----------------------------------------------------------------------
      integer, intent(in) ::  &
       ncol, &                           ! column count
       lchnk                             ! chunk index
      real(r8), intent(in) :: &
        calday                           ! calendar day of year
      real(r8), intent(in) :: &
        tn(pcols,pver), &                ! neutral gas temperature (K)
        o2(ncol,pver), &                 ! O2 concentration (kg/kg)
        o1(ncol,pver), &                 ! O concentration (kg/kg)
        mbar(ncol,pver)                  ! mean molecular weight (g/mole)
      real(r8), intent(in) :: &
        pmid(pcols,pver)                 ! midpoint pressure (Pa)
      real(r8), intent(in) :: &
        rlats(ncol), &                   ! column latitudes (radians)
        rlons(ncol)
      real(r8), intent(out) :: &
        qo2p(ncol,pver), &               ! o2+ production
        qop(ncol,pver), &                ! o+ production
        qn2p(ncol,pver), &               ! n2+ production
        qnp(ncol,pver)                   ! n+ production

      type(physics_buffer_desc),pointer :: pbuf(:)

!-----------------------------------------------------------------------
! 	... local variables
!-----------------------------------------------------------------------
      integer  :: i, k
      integer  :: hemis(ncol)
      real(r8) :: ofda, cosofa, sinofa, aslona
      real(r8) :: dlat_aur(ncol)
      real(r8) :: dlon_aur(ncol)
      real(r8) :: colat(ncol)
      real(r8) :: sinlat(ncol)
      real(r8) :: coslat(ncol)
      real(r8) :: coslon(ncol)
      real(r8) :: sinlon(ncol)
      real(r8) :: alon(ncol)
      real(r8) :: cusp(ncol)
      real(r8) :: alfa(ncol)
      real(r8) :: alfa2(ncol)
      real(r8) :: flux(ncol)
      real(r8) :: flux2(ncol)
      real(r8) :: drizl(ncol)
      logical  :: do_aurora(ncol)

      real(r8) :: dayfrac, rotation
      real(r8), pointer :: qteaur(:)    ! for electron temperature

      if (.not. aurora_active) return

!-----------------------------------------------------------------------
! 	... initialize ion production
!-----------------------------------------------------------------------
      do k = 1,pver
        qo2p(:,k) = 0._r8
        qop(:,k)  = 0._r8
        qn2p(:,k) = 0._r8
        qnp(:,k)  = 0._r8
      end do

!-----------------------------------------------------------------------
! 	... output mag lons, lats
!-----------------------------------------------------------------------
      call outfld( 'ALONM', r2d*alonm(:ncol,lchnk), pcols, lchnk )
      call outfld( 'ALATM', r2d*alatm(:ncol,lchnk), pcols, lchnk )

      if (indxQTe>0) then
        call pbuf_get_field(pbuf, indxQTe, qteaur)
        qteaur(:) = 0._r8
     endif
     
!-----------------------------------------------------------------------
!    aurora is active for columns poleward of 30 deg
!-----------------------------------------------------------------------
      do_aurora(:) = abs( rlats(:) ) > pi/6._r8
      if( all( .not. do_aurora(:) ) ) then
         return
      end if

!-----------------------------------------------------------------------
! 	... set rotation based on sun location
!-----------------------------------------------------------------------
      dayfrac = (calday - int(calday))
      rotation = maglon0 + dayfrac*twopi-pi

      do i = 1,ncol
        if( do_aurora(i) ) then
          dlat_aur(i) = alatm(i,lchnk)
          dlon_aur(i) = alonm(i,lchnk) + rotation ! rotate it 
          if( dlon_aur(i) > pi ) then
             dlon_aur(i) = dlon_aur(i) - twopi
          else if( dlon_aur(i) < -pi ) then
             dlon_aur(i) = dlon_aur(i) + twopi
          end if
          if( dlat_aur(i) > 0._r8 ) then
            hemis(i) = 2
          else
            hemis(i) = 1
          end if
!-----------------------------------------------------------------------
! 	... find auroral circle coordinates
!-----------------------------------------------------------------------
            ofda      = sqrt( offa(hemis(i))**2 + dskofa(hemis(i))**2)
            cosofa    = cos( ofda )
            sinofa    = sin( ofda )
            aslona    = asin( dskofa(hemis(i))/ofda )
            sinlat(i) = sin( abs( dlat_aur(i) ) )
            coslat(i) = cos( dlat_aur(i) )
            sinlon(i) = sin( dlon_aur(i) + aslona )
            coslon(i) = cos( dlon_aur(i) + aslona )
            colat(i)  = acos( cosofa*sinlat(i) - sinofa*coslat(i)*coslon(i))
            alon(i)   = mod( atan2( sinlon(i)*coslat(i),sinlat(i)*sinofa &
                                    + cosofa*coslat(i)*coslon(i) ) - aslona + 3._r8*pi,twopi) - pi
        end if
      end do
#ifdef AURORA_DIAGS
      write(iulog,*) '-----------------------------------------------------'
      write(iulog,*) 'aurora: diagnostics for lchnk = ',lchnk
      write(iulog,*) '        geo lats'
      write(iulog,'(1p,5g15.7)') r2d*rlats(:ncol)
      write(iulog,*) '        geo lons'
      write(iulog,'(1p,5g15.7)') r2d*rlons(:ncol)
      write(iulog,*) '        mag lats'
      write(iulog,'(1p,5g15.7)') r2d*dlat_aur(:ncol)
      write(iulog,*) '        mag lons'
      write(iulog,'(1p,5g15.7)') r2d*alonm(:ncol,lchnk)
      write(iulog,*) '        mag table lons'
      write(iulog,'(1p,5g15.7)') r2d*dlon_aur(:ncol)
      write(iulog,*) '     min,max mag lons = ',r2d*minval(alonm(:ncol,lchnk)),r2d*maxval(alonm(:ncol,lchnk))
      write(iulog,*) '-----------------------------------------------------'
#endif

!-----------------------------------------------------------------------
! 	... make cusp
!-----------------------------------------------------------------------
      call aurora_cusp( cusp, do_aurora, hemis, colat, alon, ncol )

!-----------------------------------------------------------------------
! 	... make alfa, flux, and drizzle
!-----------------------------------------------------------------------
      call aurora_heat( flux, flux2, alfa, alfa2, &
                        drizl, do_aurora, hemis, &
                        alon, colat, ncol, pbuf )

!-----------------------------------------------------------------------
! 	... auroral additions to ionization rates
!-----------------------------------------------------------------------
      call aurora_ions( drizl, cusp, alfa, alfa2, &
                        flux, flux2, tn, o2, &
                        o1, mbar, qo2p, qop, qn2p, &
                        qnp, pmid, do_aurora, ncol, lchnk, pbuf )

      end subroutine aurora_prod

      subroutine aurora_hrate( tn, mbar, rlats, &
                               aur_hrate, cpair, pmid, lchnk, calday, &
                               ncol, rlons, pbuf  )
!-----------------------------------------------------------------------
! 	... auroral parameterization driver
!-----------------------------------------------------------------------

      use mo_apex, only : alatm, alonm                      ! magnetic latitude,longitude grid (radians)
      use mo_apex, only : maglon0
      use ppgrid,  only : pcols, pver
      use physics_buffer,only: physics_buffer_desc

      implicit none

!-----------------------------------------------------------------------
! 	... dummy arguments
!-----------------------------------------------------------------------
      integer, intent(in) ::  &
       ncol, &                           ! column count
       lchnk                             ! chunk index
      real(r8), intent(in) :: &
        calday                           ! calendar day of year
      real(r8), intent(in) :: &
        tn(pcols,pver), &                ! neutral gas temperature (K)
        mbar(ncol,pver)                  ! mean molecular weight (g/mole)
      real(r8), intent(in) :: &
        cpair(ncol,pver)                 ! specific heat capacity (J/K/kg)
      real(r8), intent(in) :: &
        pmid(pcols,pver)                 ! midpoint pressure (Pa)
      real(r8), intent(in) :: &
        rlats(ncol), &                   ! column latitudes (radians)
        rlons(ncol)
      real(r8), intent(out) :: &
        aur_hrate(ncol,pver)             ! auroral heating rate
      type(physics_buffer_desc),pointer :: pbuf(:)

!-----------------------------------------------------------------------
! 	... local variables
!-----------------------------------------------------------------------
      real(r8), parameter :: aur_therm     = 807._r8
      real(r8), parameter :: jkcal         = 4184._r8
      real(r8), parameter :: aur_heat_eff  = .05_r8
      real(r8), parameter :: aur_hconst    = 1.e3_r8*jkcal*aur_therm*aur_heat_eff

      integer  :: i, k
      integer  :: hemis(ncol)
      real(r8) :: ofda, cosofa, sinofa, aslona
      real(r8) :: dlat_aur(ncol)
      real(r8) :: dlon_aur(ncol)
      real(r8) :: colat(ncol)
      real(r8) :: sinlat(ncol)
      real(r8) :: coslat(ncol)
      real(r8) :: coslon(ncol)
      real(r8) :: sinlon(ncol)
      real(r8) :: alon(ncol)
      real(r8) :: cusp(ncol)
      real(r8) :: alfa(ncol)
      real(r8) :: alfa2(ncol)
      real(r8) :: flux(ncol)
      real(r8) :: flux2(ncol)
      real(r8) :: drizl(ncol)
      real(r8) :: qsum(ncol,pver)                      ! total ion production (1/s)
      logical  :: do_aurora(ncol)

      real(r8) :: dayfrac, rotation

!-----------------------------------------------------------------------
! 	... initialize ion production
!-----------------------------------------------------------------------
      do k = 1,pver
        aur_hrate(:,k) = 0._r8
      end do

      if (.not. aurora_active) return

!-----------------------------------------------------------------------
! 	... check latitudes, and return if all below 32.5 deg
!-----------------------------------------------------------------------
      do_aurora(:) = abs( rlats(:) ) > pi/6._r8
      if( all( .not. do_aurora(:) ) ) then
         return
      end if

!-----------------------------------------------------------------------
! 	... set rotation based on sun location
!-----------------------------------------------------------------------
      dayfrac = (calday - int(calday))
      rotation = maglon0 + dayfrac*twopi-pi

      do i = 1,ncol
        if( do_aurora(i) ) then
          dlat_aur(i) = alatm(i,lchnk)
          dlon_aur(i) = alonm(i,lchnk) + rotation ! rotate it 
          if( dlon_aur(i) > pi ) then
             dlon_aur(i) = dlon_aur(i) - twopi
          else if( dlon_aur(i) < -pi ) then
             dlon_aur(i) = dlon_aur(i) + twopi
          end if
          if( dlat_aur(i) > 0._r8 ) then
            hemis(i) = 2
          else
            hemis(i) = 1
          end if
!-----------------------------------------------------------------------
! 	... find auroral circle coordinates
!-----------------------------------------------------------------------
            ofda      = sqrt( offa(hemis(i))**2 + dskofa(hemis(i))**2)
            cosofa    = cos( ofda )
            sinofa    = sin( ofda )
            aslona    = asin( dskofa(hemis(i))/ofda )
            sinlat(i) = sin( abs( dlat_aur(i) ) )
            coslat(i) = cos( dlat_aur(i) )
            sinlon(i) = sin( dlon_aur(i) + aslona )
            coslon(i) = cos( dlon_aur(i) + aslona )
            colat(i)  = acos( cosofa*sinlat(i) - sinofa*coslat(i)*coslon(i))
            alon(i)   = mod( atan2( sinlon(i)*coslat(i),sinlat(i)*sinofa &
                                    + cosofa*coslat(i)*coslon(i) ) - aslona + 3._r8*pi,twopi) - pi
        end if
      end do
#ifdef AURORA_DIAGS
      write(iulog,*) '-----------------------------------------------------'
      write(iulog,*) 'aurora: diagnostics for lchnk = ',lchnk
      write(iulog,*) '        geo lats'
      write(iulog,'(1p,5g15.7)') r2d*rlats(:ncol)
      write(iulog,*) '        geo lons'
      write(iulog,'(1p,5g15.7)') r2d*rlons(:ncol)
      write(iulog,*) '        mag lats'
      write(iulog,'(1p,5g15.7)') r2d*dlat_aur(:ncol)
      write(iulog,*) '        mag lons'
      write(iulog,'(1p,5g15.7)') r2d*alonm(:ncol,lchnk)
      write(iulog,*) '        mag table lons'
      write(iulog,'(1p,5g15.7)') r2d*dlon_aur(:ncol)
      write(iulog,*) '     min,max mag lons = ',r2d*minval(alonm(:ncol,lchnk)),r2d*maxval(alonm(:ncol,lchnk))
      write(iulog,*) '-----------------------------------------------------'
#endif

!-----------------------------------------------------------------------
! 	... make cusp
!-----------------------------------------------------------------------
      call aurora_cusp( cusp, do_aurora, hemis, colat, alon, ncol )

!-----------------------------------------------------------------------
! 	... make alfa, flux, and drizzle
!-----------------------------------------------------------------------
      call aurora_heat( flux, flux2, alfa, alfa2, &
                        drizl, do_aurora, hemis, &
                        alon, colat, ncol, pbuf )
 
!-----------------------------------------------------------------------
! 	... auroral additions to ionization rates
!-----------------------------------------------------------------------
      call total_ion_prod( drizl, cusp, alfa, alfa2, &
                           flux, flux2, tn, &
                           mbar, qsum, pmid, do_aurora, &
                           ncol )

!-----------------------------------------------------------------------
! 	... form auroral heating rate
!-----------------------------------------------------------------------
      do k = 1,pver
         aur_hrate(:,k) = aur_hconst * qsum(:,k) / (cpair(:,k) * mbar(:,k))
      end do

      end subroutine aurora_hrate

      subroutine aurora_cusp( cusp, do_aurora, hemis, colat, alon, ncol )
!-----------------------------------------------------------------------
! 	... calculate horizontal variation of polar cusp heating
!-----------------------------------------------------------------------

      implicit none

!-----------------------------------------------------------------------
! 	... dummy arguments
!-----------------------------------------------------------------------
      integer, intent(in)     :: ncol
      integer, intent(in)     :: hemis(ncol)
      real(r8), intent(in)    :: colat(ncol)
      real(r8), intent(in)    :: alon(ncol)
      real(r8), intent(out)   :: cusp(ncol)
      logical, intent(in)     :: do_aurora(ncol)

!-----------------------------------------------------------------------
! 	... local variables
!-----------------------------------------------------------------------
      real(r8), parameter :: s5  =.08726646_r8, &
                             s20 =.34906585_r8

      where( do_aurora(:) )
         cusp(:) = (exp( -((theta0(hemis(:)) - colat(:))/s5)**2 ) &
                      + exp( -((pi - theta0(hemis(:)) - colat(:))/s5)**2) ) &
                        *exp( -(atan2( sin(alon(:) - phid(hemis(:))), cos(alon(:) - phid(hemis(:))) )/s20)**2 )
      elsewhere
         cusp(:) = 0._r8
      endwhere

      end subroutine aurora_cusp 

      subroutine aurora_heat( flux, flux2, alfa, alfa2, &
                              drizl, do_aurora, hemis, &
                              alon, colat, ncol, pbuf )
!-----------------------------------------------------------------------
! 	... calculate alfa, flux, and drizzle
!-----------------------------------------------------------------------
      use physics_buffer,only: physics_buffer_desc,pbuf_get_field
      
      implicit none

!-----------------------------------------------------------------------
! 	... dummy arguments
!-----------------------------------------------------------------------
      integer, intent(in)     :: ncol
      integer, intent(in)     :: hemis(ncol)
      real(r8), intent(in)    :: colat(ncol)
      real(r8), intent(in)    :: alon(ncol)
      real(r8), intent(inout) :: flux(ncol)
      real(r8), intent(inout) :: flux2(ncol)
      real(r8), intent(inout) :: drizl(ncol)
      real(r8), intent(inout) :: alfa(ncol)
      real(r8), intent(inout) :: alfa2(ncol)
      logical, intent(in)     :: do_aurora(ncol)
      type(physics_buffer_desc),pointer :: pbuf(:)

!-----------------------------------------------------------------------
! 	... local variables
!-----------------------------------------------------------------------
      real(r8), dimension(ncol) :: &
        coslamda, &                             ! cos(angle from throat)
        halfwidth, &                            ! oval half-width
        wrk, &                                  ! temp wrk array
        dtheta                                  ! latitudinal variation (Gaussian)
      real(r8) :: ekev 
      real(r8), pointer   :: amie_efxg(:) ! Pointer to pbuf AMIE energy flux (mW m-2)
      real(r8), pointer   :: amie_kevg(:) ! Pointer to pbuf AMIE mean energy (keV)
      real(r8), pointer   :: qteaur(:)    ! for electron temperature
      integer  :: n
      
!-----------------------------------------------------------------------
! Low-energy protons:
!
!     alfap0 = 0.5*(alfap1+alfap2)
!     e0p = 0.5*(pe1+pe2)
!
! coslamda  = cos(lamda)
! halfwidth = auroral half width
! dtheta    = colat-theta0(ihem)
! alfa      = electron energy
!-----------------------------------------------------------------------
      where( do_aurora(:) )
         coslamda(:) = cos( atan2( sin( alon(:) - rrote ),cos( alon(:) - rrote ) ) )
!-----------------------------------------------------------------------
! 	... auroral oval half-width (equation (1) in Roble,1987):
!-----------------------------------------------------------------------
         halfwidth(:) = h0*(1._r8 - rh*cos( atan2( sin(alon(:) - rroth),cos( alon(:) - rroth ) ) ) )
         dtheta(:)    = colat(:) - rrad(hemis(:))
      endwhere
!-----------------------------------------------------------------------
! 	... characteristic energy (equation (2) in Roble,1987):
!-----------------------------------------------------------------------
      if( alfa0 > .01_r8 ) then
         where( do_aurora(:) )
            alfa(:) =  alfa0*(1._r8 - ralfa*coslamda(:))
         endwhere
      else
         alfa(:) =  0._r8
      end if

      wrk = 0._r8
      where( do_aurora(:) )
!-----------------------------------------------------------------------
! 	... flux, drizzle, alfa2, flux2
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! 	... energy flux (equation (3) in Roble,1987):
!-----------------------------------------------------------------------
         wrk(:)   = exp( -(dtheta(:)/halfwidth(:))**2 )
         flux(:)  = e0*(1._r8 - ree*coslamda(:))*wrk(:) / (2._r8*alfa(:)*1.602e-9_r8)
         drizl(:) = exp( -((dtheta(:) + abs(dtheta(:)))/(2._r8*h0))**2 )
         alfa2(:) = alfa20*(1._r8 - ralfa2*coslamda(:))
         flux2(:) = e20*(1._r8 - re2*coslamda(:))*wrk(:) / (2._r8*alfa2(:)*1.602e-9_r8)
      endwhere

!-----------------------------------------------------------------------
! 	... for electron temperature (used in settei):  
!-----------------------------------------------------------------------
      if (indxQTe>0) then
         call pbuf_get_field(pbuf, indxQTe, qteaur)
         ! The factor of -7e8 is based on the energy per ion pair and heating efficiency from Roble (1987)
         qteaur(:ncol) = -7.e8_r8*wrk(:ncol)
      end if

!----------------------------------------------------------------------------------------------
!       ... If turned on, use amie energy flux and mean energy to replace flux(:) and alfa(:)
!----------------------------------------------------------------------------------------------
      if (amie_period .and. indxAMIEefxg>0 .and. indxAMIEkevg>0) then
         !---------------------------------------------------------------------------
         ! Overwrite with AMIE mean energy and energy flux in physics buffer
         !---------------------------------------------------------------------------
         call pbuf_get_field(pbuf, indxAMIEefxg, amie_efxg)
         call pbuf_get_field(pbuf, indxAMIEkevg, amie_kevg)
         do n=1,ncol
            ekev = max(amie_kevg(n),1._r8)
            alfa(n) = ekev/2._r8
            flux(n) = max(amie_efxg(n)/(ekev*1.602e-9_r8),1.e-20_r8)
         enddo
      endif

      end subroutine aurora_heat

      subroutine aurora_ions( drizl, cusp, alfa1, alfa2, &
                              flux1, flux2, tn, o2, &
                              o1, mbar, qo2p, qop, qn2p, &
                              qnp, pmid, do_aurora, ncol, lchnk, pbuf )
!-----------------------------------------------------------------------
! 	... calculate auroral additions to ionization rates
!-----------------------------------------------------------------------

      use ppgrid,      only : pcols, pver
      use cam_history, only : outfld

      use physics_buffer,only: physics_buffer_desc, pbuf_get_field

      implicit none

!-----------------------------------------------------------------------
! 	... dummy arguments
!-----------------------------------------------------------------------
      integer, intent(in) :: ncol
      integer, intent(in) :: lchnk
      real(r8), intent(in), dimension(ncol) :: &
                             drizl, &
                             cusp, &
                             alfa1, &
                             alfa2, &
                             flux1, &
                             flux2
      real(r8), dimension(pcols,pver), intent(in) :: &
                             tn, &                     ! midpoint neutral temperature (K)
                             pmid                      ! midpoint pressure (Pa)
      real(r8), dimension(ncol,pver), intent(in) :: &
                             o2, &                     ! midpoint o2 concentration (kg/kg)
                             o1, &                     ! midpoint o  concentration (kg/kg)
                             mbar                      ! mean molecular mass (g/mole)
      real(r8), dimension(ncol,pver), intent(inout) :: &
                             qo2p, &                   ! o2p prod from aurora (molecules/cm^3/s)
                             qop, &                    ! op prod from aurora (molecules/cm^3/s)
                             qn2p, &                   ! n2p prod from aurora (molecules/cm^3/s)
                             qnp                       ! np prod from aurora (molecules/cm^3/s)
      logical, intent(in) :: do_aurora(ncol)

      type(physics_buffer_desc),pointer :: pbuf(:)

!-----------------------------------------------------------------------
! 	... local variables
!-----------------------------------------------------------------------
      real(r8), parameter :: const0        = 1.e-20_r8

      integer  :: k
      real(r8), dimension(ncol) :: &
        p0ez, &
        press, &                                   ! pressure at interface levels (dyne/cm^2)
        tempi, &                                   ! temperature at interface levels (K)
        xalfa1, &
        xalfa2, &
        xcusp, &
        xdrizl, &                                  ! input to sub aion
        cusp_ion, &
        drizl_ion, &                               ! output from sub aion
        alfa1_ion, &
        alfa2_ion, &
        barm_t, &
        qsum, &
        denom, &
        barm, &
        falfa1, &
        falfa2, &
        fcusp, &
        fdrizl, &
        xn2
      real(r8), dimension(ncol) :: &
        qo2p_aur, &
        qop_aur, &
        qn2p_aur                                   ! auroral ionization for O2+, O+, N2+
      real(r8) :: qia(5)                           ! low energy proton source (not in use, 1/02)
      real(r8) :: wrk(ncol,pver)

      real(r8), pointer   :: aurIPRateSum(:,:) ! Pointer to pbuf auroral ion production sum for O2+,O+,N2+ (s-1 cm-3)
      
      qia(:) = 0._r8
      wrk(:,:) = 0._r8

      !-----------------------------------------------------------
      !  Point to production rates array in physics buffer where 
      !  rates will be stored for ionosphere module access.  Also, 
      !  initialize rates to zero before column loop since only 
      !  daylight values are filled
      !-----------------------------------------------------------
      if (indxAIPRS>0) then
        call pbuf_get_field(pbuf, indxAIPRS, aurIPRateSum)
        aurIPRateSum(:,:) = 0._r8
      endif

level_loop : &
      do k = 1,lev1
          where( do_aurora(:) )
             press(:ncol) = 10._r8*pmid(:ncol,k)              ! from Pa to dyne/cm^2
             tempi(:ncol) = tn(:ncol,k)
             barm(:)      = mbar(:,k)
             p0ez(:)      = (press(:)/(grav*4.e-6_r8))**.606_r8
             xalfa1(:)    = p0ez(:)/alfa1(:)
             xalfa2(:)    = p0ez(:)/alfa2(:)
             xcusp (:)    = p0ez(:)/alfac
             xdrizl(:)    = p0ez(:)/alfad

!-----------------------------------------------------------------------
! 	... initialize (whole array operations):
!-----------------------------------------------------------------------
             alfa1_ion(:) = const0
             alfa2_ion(:) = const0
             cusp_ion(:)  = const0
             drizl_ion(:) = const0
          endwhere
!-----------------------------------------------------------------------
! 	... auroral electrons
!-----------------------------------------------------------------------
          call aion( xalfa1, alfa1_ion, do_aurora, ncol )
          call aion( xalfa2, alfa2_ion, do_aurora, ncol )
          call aion( xcusp , cusp_ion, do_aurora, ncol  )
          call aion( xdrizl, drizl_ion, do_aurora, ncol )
          where( do_aurora(:) )
             falfa1(:) = alfa1(:)*flux1(:)  ! s7
             falfa2(:) = alfa2(:)*flux2(:)  ! s8
             fcusp (:) = cusp(:)*alfac*fc   ! s9
             fdrizl(:) = drizl(:)*alfad*fd  ! s10
             qsum(:)   = falfa1(:)*alfa1_ion(:) &    ! s7*s3
                       + falfa2(:)*alfa2_ion(:) &    ! s8*s4
                       + fcusp(:)*cusp_ion (:) &     ! s9*s5
                       + fdrizl(:)*drizl_ion(:)       ! s10*s6
          endwhere

!-----------------------------------------------------------------------
! 	... form production
!-----------------------------------------------------------------------
          where( do_aurora(:) )
             barm_t(:) = grav*barm(:)/(35.e-3_r8*gask*tempi(:))
             qsum(:)   = qsum(:)*barm_t(:)               ! s1 = s1*s11
             wrk(:,k)  = qsum(:)
!-----------------------------------------------------------------------
! 	... denominator of equations (13-16) in Roble,1987.
!-----------------------------------------------------------------------
             xn2(:)   = max( (1._r8 - o2(:,k) - o1(:,k)),1.e-8_r8 )
             denom(:) = 0.92_r8*xn2(:)*rmassinv_n2 &
                      + 1.5_r8*o2(:,k) *rmassinv_o2 + 0.56_r8*o1(:,k) *rmassinv_o1
!-----------------------------------------------------------------------
! 	... production of O2+ (equation (15) in Roble,1987):
!-----------------------------------------------------------------------
             qo2p_aur(:) = qsum(:)*o2(:,k)/(rmass_o2*denom(:)) + qia(2)
!-----------------------------------------------------------------------
! 	... production of O+ (equation (16) in Roble,1987):
!-----------------------------------------------------------------------
             qop_aur(:) = qsum(:)*(.5_r8 *o2(:,k)*rmassinv_o2 &
                                   + .56_r8*o1(:,k)*rmassinv_o1)/denom(:) + qia(3)
!-----------------------------------------------------------------------
! 	... production of N2+ (equation (13) in Roble,1987)
!-----------------------------------------------------------------------
             qn2p_aur(:) = qsum(:)*.7_r8*xn2(:)/(rmass_n2*denom(:)) + qia(1)
             qo2p(:,k)   = qo2p(:,k) + qo2p_aur(:)
             qop(:,k)    = qop(:,k) + qop_aur(:)
             qn2p(:,k)   = qn2p(:,k) + qn2p_aur(:)
             qnp(:,k)    = qnp (:,k) + .22_r8/.7_r8 * qn2p_aur(:)
          endwhere
      end do level_loop

      !----------------------------------------------------------------
      !  Store the sum of the ion production rates in pbuf to be used 
      !  in the ionosx module 
      !----------------------------------------------------------------
      if (indxAIPRS>0) then
      
        aurIPRateSum(1:ncol,1:pver) = wrk(1:ncol,1:pver) 
      
      endif
  
      call outfld( 'QSUM', wrk, ncol, lchnk )

      end subroutine aurora_ions

      subroutine total_ion_prod( drizl, cusp, alfa1, alfa2, &
                                 flux1, flux2, tn, &
                                 mbar, tpions, pmid, do_aurora, &
                                 ncol )
!-----------------------------------------------------------------------
! 	... calculate auroral additions to ionization rates
!-----------------------------------------------------------------------

      use ppgrid,      only : pcols, pver

      implicit none

!-----------------------------------------------------------------------
! 	... dummy arguments
!-----------------------------------------------------------------------
      integer, intent(in) :: ncol
      real(r8), intent(in), dimension(ncol) :: &
                             drizl, &
                             cusp, &
                             alfa1, &
                             alfa2, &
                             flux1, &
                             flux2
      real(r8), dimension(pcols,pver), intent(in) :: &
                             tn, &                     ! midpoint neutral temperature (K)
                             pmid                      ! midpoint pressure (Pa)
      real(r8), dimension(ncol,pver), intent(in) :: &
                             mbar                      ! mean molecular mass (g/mole)
      real(r8), dimension(ncol,pver), intent(inout) :: &
                             tpions                    ! total ion production (1/s)
      logical, intent(in) :: do_aurora(ncol)

!-----------------------------------------------------------------------
! 	... local variables
!-----------------------------------------------------------------------
      real(r8), parameter :: const0        = 1.e-20_r8

      integer  :: k
      real(r8), dimension(ncol) :: &
        p0ez, &
        press, &                                   ! pressure at interface levels (dyne/cm^2)
        tempi, &                                   ! temperature at interface levels (K)
        xalfa1, &
        xalfa2, &
        xcusp, &
        xdrizl, &                                  ! input to sub aion
        cusp_ion, &
        drizl_ion, &                               ! output from sub aion
        alfa1_ion, &
        alfa2_ion, &
        barm_t, &
        qsum, &
        barm, &
        falfa1, &
        falfa2, &
        fcusp

      tpions(:,:) = 0._r8

level_loop : &
      do k = 1,lev1
          where( do_aurora(:) )
             press(:ncol) = 10._r8*pmid(:ncol,k)              ! from Pa to dyne/cm^2
             tempi(:ncol) = tn(:ncol,k)
             barm(:)      = mbar(:,k)
             p0ez(:)      = (press(:)/(grav*4.e-6_r8))**.606_r8
             xalfa1(:)    = p0ez(:)/alfa1(:)
             xalfa2(:)    = p0ez(:)/alfa2(:)
             xcusp (:)    = p0ez(:)/alfac
             xdrizl(:)    = p0ez(:)/alfad

!-----------------------------------------------------------------------
! 	... initiliaze (whole array operations):
!-----------------------------------------------------------------------
             alfa1_ion(:) = const0
             alfa2_ion(:) = const0
             cusp_ion(:)  = const0
             drizl_ion(:) = const0
          endwhere
!-----------------------------------------------------------------------
! 	... auroral electrons
!-----------------------------------------------------------------------
          call aion( xalfa1, alfa1_ion, do_aurora, ncol )
          call aion( xalfa2, alfa2_ion, do_aurora, ncol )
          call aion( xcusp , cusp_ion, do_aurora, ncol  )
          call aion( xdrizl, drizl_ion, do_aurora, ncol )
          where( do_aurora(:) )
             falfa1(:) = alfa1(:)*flux1(:)  ! s7
             falfa2(:) = alfa2(:)*flux2(:)  ! s8
             fcusp (:) = cusp(:)*alfac*fc   ! s9
             qsum(:)   = falfa1(:)*alfa1_ion(:) &    ! s7*s3
                       + falfa2(:)*alfa2_ion(:) &    ! s8*s4
                       + fcusp(:)*cusp_ion (:) &     ! s9*s5
                       + drizl(:)*drizl_ion(:)       ! s10*s6
          endwhere

!-----------------------------------------------------------------------
! 	... form production
!-----------------------------------------------------------------------
          where( do_aurora(:) )
             barm_t(:)   = grav*barm(:)/(35.e-3_r8*gask*tempi(:))
             tpions(:,k) = qsum(:)*barm_t(:)               ! s1 = s1*s11
          endwhere
      end do level_loop

      end subroutine total_ion_prod

      subroutine aion( si, so, do_aurora, ncol )
!-----------------------------------------------------------------------
! Calculates integrated f(x) needed for total auroral ionization.
! See equations (10-12) in Roble,1987.
! Coefficients for equation (12) of Roble,1987 are in variable cc 
! (revised since 1987):
! Uses the identity x**y = exp(y*ln(x)) for performance 
! (fewer (1/2) trancendental functions are required).
!------------------------------------------------------------------------

      implicit none

!------------------------------------------------------------------------
! 	... dummy arguments
!------------------------------------------------------------------------
      integer,  intent(in)  :: ncol
      real(r8), intent(in)  :: si(ncol)
      real(r8), intent(out) :: so(ncol)
      logical,  intent(in)  :: do_aurora(ncol)

!------------------------------------------------------------------------
! 	... local variables
!------------------------------------------------------------------------
      real(r8), parameter :: cc(8) = &
       (/ 3.2333134511131_r8 ,  2.5658873458085_r8 ,  2.2540957232641_r8 , &
          0.72971983372673_r8,  1.1069072431948_r8 ,  1.7134937681128_r8 , &
          1.8835442312993_r8 ,  0.86472135072090_r8 /)

      real(r8) :: xlog(ncol)

      where( do_aurora(:) )
         xlog(:) = log( si(:) )
         so(:)   = cc(1)*exp( cc(2)*xlog(:) - cc(3)*exp( cc(4)*xlog(:) ) ) &
                   + cc(5)*exp( cc(6)*xlog(:) - cc(7)*exp( cc(8)*xlog(:) ) )
      elsewhere
         so(:) = 0._r8
      endwhere

      end subroutine aion

      end module mo_aurora
