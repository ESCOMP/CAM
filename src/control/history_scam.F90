module history_scam
!-----------------------------------------------------------------------
!
! Purpose: SCAM specific history code.
!
! Public functions/subroutines:
!   bldfld, h_default
!
! Author: anonymous from code in cam_history.F90
!-----------------------------------------------------------------------
   use shr_kind_mod, only: r8 => shr_kind_r8
   use cam_history, only: addfld, add_default, horiz_only

   implicit none

PRIVATE

   public :: scm_intht
   public :: initialize_iop_history

!#######################################################################
CONTAINS
   subroutine scm_intht()
!-----------------------------------------------------------------------
!
! Purpose:
!
! add master list fields to scm
!
! Method: Call a subroutine to add each field
!
! Author: CCM Core Group
!
!-----------------------------------------------------------------------
      use dycore,      only: dycore_is
      use cam_history, only: write_camiop
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
! Local variables
!
      character(len=100) outgrid

      if (dycore_is('SE')) then
         ! for camiop mode use the GLL grid otherwise use physics grids for SCM mode output
         if (write_camiop) then
            outgrid = 'GLL'
         else
            outgrid = 'physgrid'
         end if
      else if (dycore_is('EUL')) then
         outgrid = 'gauss_grid'
      else
         outgrid = 'unknown'
      end if
!
! Call addfld to add each field to the Master Field List.
!
      call addfld ('TDIFF',    (/ 'lev' /), 'A', 'K','difference from observed temp',                    gridname=trim(outgrid))
      call addfld ('UDIFF',    (/ 'lev' /), 'A', 'K','difference from observed u wind',                  gridname=trim(outgrid))
      call addfld ('VDIFF',    (/ 'lev' /), 'A', 'K','difference from observed v wind',                  gridname=trim(outgrid))

      call addfld ('TOBS',     (/ 'lev' /), 'A', 'K','observed temp')
      call addfld ('QDIFF',    (/ 'lev' /), 'A', 'kg/kg','difference from observed water',               gridname=trim(outgrid))

      call addfld ('QOBS',     (/ 'lev' /), 'A', 'kg/kg','observed water',                               gridname='physgrid')
      call addfld ('PRECOBS',  (/ 'lev' /), 'A', 'mm/day','Total (convective and large-scale) precipitation rate',             &
                                                                                                         gridname='physgrid')
      call addfld ('DIVQ',     (/ 'lev' /), 'A', 'kg/kg/s','Q advection tendency (horizontal)',          gridname='physgrid')
      call addfld ('DIVQ3D',   (/ 'lev' /), 'A', 'kg/kg/s','Q advection tendency (horiz/vert combined)', gridname=trim(outgrid))
      call addfld ('DIVV',     (/ 'lev' /), 'A', 'm/s2','V advection tendency (horizontal)',             gridname='physgrid')
      call addfld ('DIVU',     (/ 'lev' /), 'A', 'm/s2','U advection tendency (horizontal)',             gridname='physgrid')
      call addfld ('DIVT',     (/ 'lev' /), 'A', 'K/s','T advection tendency (horizontal)',              gridname='physgrid')
      call addfld ('DIVT3D',   (/ 'lev' /), 'A', 'K/s','T advection tendency (horiz/vert combined)',     gridname=trim(outgrid))
      call addfld ('DIVU3D',   (/ 'lev' /), 'A', 'K/s','U advection tendency (horiz/vert combined)',       gridname=trim(outgrid))
      call addfld ('DIVV3D',   (/ 'lev' /), 'A', 'K/s','V advection tendency (horiz/vert combined)',       gridname=trim(outgrid))

      call addfld ('SHFLXOBS', horiz_only,  'A', 'W/m2','Obs Surface sensible heat flux',                gridname='physgrid')
      call addfld ('LHFLXOBS', horiz_only,  'A', 'W/m2','Obs Surface latent heat flux',                  gridname='physgrid')
      call addfld ('TRELAX',   (/ 'lev' /), 'A', 'K','t relaxation amount',                              gridname=trim(outgrid))
      call addfld ('QRELAX',   (/ 'lev' /), 'A', 'kg/kg','q relaxation amount',                          gridname=trim(outgrid))
      call addfld ('TAURELAX', (/ 'lev' /), 'A', 'seconds','relaxation time constant',                   gridname=trim(outgrid))
      call add_default ('TDIFF', 1, ' ')
      call add_default ('QDIFF', 1, ' ')

    ! Vertical advective forcing of 'T,u,v,qv,ql,qi,nl,ni' in forecast.F90

      call addfld ('TTEN_XYADV',   (/ 'lev' /), 'I', 'K/s',    'T  horizontal advective forcing',            gridname=trim(outgrid) )
      call addfld ('UTEN_XYADV',   (/ 'lev' /), 'I', 'm/s^2',  'U  horizontal advective forcing',            gridname=trim(outgrid) )
      call addfld ('VTEN_XYADV',   (/ 'lev' /), 'I', 'm/s^2',  'V  horizontal advective forcing',            gridname=trim(outgrid) )
      call addfld ('QVTEN_XYADV',  (/ 'lev' /), 'I', 'kg/kg/s','QV horizontal advective forcing',            gridname=trim(outgrid) )
      call addfld ('QLTEN_XYADV',  (/ 'lev' /), 'I', 'kg/kg/s','QL horizontal advective forcing',            gridname=trim(outgrid) )
      call addfld ('QITEN_XYADV',  (/ 'lev' /), 'I', 'kg/kg/s','QI horizontal advective forcing',            gridname=trim(outgrid) )
      call addfld ('NLTEN_XYADV',  (/ 'lev' /), 'I', '#/kg/s', 'NL horizontal advective forcing',            gridname=trim(outgrid) )
      call addfld ('NITEN_XYADV',  (/ 'lev' /), 'I', '#/kg/s', 'NI horizontal advective forcing',            gridname=trim(outgrid) )

!      call addfld ('T3D_ADV_SLT', 'K/s'     , pver, 'I', 'T  3d slt advective forcing',                 gridname='physgrid')
!      call addfld ('U3D_ADV_SLT', 'm/s^2'   , pver, 'I', 'U  3d slt advective forcing',                 gridname='physgrid')
!      call addfld ('V3D_ADV_SLT', 'm/s^2'   , pver, 'I', 'V  3d slt advective forcing',                 gridname='physgrid')
      call addfld ('TTEN_ZADV',   (/ 'lev' /), 'I', 'K/s',    'T  vertical   advective forcing',            gridname=trim(outgrid) )
      call addfld ('UTEN_ZADV',   (/ 'lev' /), 'I', 'm/s^2',  'U  vertical   advective forcing',            gridname=trim(outgrid) )
      call addfld ('VTEN_ZADV',   (/ 'lev' /), 'I', 'm/s^2',  'V  vertical   advective forcing',            gridname=trim(outgrid) )
      call addfld ('QVTEN_ZADV',  (/ 'lev' /), 'I', 'kg/kg/s','QV vertical   advective forcing',            gridname=trim(outgrid) )
      call addfld ('QLTEN_ZADV',  (/ 'lev' /), 'I', 'kg/kg/s','QL vertical   advective forcing',            gridname=trim(outgrid) )
      call addfld ('QITEN_ZADV',  (/ 'lev' /), 'I', 'kg/kg/s','QI vertical   advective forcing',            gridname=trim(outgrid) )
      call addfld ('NLTEN_ZADV',  (/ 'lev' /), 'I', '#/kg/s', 'NL vertical   advective forcing',            gridname=trim(outgrid) )
      call addfld ('NITEN_ZADV',  (/ 'lev' /), 'I', '#/kg/s', 'NI vertical   advective forcing',            gridname=trim(outgrid) )

      call addfld ('TTEN_PHYS',   (/ 'lev' /), 'I', 'K/s',    'T  vertical   advective forcing',            gridname=trim(outgrid) )
      call addfld ('UTEN_PHYS',   (/ 'lev' /), 'I', 'm/s^2',  'U  vertical   advective forcing',            gridname=trim(outgrid) )
      call addfld ('VTEN_PHYS',   (/ 'lev' /), 'I', 'm/s^2',  'V  vertical   advective forcing',            gridname=trim(outgrid) )
      call addfld ('QVTEN_PHYS',   (/ 'lev' /), 'I','kg/kg/s','QV vertical   advective forcing',            gridname=trim(outgrid) )
      call addfld ('QLTEN_PHYS',   (/ 'lev' /), 'I','kg/kg/s','QL vertical   advective forcing',            gridname=trim(outgrid) )
      call addfld ('QITEN_PHYS',   (/ 'lev' /), 'I','kg/kg/s','QI vertical   advective forcing',            gridname=trim(outgrid) )
      call addfld ('NLTEN_PHYS',   (/ 'lev' /), 'I','#/kg/s', 'NL vertical   advective forcing',            gridname=trim(outgrid) )
      call addfld ('NITEN_PHYS',   (/ 'lev' /), 'I','#/kg/s', 'NI vertical   advective forcing',            gridname=trim(outgrid) )

   end subroutine scm_intht
!#######################################################################
   subroutine initialize_iop_history()
!-----------------------------------------------------------------------
!
! Purpose: Add fields and set defaults for SCAM CAM BFB IOP initial file
! as well as single column output history
!
! Method: Call a subroutine to add each field
!
!-----------------------------------------------------------------------
!
! !USES:
    use constituents,     only: pcnst, cnst_name
    use dycore,           only: dycore_is
! !ARGUMENTS:
    implicit none

! !LOCAL VARIABLES:
    integer m
    character(len=100) outgrid

!-----------------------------------------------------------------------

    if (dycore_is('SE')) then
       outgrid = 'GLL'
    else if (dycore_is('EUL')) then
       outgrid = 'gauss_grid'
    else if (dycore_is('EUL')) then
       outgrid = 'unknown'
    end if

    if (trim(outgrid) == 'gauss_grid') then
       call addfld ('CLAT1&IC',  horiz_only,  'I', ' ','cos lat for bfb testing', gridname=trim(outgrid))
       call add_default ('CLAT1&IC',0,'I')
       call addfld ('CLON1&IC',  horiz_only,  'I', ' ','cos lon for bfb testing', gridname=trim(outgrid))
       call add_default ('CLON1&IC',0,'I')
       call addfld ('PHI&IC',    horiz_only,  'I', ' ','lat for bfb testing', gridname=trim(outgrid))
       call add_default ('PHI&IC',0,  'I')
       call addfld ('LAM&IC',    horiz_only,  'I', ' ','lon for bfb testing', gridname=trim(outgrid))
       call add_default ('LAM&IC',0,  'I')

       call addfld ('CLAT',    horiz_only,   'A', ' ',   'cos lat for bfb testing', gridname=trim(outgrid))
       call add_default ('CLAT',2,' ')

       call addfld ('fixmas',  horiz_only,   'A', 'percent','Mass fixer',gridname=trim(outgrid))
       call add_default ('fixmas',2,' ')
       call addfld ('beta',    horiz_only,   'A', 'percent','Mass fixer',gridname=trim(outgrid))
       call add_default ('beta',2,' ')
    end if

    call addfld ('q',       (/ 'lev' /),  'A', 'kg/kg',  'Q for scam',gridname=trim(outgrid))
    call add_default ('q',2, ' ')
    call addfld ('u',       (/ 'lev' /),  'A', 'm/s',    'U for scam',gridname=trim(outgrid))
    call add_default ('u',2,' ')
    call addfld ('v',       (/ 'lev' /),  'A', 'm/s',    'V for scam',gridname=trim(outgrid))
    call add_default ('v',2,' ')
    call addfld ('t',       (/ 'lev' /),  'A', 'K',      'Temperature for scam',gridname=trim(outgrid))
    call add_default ('t',2,' ')
    call addfld ('Tg',      horiz_only,   'A', 'K',      'Surface temperature (radiative) for scam',gridname='physgrid')
    call add_default ('Tg',2,' ')
    call addfld ('Ps',      horiz_only,   'A', 'Pa',     'Ps for scam',gridname=trim(outgrid))
    call add_default ('Ps',2,' ')
    call addfld ('divT3d',  (/ 'lev' /),  'A', 'K',      'Dynamics Residual for T',gridname=trim(outgrid))
    call add_default ('divT3d',2,' ')
    call addfld ('divU3d',  (/ 'lev' /),  'A', 'K',      'Dynamics Residual for U',gridname=trim(outgrid))
    call add_default ('divU3d',2,' ')
    call addfld ('divV3d',  (/ 'lev' /),  'A', 'K',      'Dynamics Residual for V',gridname=trim(outgrid))
    call add_default ('divV3d',2,' ')
    call addfld ('heat_glob',horiz_only, 'A', 'K/s', 'Global mean total energy difference')
    call add_default ('heat_glob',2,' ')
    do m=1,pcnst
       call addfld (trim(cnst_name(m))//'_dten', (/ 'lev' /), 'A', 'kg/kg', &
            trim(cnst_name(m))//' IOP Dynamics Residual for '//trim(cnst_name(m)),gridname=trim(outgrid))
       call add_default (trim(cnst_name(m))//'_dten',2,' ')
       if (trim(outgrid) == 'gauss_grid') then
          call addfld (trim(cnst_name(m))//'_alph', horiz_only, 'A', 'kg/kg',trim(cnst_name(m))//' alpha constituent fixer', &
               gridname=trim(outgrid))
          call add_default (trim(cnst_name(m))//'_alph',2,' ')
          call addfld (trim(cnst_name(m))//'_dqfx', (/ 'lev' /), 'A', 'kg/kg',trim(cnst_name(m))//' dqfx3 fixer',            &
               gridname=trim(outgrid))
          call add_default (trim(cnst_name(m))//'_dqfx',2,' ')
       end if
    end do
    call addfld ('shflx',  horiz_only,  'A', 'W/m2', 'Surface sensible heat flux for scam',gridname='physgrid')
    call add_default ('shflx',2,' ')
    call addfld ('lhflx',  horiz_only,  'A', 'W/m2', 'Surface latent heat flux for scam',gridname='physgrid')
    call add_default ('lhflx',2,' ')
    call addfld ('trefht', horiz_only,  'A', 'K',    'Reference height temperature',gridname='physgrid')
    call add_default ('trefht',2,' ')
    call addfld ('Tsair',  horiz_only,  'A', 'K',    'Reference height temperature for scam',gridname='physgrid')
    call add_default ('Tsair',2,' ')
    call addfld ('phis',   horiz_only,  'I', 'm2/s2','Surface geopotential for scam',gridname='physgrid')
    call add_default ('phis',2,' ')
    call addfld ('Prec',   horiz_only,  'A', 'm/s',  'Total (convective and large-scale) precipitation rate for scam',   &
         gridname='physgrid')
    call add_default ('Prec',2,' ')
    call addfld ('omega',  (/ 'lev' /), 'A', 'Pa/s', 'Vertical velocity (pressure)',gridname='physgrid')
    call add_default ('omega',2,' ')

  end subroutine initialize_iop_history

!#######################################################################
 end module history_scam
