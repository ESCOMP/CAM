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

   implicit none

PRIVATE

   public :: scm_intht

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
      use cam_history, only: addfld, add_default, horiz_only
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
! Local variables
!
      integer m,j        ! Indices
      real(r8) dummy
!
! Call addfld to add each field to the Master Field List.
!
      call addfld ('TDIFF',    (/ 'lev' /), 'A', 'K','difference from observed temp',                    gridname='gauss_grid')
      call addfld ('UDIFF',    (/ 'lev' /), 'A', 'K','difference from observed u wind',                  gridname='gauss_grid')
      call addfld ('VDIFF',    (/ 'lev' /), 'A', 'K','difference from observed v wind',                  gridname='gauss_grid')

      call addfld ('TOBS',     (/ 'lev' /), 'A', 'K','observed temp')
      call addfld ('QDIFF',    (/ 'lev' /), 'A', 'kg/kg','difference from observed water',               gridname='gauss_grid')

      call addfld ('QOBS',     (/ 'lev' /), 'A', 'kg/kg','observed water',                               gridname='physgrid')
      call addfld ('PRECOBS',  (/ 'lev' /), 'A', 'mm/day','Total (convective and large-scale) precipitation rate',             &
                                                                                                         gridname='physgrid')
      call addfld ('DIVQ',     (/ 'lev' /), 'A', 'kg/kg/s','Q advection tendency (horizontal)',          gridname='physgrid')
      call addfld ('DIVQ3D',   (/ 'lev' /), 'A', 'kg/kg/s','Q advection tendency (horiz/vert combined)', gridname='gauss_grid')
      call addfld ('DIVV',     (/ 'lev' /), 'A', 'm/s2','V advection tendency (horizontal)',             gridname='physgrid')
      call addfld ('DIVU',     (/ 'lev' /), 'A', 'm/s2','U advection tendency (horizontal)',             gridname='physgrid')
      call addfld ('DIVT',     (/ 'lev' /), 'A', 'K/s','T advection tendency (horizontal)',              gridname='physgrid')
      call addfld ('DIVT3D',   (/ 'lev' /), 'A', 'K/s','T advection tendency (horiz/vert combined)',     gridname='gauss_grid')
      call addfld ('DIVU3D',   (/ 'lev' /), 'A', 'K/s','U advection tendency (horiz/vert combined)',       gridname='gauss_grid')
      call addfld ('DIVV3D',   (/ 'lev' /), 'A', 'K/s','V advection tendency (horiz/vert combined)',       gridname='gauss_grid')

      call addfld ('SHFLXOBS', horiz_only,  'A', 'W/m2','Obs Surface sensible heat flux',                gridname='physgrid')
      call addfld ('LHFLXOBS', horiz_only,  'A', 'W/m2','Obs Surface latent heat flux',                  gridname='physgrid')
      call addfld ('TRELAX',   (/ 'lev' /), 'A', 'K','t relaxation amount',                              gridname='gauss_grid')
      call addfld ('QRELAX',   (/ 'lev' /), 'A', 'kg/kg','q relaxation amount',                          gridname='gauss_grid')
      call addfld ('TAURELAX', (/ 'lev' /), 'A', 'seconds','relaxation time constant',                   gridname='gauss_grid')
      call add_default ('TDIFF', 1, ' ')
      call add_default ('QDIFF', 1, ' ')

    ! Vertical advective forcing of 'T,u,v,qv,ql,qi,nl,ni' in forecast.F90

      call addfld ('TTEN_XYADV',   (/ 'lev' /), 'I', 'K/s',    'T  horizontal advective forcing',            gridname='gauss_grid' )
      call addfld ('UTEN_XYADV',   (/ 'lev' /), 'I', 'm/s^2',  'U  horizontal advective forcing',            gridname='gauss_grid' )
      call addfld ('VTEN_XYADV',   (/ 'lev' /), 'I', 'm/s^2',  'V  horizontal advective forcing',            gridname='gauss_grid' )
      call addfld ('QVTEN_XYADV',  (/ 'lev' /), 'I', 'kg/kg/s','QV horizontal advective forcing',            gridname='gauss_grid' )
      call addfld ('QLTEN_XYADV',  (/ 'lev' /), 'I', 'kg/kg/s','QL horizontal advective forcing',            gridname='gauss_grid' )
      call addfld ('QITEN_XYADV',  (/ 'lev' /), 'I', 'kg/kg/s','QI horizontal advective forcing',            gridname='gauss_grid' )
      call addfld ('NLTEN_XYADV',  (/ 'lev' /), 'I', '#/kg/s', 'NL horizontal advective forcing',            gridname='gauss_grid' )
      call addfld ('NITEN_XYADV',  (/ 'lev' /), 'I', '#/kg/s', 'NI horizontal advective forcing',            gridname='gauss_grid' )

!      call addfld ('T3D_ADV_SLT', 'K/s'     , pver, 'I', 'T  3d slt advective forcing',                 gridname='physgrid')
!      call addfld ('U3D_ADV_SLT', 'm/s^2'   , pver, 'I', 'U  3d slt advective forcing',                 gridname='physgrid')
!      call addfld ('V3D_ADV_SLT', 'm/s^2'   , pver, 'I', 'V  3d slt advective forcing',                 gridname='physgrid')
      call addfld ('TTEN_ZADV',   (/ 'lev' /), 'I', 'K/s',    'T  vertical   advective forcing',            gridname='gauss_grid' )
      call addfld ('UTEN_ZADV',   (/ 'lev' /), 'I', 'm/s^2',  'U  vertical   advective forcing',            gridname='gauss_grid' )
      call addfld ('VTEN_ZADV',   (/ 'lev' /), 'I', 'm/s^2',  'V  vertical   advective forcing',            gridname='gauss_grid' )
      call addfld ('QVTEN_ZADV',  (/ 'lev' /), 'I', 'kg/kg/s','QV vertical   advective forcing',            gridname='gauss_grid' )
      call addfld ('QLTEN_ZADV',  (/ 'lev' /), 'I', 'kg/kg/s','QL vertical   advective forcing',            gridname='gauss_grid' )
      call addfld ('QITEN_ZADV',  (/ 'lev' /), 'I', 'kg/kg/s','QI vertical   advective forcing',            gridname='gauss_grid' )
      call addfld ('NLTEN_ZADV',  (/ 'lev' /), 'I', '#/kg/s', 'NL vertical   advective forcing',            gridname='gauss_grid' )
      call addfld ('NITEN_ZADV',  (/ 'lev' /), 'I', '#/kg/s', 'NI vertical   advective forcing',            gridname='gauss_grid' )

      call addfld ('TTEN_PHYS',   (/ 'lev' /), 'I', 'K/s',    'T  vertical   advective forcing',            gridname='gauss_grid' )
      call addfld ('UTEN_PHYS',   (/ 'lev' /), 'I', 'm/s^2',  'U  vertical   advective forcing',            gridname='gauss_grid' )
      call addfld ('VTEN_PHYS',   (/ 'lev' /), 'I', 'm/s^2',  'V  vertical   advective forcing',            gridname='gauss_grid' )
      call addfld ('QVTEN_PHYS',   (/ 'lev' /), 'I','kg/kg/s','QV vertical   advective forcing',            gridname='gauss_grid' )
      call addfld ('QLTEN_PHYS',   (/ 'lev' /), 'I','kg/kg/s','QL vertical   advective forcing',            gridname='gauss_grid' )
      call addfld ('QITEN_PHYS',   (/ 'lev' /), 'I','kg/kg/s','QI vertical   advective forcing',            gridname='gauss_grid' )
      call addfld ('NLTEN_PHYS',   (/ 'lev' /), 'I','#/kg/s', 'NL vertical   advective forcing',            gridname='gauss_grid' )
      call addfld ('NITEN_PHYS',   (/ 'lev' /), 'I','#/kg/s', 'NI vertical   advective forcing',            gridname='gauss_grid' )

   end subroutine scm_intht

!#######################################################################
 end module history_scam
