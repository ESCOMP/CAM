!-----------------------------------------------------------------------
! Reads namelist options for slh recycling efficiency using the Free-Regime Aprox. (FRA)
!
! Created by Rafa Fernandez -- 30 June 2023
!-----------------------------------------------------------------------
module mo_slh_routines

  use spmd_utils,     only : masterproc
  use cam_abortutils, only : endrun
  use shr_kind_mod,   only : r8 => shr_kind_r8
  use constituents,   only : cnst_mw, cnst_get_ind
  use chem_mods,      only : gas_pcnst
  use ppgrid,         only : pver
  use cam_history,    only : addfld,add_default,horiz_only, outfld
  use infnan,         only : nan, assignment(=)
  use physconst,      only : mwdry ! molecular weight of dry air
  use physics_types,  only : physics_state
  use camsrfexch,     only : cam_in_t
  use cam_logfile,    only : iulog

  implicit none

  private

  real(r8), public, protected :: SSAdehal_ScalingFactor  = 1._r8
  real(r8), public, protected :: SSAhno3_ScalingFactor   = 1._r8
  real(r8), public, protected :: SSAn2o5_ScalingFactor   = 1._r8
  real(r8), public, protected :: LIQfraprx_ScalingFactor_I   = 1._r8
  real(r8), public, protected :: ICEfraprx_ScalingFactor_I   = 1._r8
  real(r8), public, protected :: ICEfraprx_ScalingFactor_Br  = 1._r8

  public :: iodine_emissions_init
  public :: iodine_emissions_srf
  public :: slh_readnl

  integer  :: hoi_ndx=-1, i2_ndx=-1
  integer  :: hoi_cnst=-1, i2_cnst=-1, o3_cnst=-1

  real(r8) :: conv_fact_i2  = -huge(1._r8)
  real(r8) :: conv_fact_hoi = -huge(1._r8)
  real(r8) :: conv_fact_o3  = -huge(1._r8)

  logical :: do_iodine_emis = .false.
  logical :: has_emis_i2    = .false.
  logical :: has_emis_hoi   = .false.

contains

  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  subroutine slh_readnl(nlfile)

    use namelist_utils, only: find_group_name
    use spmd_utils,     only: mpicom, masterprocid, mpi_real8, mpi_success

    character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

    integer :: unitn, ierr
    character(len=*), parameter :: subname = 'slh_readnl'

    namelist /slh_nl/ SSAdehal_ScalingFactor, SSAhno3_ScalingFactor, SSAn2o5_ScalingFactor, &
                      LIQfraprx_ScalingFactor_I, ICEfraprx_ScalingFactor_I, ICEfraprx_ScalingFactor_Br

    if (masterproc) then
       open( newunit=unitn, file=trim(nlfile), status='old' )
       call find_group_name(unitn, 'slh_nl', status=ierr)
       if (ierr == 0) then
          read(unitn, slh_nl, iostat=ierr)
          if (ierr /= 0) then
             call endrun('mo_slh_routines->slh_readnl: ERROR reading slh_nl namelist')
          end if
       end if
       close(unitn)
    end if

    call mpi_bcast(SSAdehal_ScalingFactor, 1, mpi_real8, masterprocid, mpicom, ierr)
    if (ierr/=mpi_success) then
       call endrun(subname//': MPI_BCAST ERROR: SSAdehal_ScalingFactor')
    end if
    call mpi_bcast(SSAhno3_ScalingFactor, 1, mpi_real8, masterprocid, mpicom, ierr)
    if (ierr/=mpi_success) then
       call endrun(subname//': MPI_BCAST ERROR: SSAhno3_ScalingFactor')
    end if
    call mpi_bcast(SSAn2o5_ScalingFactor, 1, mpi_real8, masterprocid, mpicom, ierr)
    if (ierr/=mpi_success) then
       call endrun(subname//': MPI_BCAST ERROR: SSAn2o5_ScalingFactor')
    end if
    call mpi_bcast(LIQfraprx_ScalingFactor_I, 1, mpi_real8, masterprocid, mpicom, ierr)
    if (ierr/=mpi_success) then
       call endrun(subname//': MPI_BCAST ERROR: LIQfraprx_ScalingFactor_I')
    end if
    call mpi_bcast(ICEfraprx_ScalingFactor_I, 1, mpi_real8, masterprocid, mpicom, ierr)
    if (ierr/=mpi_success) then
       call endrun(subname//': MPI_BCAST ERROR: ICEfraprx_ScalingFactor_I')
    end if
    call mpi_bcast(ICEfraprx_ScalingFactor_Br, 1, mpi_real8, masterprocid, mpicom, ierr)
    if (ierr/=mpi_success) then
       call endrun(subname//': MPI_BCAST ERROR: ICEfraprx_ScalingFactor_Br')
    end if

    if (masterproc) then
       write(iulog,*) subname//': SSAdehal_ScalingFactor = ',SSAdehal_ScalingFactor
       write(iulog,*) subname//': SSAhno3_ScalingFactor = ',SSAhno3_ScalingFactor
       write(iulog,*) subname//': SSAn2o5_ScalingFactor = ',SSAn2o5_ScalingFactor
       write(iulog,*) subname//': LIQfraprx_ScalingFactor_I = ',LIQfraprx_ScalingFactor_I
       write(iulog,*) subname//': ICEfraprx_ScalingFactor_I = ',ICEfraprx_ScalingFactor_I
       write(iulog,*) subname//': ICEfraprx_ScalingFactor_Br = ',ICEfraprx_ScalingFactor_Br
    end if

  end subroutine slh_readnl

  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  subroutine iodine_emissions_init( srf_emis_specifier )

    use mo_chem_utls,  only : get_spc_ndx

    !-----------------------------------------------------------------------
    ! 	... dummy arguments
    !-----------------------------------------------------------------------
    character(len=*), intent(in) :: srf_emis_specifier(:)

    !-----------------------------------------------------------------------
    ! 	... local variables
    !-----------------------------------------------------------------------
    integer  :: m, n, i                     ! Indices
    character(len=16)  :: spc_name
    real(r8) :: i2_mw, hoi_mw, o3_mw    ! Molecular weight of species

    call cnst_get_ind( 'I2',  i2_cnst,  abort=.false.  )
    call cnst_get_ind( 'HOI', hoi_cnst, abort=.false.  )
    i2_ndx      = get_spc_ndx('I2'     )
    hoi_ndx     = get_spc_ndx('HOI'    )

    has_emis_i2  = .false.
    has_emis_hoi = .false.

    count_emis: do n=1,gas_pcnst
       if ( len_trim(srf_emis_specifier(n) ) == 0 ) then
          exit count_emis
       endif
       i = scan(srf_emis_specifier(n),'->')
       spc_name = trim(adjustl(srf_emis_specifier(n)(:i-1)))

       m = get_spc_ndx(spc_name)

       if (m == i2_ndx) then
          has_emis_i2  = .true.
       endif
       if (m == hoi_ndx) then
          has_emis_hoi = .true.
       endif
    enddo count_emis

!   it is a requirement to include I2 and HOI emissions in user_nl_cam, not only to
!   have defined iodine species in chem_mech.in
    do_iodine_emis = has_emis_i2 .and. has_emis_hoi .and. (i2_ndx>0) .and. (hoi_ndx>0)

    if (.not. do_iodine_emis) return

    call cnst_get_ind( 'O3', o3_cnst )

    o3_mw = cnst_mw(o3_cnst)
    conv_fact_o3 = 1.0e9_r8 * mwdry/o3_mw

    i2_mw  = cnst_mw(i2_cnst)
    hoi_mw = cnst_mw(hoi_cnst)

    conv_fact_i2  = (1._r8 / 86400._r8) * 1.0e-9_r8 * i2_mw  * 1.e-3_r8    ! converts units of nmol/m2/day to kg/m2/s
    conv_fact_hoi = (1._r8 / 86400._r8) * 1.0e-9_r8 * hoi_mw * 1.e-3_r8    ! converts units of nmol/m2/day to kg/m2/s

    call addfld( 'IODIDE_AQ',   horiz_only,  'I', 'mol/dm3', 'Oceanic Iodide concentration (formula)'  )
    call addfld( 'WS_msk',      horiz_only,  'I', 'm/s',     'Wind Surface Speed with mask' )
    call addfld( 'O3_srf',      horiz_only,  'I', 'mol/mol', 'O3 surface volume mixing ratio' )
    call addfld( 'FLX_I2',      horiz_only,  'I', 'nmol/m2/day',  'Total Iodine Flux from Ocean (formula)' )
    call addfld( 'FLX_HOI',     horiz_only,  'I', 'nmol/m2/day',  'Total Iodine Flux from Ocean (formula)' )
    call addfld( 'FLX_ITOT',    horiz_only,  'I', 'nmol/m2/day',  'Total Iodine Flux from Ocean (formula)' )

  end subroutine iodine_emissions_init

  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  subroutine iodine_emissions_srf( state, cam_in )

    type(physics_state),    intent(in)    :: state   ! Physics state variables
    type(cam_in_t),         intent(inout) :: cam_in  ! import state

    integer :: lchnk, ncol
    real(r8) :: o3_surf(state%ncol)                  ! surface ozone vmr
    real(r8) :: ws_mask(state%ncol)                  ! Wind Speed with mask for ws > 3.0 m/s
    real(r8) :: iodide_aq(state%ncol)                ! Iodide disolved in ocean surface waters (mol dm-3)
    real(r8) :: flx_i2(state%ncol)                   ! Flux of I2  form formula in paper
    real(r8) :: flx_hoi(state%ncol)                  ! Flux of I2  form formula in paper
    real(r8) :: flx_itot(state%ncol)                 ! Flux of I2  form formula in paper

    if (.not. do_iodine_emis) return

    lchnk= state%lchnk
    ncol = state%ncol
    o3_surf (:ncol) = state%q(:ncol,pver,o3_cnst) * conv_fact_o3
    ws_mask(:ncol) = sqrt( state%u(:ncol,pver)*state%u(:ncol,pver) &
                         + state%v(:ncol,pver)*state%v(:ncol,pver) )

!rpf and cac: This threshold is required to aviod huge iodine emissions for low wind speeds
    where ( ws_mask(:ncol) < 3._r8 )
      ws_mask(:ncol) = 3._r8
    endwhere

    call iodine_emissions_set( ncol, cam_in%ocnfrac(:ncol), o3_surf, cam_in%sst(:ncol), ws_mask,  &
                              iodide_aq, flx_i2, flx_hoi )

    flx_itot(:)    = 2._r8 * flx_i2(:) + flx_hoi(:)

    call outfld( 'IODIDE_AQ',  iodide_aq (:ncol), ncol, lchnk )
    call outfld( 'WS_msk',     ws_mask   (:ncol), ncol, lchnk )
    call outfld( 'O3_srf',     o3_surf   (:ncol), ncol, lchnk )

    call outfld( 'FLX_I2',     flx_i2    (:ncol), ncol, lchnk )
    call outfld( 'FLX_HOI',    flx_hoi   (:ncol), ncol, lchnk )
    call outfld( 'FLX_ITOT',   flx_itot  (:ncol), ncol, lchnk )

    cam_in%cflx(:ncol,i2_cnst) = cam_in%cflx(:ncol,i2_cnst) + flx_i2(:ncol) * conv_fact_i2
    cam_in%cflx(:ncol,hoi_cnst) = cam_in%cflx(:ncol,hoi_cnst) + flx_hoi(:ncol) * conv_fact_hoi

  end subroutine iodine_emissions_srf

  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  subroutine iodine_emissions_set( ncol, ocnfrac, o3_surf, sst, ws_mask, iodide_aq, flx_i2, flx_hoi )

    integer,  intent(in) :: ncol
    real(r8), intent(in)    :: ocnfrac(ncol)         ! ocean fraction
    real(r8), intent(in)    :: o3_surf(ncol)         ! surface ozone vmr
    real(r8), intent(in)    :: sst(ncol)             ! SST
    real(r8), intent(in)    :: ws_mask(ncol)         ! Wind Speed with mask for ws > 3 m/s
    real(r8), intent(out)   :: iodide_aq(ncol)       ! Iodide disolved in ocean surface waters (mol dm-3)
    real(r8), intent(out)   :: flx_i2(ncol)          ! Oceanic flux of hoi (nmol/m2/day)
    real(r8), intent(out)   :: flx_hoi(ncol)         ! Oceanic flux of i2  (nmol/m2/day)

    integer :: i

    if (.not. do_iodine_emis) return

    flx_i2(:) = nan
    flx_hoi(:) = nan

    do i = 1,ncol
       !!! Be CAREFULL that if SST = 0 then the division is not possible. Should it be forced to zero?
       if (ocnfrac(i) == 0) then
          ! Actually it should be landfrac = 1 ... but we do not want to bring landfrac as an additional agrument
          iodide_aq(i) = 0._r8
       else
          iodide_aq(i) = 1.46e6_r8 * exp( -9134._r8 / sst(i) )
       end if

       if (ws_mask(i) > 0.0_r8) then
          ! rpf oceanic iodine emissions were extremely large for very small windspeeds.
          ! ws_mask = 3.0 m/s for low wind-speeds imposed above
          flx_i2(i) = o3_surf(i) * ( iodide_aq(i) ) ** (1.3_r8) * ( 1.74e9_r8 - ( 6.54e8_r8 * log ( ws_mask(i) ) ) )
          flx_i2(i) = flx_i2(i) * ocnfrac(i)
       else
          flx_i2(i) = 0._r8
       endif

       if ( flx_i2(i) < 0._r8 ) flx_i2(i) = 0._r8

       if (ws_mask(i) > 0.0_r8) then
          ! rpf oceanic iodine emissions were extremely large for very small windspeeds.
          ! ws_mask = 3.0 m/s for low wind-speeds imposed above
          flx_hoi(i) = o3_surf(i) * ( 4.15e5_r8 * ( sqrt( iodide_aq(i) ) / ws_mask(i) ) &
               - ( 20.6_r8 / ws_mask(i) ) - 23600._r8 * sqrt ( iodide_aq(i) ) )
          flx_hoi(i) = flx_hoi(i) * ocnfrac(i)
       else
          flx_hoi(i) = 0._r8
       endif

       if ( flx_hoi(i) < 0._r8 ) flx_hoi(i) = 0._r8    ! for low Iodide_aq values it gives negative fluxes
    end do

  end subroutine iodine_emissions_set

end module mo_slh_routines
