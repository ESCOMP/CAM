module mo_setsox

  use shr_kind_mod, only : r8 => shr_kind_r8
  use cam_logfile,  only : iulog
  use physics_types,only : physics_state

  implicit none

  private
  public :: sox_inti, setsox
  public :: has_sox

  logical            ::  inv_o3
  integer            ::  id_msa

  integer :: id_so2, id_nh3, id_hno3, id_h2o2, id_o3, id_ho2
  integer :: id_so4, id_h2so4

  logical :: has_sox = .true.
  logical :: inv_so2, inv_nh3, inv_hno3, inv_h2o2, inv_ox, inv_nh4no3, inv_ho2

  logical :: cloud_borne = .false.

  ! Indices for species in the shared array of Henry's Law constant parameters
  integer :: heff_id_hno3, heff_id_so2, heff_id_nh3, heff_id_co2, heff_id_h2o2, heff_id_o3

contains

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
  subroutine sox_inti
    !-----------------------------------------------------------------------
    !	... initialize the hetero sox routine
    !-----------------------------------------------------------------------

    use mo_chem_utls, only : get_spc_ndx, get_inv_ndx
    use spmd_utils,   only : masterproc
    use phys_control, only : phys_getopts
    use carma_flags_mod, only : carma_do_cloudborne
    use sox_cldaero_mod, only : sox_cldaero_init

    logical :: modal_aerosols

    call phys_getopts( prog_modal_aero_out=modal_aerosols )
    cloud_borne = modal_aerosols .or. carma_do_cloudborne

    !-----------------------------------------------------------------
    !       ... get species indicies
    !-----------------------------------------------------------------

    if (cloud_borne) then
       id_h2so4 = get_spc_ndx( 'H2SO4' )
    else
       id_so4 = get_spc_ndx( 'SO4' )
    endif
    id_msa = get_spc_ndx( 'MSA' )

    inv_so2 = .false.
    id_so2 = get_inv_ndx( 'SO2' )
    inv_so2 = id_so2 > 0
    if ( .not. inv_so2 ) then
       id_so2 = get_spc_ndx( 'SO2' )
    endif

    inv_NH3 = .false.
    id_NH3 = get_inv_ndx( 'NH3' )
    inv_NH3 = id_NH3 > 0
    if ( .not. inv_NH3 ) then
       id_NH3 = get_spc_ndx( 'NH3' )
    endif

    inv_HNO3 = .false.
    id_HNO3 = get_inv_ndx( 'HNO3' )
    inv_HNO3 = id_hno3 > 0
    if ( .not. inv_HNO3 ) then
       id_HNO3 = get_spc_ndx( 'HNO3' )
    endif

    inv_H2O2 = .false.
    id_H2O2 = get_inv_ndx( 'H2O2' )
    inv_H2O2 = id_H2O2 > 0
    if ( .not. inv_H2O2 ) then
       id_H2O2 = get_spc_ndx( 'H2O2' )
    endif

    inv_HO2 = .false.
    id_HO2 = get_inv_ndx( 'HO2' )
    inv_HO2 = id_HO2 > 0
    if ( .not. inv_HO2 ) then
       id_HO2 = get_spc_ndx( 'HO2' )
    endif

    inv_o3 = get_inv_ndx( 'O3' ) > 0
    if (inv_o3) then
       id_o3 = get_inv_ndx( 'O3' )
    else
       id_o3 = get_spc_ndx( 'O3' )
    endif
    inv_ho2 = get_inv_ndx( 'HO2' ) > 0
    if (inv_ho2) then
       id_ho2 = get_inv_ndx( 'HO2' )
    else
       id_ho2 = get_spc_ndx( 'HO2' )
    endif

    has_sox = (id_so2>0) .and. (id_h2o2>0) .and. (id_o3>0) .and. (id_ho2>0)
    if (cloud_borne) then
       has_sox = has_sox .and. (id_h2so4>0)
    else
       has_sox = has_sox .and. (id_so4>0) .and. (id_nh3>0)
    endif

    ! Lookup Effective Henry's Law Constant parameters from the common
    ! data file read in the shared code.
    heff_id_hno3 = get_heff_index( 'HNO3' )
    heff_id_so2  = get_heff_index( 'SO2'  )
    heff_id_nh3  = get_heff_index( 'NH3'  )
    heff_id_co2  = get_heff_index( 'CO2'  )
    heff_id_h2o2 = get_heff_index( 'H2O2' )
    heff_id_o3   = get_heff_index( 'OX'   )

    has_sox = has_sox .and. (heff_id_hno3 > 0) .and. (heff_id_so2 > 0) &
               .and. (heff_id_nh3 > 0) .and. (heff_id_co2 > 0) &
               .and. (heff_id_h2o2 > 0) .and. (heff_id_o3 > 0)

    if (masterproc) then
       write(iulog,*) 'sox_inti: has_sox = ',has_sox
    endif

    if( has_sox ) then
       if (masterproc) then
          write(iulog,*) '-----------------------------------------'
          write(iulog,*) ' mo_setsox will do sox aerosols'
          write(iulog,*) '-----------------------------------------'
       endif
    else
       if (masterproc) then
          write(iulog,*) '-----------------------------------------'
          write(iulog,*) ' mo_setsox will not do sox aerosols'
          write(iulog,*) '-----------------------------------------'
       endif
       return
    end if

    call sox_cldaero_init()

  end subroutine sox_inti

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
  subroutine setsox( state, &
       pbuf,   &
       ncol,   &
       lchnk,  &
       loffset,&
       dtime,  &
       press,  &
       pdel,   &
       tfld,   &
       mbar,   &
       lwc,    &
       cldfrc, &
       cldnum, &
       invariants, &
       qcw,    &
       qin,    &
       xphlwc, &
       aqso4,  &
       aqh2so4,&
       aqso4_h2o2, &
       aqso4_o3,   &
       yph_in,  &
       aqso4_h2o2_3d, &
       aqso4_o3_3d &
       )

    !-----------------------------------------------------------------------
    !          ... Compute heterogeneous reactions of SOX
    !
    !       (0) using initial PH to calculate PH
    !           (a) HENRYs law constants
    !           (b) PARTIONING
    !           (c) PH values
    !
    !       (1) using new PH to repeat
    !           (a) HENRYs law constants
    !           (b) PARTIONING
    !           (c) REACTION rates
    !           (d) PREDICTION
    !
    ! NOTE: This routine assumes an Ideal Gas.
    !-----------------------------------------------------------------------
    !
    use physconst,    only : AVOGADRO_KMOL => avogad, &
                             BOLTZMANN => boltz, &
                             GAS_CONSTANT_KMOL => r_universal, &
                             MOLECULAR_WEIGHT_CO2_G_MOL => mwco2, &
                             MOLECULAR_WEIGHT_DRY_AIR_G_MOL => mwdry
    use ppgrid,       only : pcols, pver
    use chem_mods,    only : gas_pcnst, nfs
    use chem_mods,    only : adv_mass
    use physconst,    only : mwdry, gravit
    use mo_constants, only : pi
    use sox_cldaero_mod, only : sox_cldaero_update, sox_cldaero_create_obj, sox_cldaero_destroy_obj
    use cldaero_mod,     only : cldaero_conc_t
    use shr_drydep_mod,  only : dheff
    use physics_buffer,  only : physics_buffer_desc
    use rad_constituents, only : rad_cnst_get_gas

    !
    !-----------------------------------------------------------------------
    !      ... Dummy arguments
    !-----------------------------------------------------------------------
    type(physics_state),                intent(in)    :: state   ! Physics state variables
    type(physics_buffer_desc), pointer, intent(inout) :: pbuf(:) ! Physics buffer
    integer,          intent(in)    :: ncol              ! num of columns in chunk
    integer,          intent(in)    :: lchnk             ! chunk id
    integer,          intent(in)    :: loffset           ! offset of chem tracers in the advected tracers array
    real(r8),         intent(in)    :: dtime             ! time step (sec)
    real(r8),         intent(in)    :: press(:,:)        ! midpoint pressure ( Pa )
    real(r8),         intent(in)    :: pdel(:,:)         ! pressure thickness of levels (Pa)
    real(r8),         intent(in)    :: tfld(:,:)         ! temperature
    real(r8),         intent(in)    :: mbar(:,:)         ! mean wet atmospheric mass ( amu )
    real(r8), target, intent(in)    :: lwc(:,:)          ! cloud liquid water content (kg/kg)
    real(r8), target, intent(in)    :: cldfrc(:,:)       ! cloud fraction
    real(r8),         intent(in)    :: cldnum(:,:)       ! droplet number concentration (#/kg)
    real(r8),         intent(in)    :: invariants(:,:,:)
    real(r8), target, intent(inout) :: qcw(:,:,:)        ! cloud-borne aerosol (vmr)
    real(r8),         intent(inout) :: qin(:,:,:)        ! transported species ( vmr )
    real(r8),         intent(out)   :: xphlwc(:,:)       ! pH value multiplied by lwc

    real(r8),         intent(out)   :: aqso4(:,:)                   ! aqueous phase chemistry
    real(r8),         intent(out)   :: aqh2so4(:,:)                 ! aqueous phase chemistry
    real(r8),         intent(out)   :: aqso4_h2o2(:)                ! SO4 aqueous phase chemistry due to H2O2 (kg/m2)
    real(r8),         intent(out)   :: aqso4_o3(:)                  ! SO4 aqueous phase chemistry due to O3 (kg/m2)
    real(r8),         intent(in), optional :: yph_in                ! ph value
    real(r8),         intent(out), optional :: aqso4_h2o2_3d(:, :)  ! 3D SO4 aqueous phase chemistry due to H2O2 (kg/m2)
    real(r8),         intent(out), optional :: aqso4_o3_3d(:, :)    ! 3D SO4 aqueous phase chemistry due to O3 (kg/m2)

    !-----------------------------------------------------------------------
    !      ... Local variables
    !
    !           xhno3 ... in mixing ratio
    !-----------------------------------------------------------------------
    integer,  parameter :: itermax = 20
    real(r8), parameter :: ph0 = 5.0_r8  ! INITIAL PH VALUES
    real(r8), parameter :: PASCAL_TO_ATM = 1.0_r8 / 101325.0_r8 ! atm Pa-1
    real(r8), parameter :: M3_TO_L = 1.0e3_r8                   ! L m-3
    real(r8), parameter :: M3_TO_CM3 = 1.0e6_r8                 ! cm3 m-3
    real(r8), parameter :: G_TO_KG = 1.0e-3_r8                  ! kg g-1
    real(r8), parameter :: KMOL_TO_MOL = 1.0e3_r8               ! mol kmol-1
    real(r8), parameter :: SMALL_NUMBER = 1.0e-30_r8
    real(r8), parameter :: AVOGADRO = AVOGADRO_KMOL / KMOL_TO_MOL ! molecule mol-1
    real(r8), parameter :: const0 = 1.e3_r8/AVOGADRO
    real(r8)            :: MOLECULAR_WEIGHT_DRY_AIR ! kg mol-1
    real(r8), parameter :: MOLECULAR_WEIGHT_CO2 = MOLECULAR_WEIGHT_CO2_G_MOL * G_TO_KG ! kg mol-1
    real(r8), parameter :: xa0 = 11._r8
    real(r8), parameter :: xb0 = -.1_r8
    real(r8), parameter :: xa1 = 1.053_r8
    real(r8), parameter :: xb1 = -4.368_r8
    real(r8), parameter :: xa2 = 1.016_r8
    real(r8), parameter :: xb2 = -2.54_r8
    real(r8), parameter :: xa3 = .816e-32_r8
    real(r8), parameter :: xb3 = .259_r8

    real(r8), parameter :: kh0 = 690.0_r8          ! HO2(g)          -> Ho2(a)        Reference: JPL 19-5
    real(r8), parameter :: kh1 = 1.6e-5_r8         ! HO2(a)          -> H+ + O2-      Reference: JPL 19-5
    real(r8), parameter :: kh2 = 8.3e5_r8          ! HO2(a) + ho2(a) -> h2o2(a) + o2  Reference: JPL; Bielski et al. 1985
    real(r8), parameter :: kh3 = 9.7e7_r8          ! HO2(a) + o2-    -> h2o2(a) + o2  Reference: JPL; Bielski et al. 1985
    real(r8), parameter :: Ra = GAS_CONSTANT_KMOL / KMOL_TO_MOL * M3_TO_L * PASCAL_TO_ATM ! universal constant   (atm)/(M-K)

    !
    real(r8) :: xdelso4hp(ncol,pver)
    real(r8) :: xhnm(ncol,pver) ! air number density (molecules cm-3)

    integer  :: k, i, iter, file
    real(r8) :: wrk, delta
    real(r8) :: xph0, aden, xk, xe, x2
    real(r8) :: tz, xl, px, qz, es, qs, patm
    real(r8) :: Eso2, Eso4, Ehno3, Eco2, Eh2o, Enh3
    real(r8) :: so2g, h2o2g, o3g
    real(r8) :: hno3a, nh3a, so2a, h2o2a, co2a, o3a
    real(r8) :: rah2o2, rao3, pso4, ccc
    real(r8) :: cnh3, chno3, com, com1, com2, xra
    real(r8) :: f_hso3 ! fraction of aqueous S(IV) that's HSO3-
    real(r8) :: f_so3  ! fraction of aqueous S(IV) that's SO3=

    real(r8) :: hno3g(ncol,pver), nh3g(ncol,pver)
    !
    !-----------------------------------------------------------------------
    !            for Ho2(g) -> H2o2(a) formation
    !            schwartz JGR, 1984, 11589
    !-----------------------------------------------------------------------
    real(r8) :: kh4    ! kh2+kh3
    real(r8) :: xam    ! air density /cm3
    real(r8) :: ho2s   ! ho2s = ho2(a)+o2-
    real(r8) :: r1h2o2 ! prod(h2o2) by ho2 in mole/L(w)/s
    real(r8) :: r2h2o2 ! prod(h2o2) by ho2 in mix/s

    real(r8), dimension(ncol,pver)  ::             &
         xhno3, xh2o2, xso2, xso4, xno3, xco2, &
         xnh3, xnh4, xo3,         &
         cfact, &
         xph, xho2,         &
         xh2so4, xmsa, xso4_init, &
         hehno3, &            ! henry law const for hno3
         heh2o2, &            ! henry law const for h2o2
         heso2,  &            ! henry law const for so2
         henh3,  &            ! henry law const for nh3
         heo3              !!,   &            ! henry law const for o3

    real(r8), pointer :: co2_mass_mixing_ratio(:,:) ! kg kg-1

    real(r8), dimension(ncol)  :: work1
    logical :: converged

    real(r8), pointer :: xso4c(:,:)
    real(r8), pointer :: xnh4c(:,:)
    real(r8), pointer :: xno3c(:,:)
    type(cldaero_conc_t), pointer :: cldconc

    real(r8) :: fact1_hno3, fact2_hno3, fact3_hno3
    real(r8) :: fact1_so2, fact2_so2, fact3_so2, fact4_so2
    real(r8) :: fact1_nh3, fact2_nh3, fact3_nh3
    real(r8) :: fact1_co2, fact2_co2, fact3_co2, fact4_co2
    real(r8) :: tmp_hp, tmp_hso3, tmp_hco3, tmp_nh4, tmp_no3
    real(r8) :: tmp_oh, tmp_so3, tmp_so4
    real(r8) :: tmp_neg, tmp_pos
    real(r8) :: yph, yph_lo, yph_hi
    real(r8) :: ynetpos, ynetpos_lo, ynetpos_hi

    MOLECULAR_WEIGHT_DRY_AIR = MOLECULAR_WEIGHT_DRY_AIR_G_MOL * G_TO_KG  ! kg mol-1

    !-----------------------------------------------------------------
    !       ... NOTE: The press array is in pascals and must be
    !                 mutiplied by 10 to yield dynes/cm**2.
    !-----------------------------------------------------------------
    !==================================================================
    !       ... First set the PH
    !==================================================================
    !      ... Initial values
    !           The values of so2, so4 are after (1) SLT, and CHEM
    !-----------------------------------------------------------------
    xhnm(:ncol,:) = press(:ncol,:) / (tfld(:ncol,:) * M3_TO_CM3 * BOLTZMANN)  ! air number density (molecules cm-3)

    call rad_cnst_get_gas(0, 'CO2', state, pbuf, co2_mass_mixing_ratio)

    xph0 = 10._r8**(-ph0)                      ! initial PH value

    do k = 1,pver
       cfact(:,k) = xhnm(:,k)     &          ! /cm3(a)
            * 1.e6_r8             &          ! /m3(a)
            * BOLTZMANN/287._r8   &          ! Kg(a)/m3(a)
            * 1.e-3_r8                       ! Kg(a)/L(a)
    end do

    cldconc => sox_cldaero_create_obj( cldfrc,qcw,lwc, cfact, ncol, loffset )
    xso4c => cldconc%so4c
    xnh4c => cldconc%nh4c
    xno3c => cldconc%no3c

    xso4(:,:) = 0._r8
    xno3(:,:) = 0._r8
    xnh4(:,:) = 0._r8

    do k = 1,pver
       xph(:,k) = xph0                                ! initial PH value

       xco2(:ncol,k) = co2_mass_mixing_ratio(:ncol,k) &
                   * (MOLECULAR_WEIGHT_DRY_AIR / MOLECULAR_WEIGHT_CO2)  ! mixing ratio

       if ( inv_so2 ) then
          xso2 (:,k) = invariants(:,k,id_so2)/xhnm(:,k)  ! mixing ratio
       else
          xso2 (:,k) = qin(:,k,id_so2)                   ! mixing ratio
       endif

       if (id_hno3 > 0) then
          xhno3(:,k) = qin(:,k,id_hno3)
       else
          xhno3(:,k) = 0.0_r8
       endif

       if ( inv_h2o2 ) then
          xh2o2 (:,k) = invariants(:,k,id_h2o2)/xhnm(:,k)  ! mixing ratio
       else
          xh2o2 (:,k) = qin(:,k,id_h2o2)                   ! mixing ratio
       endif

       if (id_nh3  > 0) then
          xnh3 (:,k) = qin(:,k,id_nh3)
       else
          xnh3 (:,k) = 0.0_r8
       endif

       if ( inv_o3 ) then
          xo3  (:,k) = invariants(:,k,id_o3)/xhnm(:,k) ! mixing ratio
       else
          xo3  (:,k) = qin(:,k,id_o3)                  ! mixing ratio
       endif
       if ( inv_ho2 ) then
          xho2 (:,k) = invariants(:,k,id_ho2)/xhnm(:,k)! mixing ratio
       else
          xho2 (:,k) = qin(:,k,id_ho2)                 ! mixing ratio
       endif

       if (cloud_borne) then
          xh2so4(:,k) = qin(:,k,id_h2so4)
       else
          xso4  (:,k) = qin(:,k,id_so4) ! mixing ratio
       endif
       if (id_msa > 0) xmsa (:,k) = qin(:,k,id_msa)

    end do

    !-----------------------------------------------------------------
    !       ... Temperature dependent Henry constants
    !-----------------------------------------------------------------
    ver_loop0: do k = 1,pver                               !! pver loop for STEP 0
       col_loop0: do i = 1,ncol

          if (cloud_borne .and. cldfrc(i,k)>0._r8) then
             xso4(i,k) = xso4c(i,k) / cldfrc(i,k)
             xnh4(i,k) = xnh4c(i,k) / cldfrc(i,k)
             xno3(i,k) = xno3c(i,k) / cldfrc(i,k)
          endif
          xl = cldconc%xlwc(i,k)

          if( xl >= 1.e-8_r8 ) then
             work1(i) = 1._r8 / tfld(i,k) - 1._r8 / 298._r8

             !-----------------------------------------------------------------
             ! 21-mar-2011 changes by rce
             ! ph calculation now uses bisection method to solve the electro-neutrality equation
             ! 3-mode aerosols (where so4 is assumed to be nh4hso4)
             !    old code set xnh4c = so4c
             !    new code sets xnh4c = 0, then uses a -1 charge (instead of -2)
             !       for so4 when solving the electro-neutrality equation
             !-----------------------------------------------------------------

             !-----------------------------------------------------------------
             !  calculations done before iterating
             !-----------------------------------------------------------------

             !-----------------------------------------------------------------
             tz = tfld(i,k)
             patm = press(i,k) * PASCAL_TO_ATM ! atm
             xam  = press(i,k)/(BOLTZMANN*tz)  ! air density /M3

             !-----------------------------------------------------------------
             !        ... hno3
             !-----------------------------------------------------------------
             ! previous code
             !    hehno3(i,k)  = xk*(1._r8 + xe/xph(i,k))
             !    px = hehno3(i,k) * Ra * tz * xl
             !    hno3g = xhno3(i,k)/(1._r8 + px)
             !    Ehno3 = xk*xe*hno3g *patm
             ! equivalent new code
             !    hehno3 = xk + xk*xe/hplus
             !    hno3g = xhno3/(1 + px)
             !          = xhno3/(1 + hehno3*ra*tz*xl)
             !          = xhno3/(1 + xk*ra*tz*xl*(1 + xe/hplus)
             !    ehno3 = hno3g*xk*xe*patm
             !          = xk*xe*patm*xhno3/(1 + xk*ra*tz*xl*(1 + xe/hplus)
             !          = ( fact1_hno3    )/(1 + fact2_hno3 *(1 + fact3_hno3/hplus)
             !    [hno3-] = ehno3/hplus
             xk = dheff(1,heff_id_hno3) * exp( dheff(2,heff_id_hno3) * work1(i) )
             xe = dheff(3,heff_id_hno3) * exp( dheff(4,heff_id_hno3) * work1(i) )
             fact1_hno3 = xk*xe*patm*xhno3(i,k)
             fact2_hno3 = xk*ra*tz*xl
             fact3_hno3 = xe

             !-----------------------------------------------------------------
             !          ... so2
             !-----------------------------------------------------------------
             ! previous code
             !    heso2(i,k)  = xk*(1._r8 + wrk*(1._r8 + x2/xph(i,k)))
             !    px = heso2(i,k) * Ra * tz * xl
             !    so2g =  xso2(i,k)/(1._r8+ px)
             !    Eso2 = xk*xe*so2g *patm
             ! equivalent new code
             !    heso2 = xk + xk*xe/hplus * xk*xe*x2/hplus**2
             !    so2g = xso2/(1 + px)
             !         = xso2/(1 + heso2*ra*tz*xl)
             !         = xso2/(1 + xk*ra*tz*xl*(1 + (xe/hplus)*(1 + x2/hplus))
             !    eso2 = so2g*xk*xe*patm
             !          = xk*xe*patm*xso2/(1 + xk*ra*tz*xl*(1 + (xe/hplus)*(1 + x2/hplus))
             !          = ( fact1_so2    )/(1 + fact2_so2 *(1 + (fact3_so2/hplus)*(1 + fact4_so2/hplus)
             !    [hso3-] + 2*[so3--] = (eso2/hplus)*(1 + 2*x2/hplus)
             xk = dheff(1,heff_id_so2) * exp( dheff(2,heff_id_so2) * work1(i) )
             xe = dheff(3,heff_id_so2) * exp( dheff(4,heff_id_so2) * work1(i) )
             x2 = dheff(5,heff_id_so2) * exp( dheff(6,heff_id_so2) * work1(i) )
             fact1_so2 = xk*xe*patm*xso2(i,k)
             fact2_so2 = xk*ra*tz*xl
             fact3_so2 = xe
             fact4_so2 = x2

             !-----------------------------------------------------------------
             !          ... nh3
             !-----------------------------------------------------------------
             ! previous code
             !    henh3(i,k)  = xk*(1._r8 + xe*xph(i,k)/xkw)
             !    px = henh3(i,k) * Ra * tz * xl
             !    nh3g = (xnh3(i,k)+xnh4(i,k))/(1._r8+ px)
             !    Enh3 = xk*xe*nh3g/xkw *patm
             ! equivalent new code
             !    henh3 = xk + xk*xe*hplus/xkw
             !    nh3g = xnh34/(1 + px)
             !         = xnh34/(1 + henh3*ra*tz*xl)
             !         = xnh34/(1 + xk*ra*tz*xl*(1 + xe*hplus/xkw)
             !    enh3 = nh3g*xk*xe*patm/xkw
             !          = ((xk*xe*patm/xkw)*xnh34)/(1 + xk*ra*tz*xl*(1 + xe*hplus/xkw)
             !          = ( fact1_nh3            )/(1 + fact2_nh3  *(1 + fact3_nh3*hplus)
             !    [nh4+] = enh3*hplus
             ! NOTE: Algorithm modified to follow that used in wet deposition.
             !       This essentially replaces xkw (1.0e-14) with a temperature
             !       dependent value for the water dissociation constant, x2.
             xk = dheff(1,heff_id_nh3) * exp( dheff(2,heff_id_nh3) * work1(i) )
             xe = dheff(3,heff_id_nh3) * exp( dheff(4,heff_id_nh3) * work1(i) )
             x2 = dheff(5,heff_id_nh3) * exp( dheff(6,heff_id_nh3) * work1(i) )
             fact1_nh3 = (xk*xe*patm/x2)*(xnh3(i,k)+xnh4(i,k))
             fact2_nh3 = xk*ra*tz*xl
             fact3_nh3 = xe/x2

             !-----------------------------------------------------------------
             !        ... h2o effects
             ! NOTE: Algorithm modified to follow that used in wet deposition.
             !       This essentially replaces xkw (1.0e-14) with a temperature
             !       dependent value for the water dissociation constant, x2
             !       (calculated above with NH4 Heff terms).
             !-----------------------------------------------------------------
             Eh2o = x2

             !-----------------------------------------------------------------
             !        ... co2 effects
             ! NOTE: Algorithm modified to follow that used in wet deposition.
             !       This now applies the same algorithm for diprotic acids used
             !       for SO2.
             !-----------------------------------------------------------------
             xk = dheff(1,heff_id_co2) * exp( dheff(2,heff_id_co2) * work1(i) )
             xe = dheff(3,heff_id_co2) * exp( dheff(4,heff_id_co2) * work1(i) )
             x2 = dheff(5,heff_id_co2) * exp( dheff(6,heff_id_co2) * work1(i) )
             fact1_co2 = xk*xe*patm*xco2(i,k)
             fact2_co2 = xk*ra*tz*xl
             fact3_co2 = xe
             fact4_co2 = x2

             !-----------------------------------------------------------------
             !         ... so4 effect
             !-----------------------------------------------------------------
             Eso4 = xso4(i,k)*xhnm(i,k)   &         ! /cm3(a)
                  *const0/xl


             !-----------------------------------------------------------------
             ! now use bisection method to solve electro-neutrality equation
             !
             ! during the iteration loop,
             !    yph_lo = lower ph value that brackets the root (i.e., correct ph)
             !    yph_hi = upper ph value that brackets the root (i.e., correct ph)
             !    yph    = current ph value
             !    yposnet_lo and yposnet_hi = net positive ions for
             !       yph_lo and yph_hi
             !-----------------------------------------------------------------
             do iter = 1,itermax

                if (.not. present(yph_in)) then
                   if (iter == 1) then
                      ! 1st iteration ph = lower bound value
                      yph_lo = 2.0_r8
                      yph_hi = yph_lo
                      yph = yph_lo
                   else if (iter == 2) then
                      ! 2nd iteration ph = upper bound value
                      yph_hi = 7.0_r8
                      yph = yph_hi
                   else
                      ! later iteration ph = mean of the two bracketing values
                      yph = 0.5_r8*(yph_lo + yph_hi)
                   end if
                else
                   yph = yph_in
                end if

                ! calc current [H+] from ph
                xph(i,k) = 10.0_r8**(-yph)


                !-----------------------------------------------------------------
                !        ... hno3
                !-----------------------------------------------------------------
                Ehno3 = fact1_hno3/(1.0_r8 + fact2_hno3*(1.0_r8 + fact3_hno3/xph(i,k)))

                !-----------------------------------------------------------------
                !          ... so2
                !-----------------------------------------------------------------
                Eso2 = fact1_so2/(1.0_r8 + fact2_so2*(1.0_r8 + (fact3_so2/xph(i,k)) &
                     *(1.0_r8 +  fact4_so2/xph(i,k))))

                !-----------------------------------------------------------------
                !          ... nh3
                !-----------------------------------------------------------------
                Enh3 = fact1_nh3/(1.0_r8 + fact2_nh3*(1.0_r8 + fact3_nh3*xph(i,k)))

                !-----------------------------------------------------------------
                !          ... co2
                !-----------------------------------------------------------------
                Eco2 = fact1_co2/(1.0_r8 + fact2_co2*(1.0_r8 + (fact3_co2/xph(i,k)) &
                     *(1.0_r8 +  fact4_co2/xph(i,k))))

                tmp_nh4  = Enh3 * xph(i,k)
                tmp_hso3 = Eso2 / xph(i,k)
                tmp_so3  = tmp_hso3 * 2.0_r8*fact4_so2/xph(i,k)
                tmp_hco3 = Eco2 / xph(i,k)
                tmp_oh   = Eh2o / xph(i,k)
                tmp_no3  = Ehno3 / xph(i,k)
                tmp_so4 = cldconc%so4_fact*Eso4
                tmp_pos = xph(i,k) + tmp_nh4
                tmp_neg = tmp_oh + tmp_hco3 + tmp_no3 + tmp_hso3 + tmp_so3 + tmp_so4

                ynetpos = tmp_pos - tmp_neg


                ! yposnet = net positive ions/charge
                ! if the correct ph is bracketed by yph_lo and yph_hi (with yph_lo < yph_hi),
                !    then you will have yposnet_lo > 0 and yposnet_hi < 0
                converged = .false.
                if (iter > 2) then
                   if (ynetpos == 0.0_r8) then
                      ! the exact solution was found (very unlikely)
                      tmp_hp = xph(i,k)
                      converged = .true.
                      exit
                   else if (ynetpos >= 0.0_r8) then
                      ! net positive ions are >= 0 for both yph and yph_lo
                      !    so replace yph_lo with yph
                      yph_lo = yph
                      ynetpos_lo = ynetpos
                   else
                      ! net positive ions are <= 0 for both yph and yph_hi
                      !    so replace yph_hi with yph
                      yph_hi = yph
                      ynetpos_hi = ynetpos
                   end if

                   if (abs(yph_hi - yph_lo) .le. 0.005_r8) then
                      ! |yph_hi - yph_lo| <= convergence criterion, so set
                      !    final ph to their midpoint and exit
                      ! (.005 absolute error in pH gives .01 relative error in H+)
                      tmp_hp = xph(i,k)
                      yph = 0.5_r8*(yph_hi + yph_lo)
                      xph(i,k) = 10.0_r8**(-yph)
                      converged = .true.
                      exit
                   else
                      ! do another iteration
                      converged = .false.
                   end if

                else if (iter == 1) then
                   if (ynetpos <= 0.0_r8) then
                      ! the lower and upper bound ph values (2.0 and 7.0) do not bracket
                      !    the correct ph, so use the lower bound
                      tmp_hp = xph(i,k)
                      converged = .true.
                      exit
                   end if
                   ynetpos_lo = ynetpos

                else ! (iter == 2)
                   if (ynetpos >= 0.0_r8) then
                      ! the lower and upper bound ph values (2.0 and 7.0) do not bracket
                      !    the correct ph, so use they upper bound
                      tmp_hp = xph(i,k)
                      converged = .true.
                      exit
                   end if
                   ynetpos_hi = ynetpos
                end if

             end do ! iter

             if( .not. converged ) then
                write(iulog,*) 'setsox: pH failed to converge @ (',i,',',k,'), % change=', &
                     100._r8*delta
             end if
          else
             xph(i,k) =  1.e-7_r8
          end if
       end do col_loop0
    end do ver_loop0 ! end pver loop for STEP 0

    !==============================================================
    !          ... Now use the actual PH
    !==============================================================
    ver_loop1: do k = 1,pver
       col_loop1: do i = 1,ncol
          work1(i) = 1._r8 / tfld(i,k) - 1._r8 / 298._r8
          tz = tfld(i,k)

          xl = cldconc%xlwc(i,k)

          patm = press(i,k) * PASCAL_TO_ATM   ! atm
          xam  = press(i,k)/(BOLTZMANN*tz)    ! air density /M3

          !-----------------------------------------------------------------------
          !        ... hno3
          !-----------------------------------------------------------------------
          xk = dheff(1,heff_id_hno3) * exp( dheff(2,heff_id_hno3) * work1(i) )
          xe = dheff(3,heff_id_hno3) * exp( dheff(4,heff_id_hno3) * work1(i) )
          hehno3(i,k)  = xk*(1._r8 + xe/xph(i,k))

          !-----------------------------------------------------------------
          !        ... h2o2
          !-----------------------------------------------------------------
          xk = dheff(1,heff_id_h2o2) * exp( dheff(2,heff_id_h2o2) * work1(i) )
          xe = dheff(3,heff_id_h2o2) * exp( dheff(4,heff_id_h2o2) * work1(i) )
          heh2o2(i,k)  = xk*(1._r8 + xe/xph(i,k))

          !-----------------------------------------------------------------
          !         ... so2
          !-----------------------------------------------------------------
          xk = dheff(1,heff_id_so2) * exp( dheff(2,heff_id_so2) * work1(i) )
          xe = dheff(3,heff_id_so2) * exp( dheff(4,heff_id_so2) * work1(i) )
          x2 = dheff(5,heff_id_so2) * exp( dheff(6,heff_id_so2) * work1(i) )

          wrk = xe/xph(i,k)
          heso2(i,k)  = xk*(1._r8 + wrk*(1._r8 + x2/xph(i,k)))

          ! Calculate fraction of total aqueous S(IV) that is HSO3- and SO3=
          wrk = xe/(xph(i,k)*xph(i,k) + xe*xph(i,k) + xe*x2)
          f_hso3 = wrk*xph(i,k) ! [HSO3-]/[S(IV)]
          f_so3  = wrk*x2       !  [SO3=]/[S(IV)]

          !-----------------------------------------------------------------
          !          ... nh3
          !-----------------------------------------------------------------
          xk = dheff(1,heff_id_nh3) * exp( dheff(2,heff_id_nh3) * work1(i) )
          xe = dheff(3,heff_id_nh3) * exp( dheff(4,heff_id_nh3) * work1(i) )
          x2 = dheff(5,heff_id_nh3) * exp( dheff(6,heff_id_nh3) * work1(i) )
          henh3(i,k)  = xk*(1._r8 + xe*xph(i,k)/x2)

          !-----------------------------------------------------------------
          !        ... o3
          !-----------------------------------------------------------------
          xk = dheff(1,heff_id_o3) * exp( dheff(2,heff_id_o3) * work1(i) )
          heo3(i,k) = xk

          !------------------------------------------------------------------------
          !       ... for Ho2(g) -> H2o2(a) formation
          !           schwartz JGR, 1984, 11589
          !------------------------------------------------------------------------
          kh4 = (kh2 + kh3*kh1/xph(i,k)) / ((1._r8 + kh1/xph(i,k))**2)
          ho2s = kh0*xho2(i,k)*patm*(1._r8 + kh1/xph(i,k))  ! ho2s = ho2(a)+o2-
          r1h2o2 = kh4*ho2s*ho2s                            ! prod(h2o2) in mole/L(w)/s
          r2h2o2 = r1h2o2 * xl  & ! (mole(h2o2)/L(w)/s) * (L(w)/L(a))
                  / xam         & ! / (molecule(a)/m3(a))
                  * AVOGADRO    & ! * (molecule(a)/mole(a))
                  * M3_TO_L       ! * (L(a)/m3(a)) = mole(h2o2)/mole(a)/s
          
          xh2o2(i,k) = xh2o2(i,k) + r2h2o2*dtime ! updated h2o2 by het production
          
          !-----------------------------------------------
          !       ... Partioning
          !-----------------------------------------------

          !-----------------------------------------------------------------
          !        ... hno3
          !-----------------------------------------------------------------
          px = hehno3(i,k) * Ra * tz * xl
          hno3g(i,k) = (xhno3(i,k)+xno3(i,k))/(1._r8 + px)

          !------------------------------------------------------------------------
          !        ... h2o2
          !------------------------------------------------------------------------
          px = heh2o2(i,k) * Ra * tz * xl
          h2o2g =  xh2o2(i,k)/(1._r8+ px)

          !------------------------------------------------------------------------
          !         ... so2
          !------------------------------------------------------------------------
          px = heso2(i,k) * Ra * tz * xl
          so2g =  xso2(i,k)/(1._r8+ px)

          !------------------------------------------------------------------------
          !         ... o3
          !------------------------------------------------------------------------
          px = heo3(i,k) * Ra * tz * xl
          o3g =  xo3(i,k)/(1._r8+ px)

          !------------------------------------------------------------------------
          !         ... nh3
          !------------------------------------------------------------------------
          px = henh3(i,k) * Ra * tz * xl
          if (id_nh3>0) then
             nh3g(i,k) = (xnh3(i,k)+xnh4(i,k))/(1._r8+ px)
          else
             nh3g(i,k) = 0._r8
          endif

          !-----------------------------------------------
          !       ... Aqueous phase reaction rates
          !           SO2 + H2O2 -> SO4
          !           SO2 + O3   -> SO4
          !-----------------------------------------------

          !------------------------------------------------------------------------
          !       ... S(IV) (HSO3) + H2O2
          ! Reference: Seinfeld and Pandis textbook (chapter 6);
          !            original source: Hoffmann and Calvert (1985)
          !------------------------------------------------------------------------
          rah2o2 = 7.45e7_r8 * EXP( -4430.0_r8*work1(i) ) * xph(i,k) &
                   / (1.0_r8 + 13.0_r8*xph(i,k)) * f_hso3

          !------------------------------------------------------------------------
          !        ... S(IV)+ O3
          ! Reference: Seinfeld and Pandis textbook (chapter 6);
          !            original source: Hoffmann and Calvert (1985)
          !------------------------------------------------------------------------
          rao3   =   3.75e5_r8 * EXP(-5530.0_r8*work1(i)) * f_hso3 &
                   + 1.59e9_r8 * EXP(-5280.0_r8*work1(i)) * f_so3

          !-----------------------------------------------------------------
          !       ... Prediction after aqueous phase
          !       so4
          !       When Cloud is present
          !
          !       S(IV) + H2O2 = S(VI)
          !       S(IV) + O3   = S(VI)
          !
          !       reference:
          !           (1) Seinfeld
          !           (2) Benkovitz
          !-----------------------------------------------------------------

          !............................
          !       S(IV) + H2O2 = S(VI)
          !............................

          IF (XL .ge. 1.e-8_r8) THEN    !! WHEN CLOUD IS PRESENTED

             pso4 = rah2o2 * heh2o2(i,k) * h2o2g * patm  &
                    * heso2(i,k)  * so2g  * patm  ! [M/s]

             pso4 = pso4 & ! [M/s] = [mole/L(w)/s]
                  * xl & ! [mole/L(a)/s]
                  / const0 & ! [/L(a)/s]
                  / xhnm(i,k)

             ! estimate the net production of so4, without exceeding reactant concentrations
             xso4_init(i,k) = xso4(i,k)
             ccc            = max(min(pso4*dtime, min(xh2o2(i,k) - SMALL_NUMBER, xso2(i,k) - SMALL_NUMBER)), 0.0_r8)
             xso4(i,k)      = xso4(i,k) + ccc
             xh2o2(i,k)     = xh2o2(i,k) - ccc
             xso2(i,k)      = xso2(i,k)  - ccc
             xdelso4hp(i,k) = ccc

             !...........................
             !       S(IV) + O3 = S(VI)
             !...........................

             pso4 = rao3 * heo3(i,k)*o3g*patm * heso2(i,k)*so2g*patm  ! [M/s]

             pso4 = pso4        &                                ! [M/s] =  [mole/L(w)/s]
                  * xl          &                                ! [mole/L(a)/s]
                  / const0      &                                ! [/L(a)/s]
                  / xhnm(i,k)                                    ! [mixing ratio/s]

             ! estimate the net production of so4, without exceeding reactant concentrations
             xso4_init(i,k) = xso4(i,k)
             ccc            = max(min(pso4*dtime, xso2(i,k) - SMALL_NUMBER), 0.0_r8)
             xso4(i,k)      = xso4(i,k) + ccc
             xso2(i,k)      = xso2(i,k) - ccc

          END IF !! WHEN CLOUD IS PRESENTED

       end do col_loop1
    end do ver_loop1

    call sox_cldaero_update( &
          state, ncol, lchnk, loffset, dtime, mbar, pdel, press, tfld, cldnum, cldfrc, cfact, cldconc%xlwc, &
          xdelso4hp, xh2so4, xso4, xso4_init, nh3g, hno3g, xnh3, xhno3, xnh4c,  xno3c, xmsa, xso2, xh2o2, qcw, qin, &
          aqso4, aqh2so4, aqso4_h2o2, aqso4_o3, aqso4_h2o2_3d=aqso4_h2o2_3d, aqso4_o3_3d=aqso4_o3_3d )

    xphlwc(:,:) = 0._r8
    do k = 1, pver
       do i = 1, ncol
          if (cldfrc(i,k)>=1.e-5_r8 .and. lwc(i,k)>=1.e-8_r8) then
             xphlwc(i,k) = -1._r8*log10(xph(i,k)) * lwc(i,k)
          endif
       end do
    end do

    call sox_cldaero_destroy_obj(cldconc)

  end subroutine setsox

   !-----------------------------------------------------------------
   !       ... looks up Effective Henry's Law Constant parameters
   !-----------------------------------------------------------------
   pure integer function get_heff_index(species_name) result(index)
      use shr_drydep_mod, only: species_name_table
     
      character(len=*), intent(in) :: species_name

      do index = 1, size(species_name_table)
         if (trim(adjustl(species_name)) == &
             trim(adjustl(species_name_table(index)))) return
      end do
      index = -1
   end function get_heff_index

end module mo_setsox
