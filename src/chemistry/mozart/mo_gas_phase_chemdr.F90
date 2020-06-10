module mo_gas_phase_chemdr

  use shr_kind_mod,     only : r8 => shr_kind_r8
  use shr_const_mod,    only : pi => shr_const_pi
  use constituents,     only : pcnst
  use cam_history,      only : fieldname_len
  use chem_mods,        only : phtcnt, rxntot, gas_pcnst
  use chem_mods,        only : rxt_tag_cnt, rxt_tag_lst, rxt_tag_map, extcnt, num_rnts
  use dust_model,       only : dust_names, ndust => dust_nbin
  use ppgrid,           only : pcols, pver
  use phys_control,     only : phys_getopts
  use carma_flags_mod,  only : carma_hetchem_feedback
  use chem_prod_loss_diags, only: chem_prod_loss_diags_init, chem_prod_loss_diags_out

  implicit none
  save

  private
  public :: gas_phase_chemdr, gas_phase_chemdr_inti 
  public :: map2chm

  integer :: map2chm(pcnst) = 0           ! index map to/from chemistry/constituents list

  integer :: synoz_ndx, so4_ndx, h2o_ndx, o2_ndx, o_ndx, hno3_ndx, hcl_ndx, dst_ndx, cldice_ndx, snow_ndx
  integer :: o3_ndx, o3s_ndx
  integer :: het1_ndx
  integer :: ndx_cldfr, ndx_cmfdqr, ndx_nevapr, ndx_cldtop, ndx_prain
  integer :: ndx_h2so4
!
! CCMI
!
  integer :: st80_25_ndx
  integer :: st80_25_tau_ndx
  integer :: aoa_nh_ndx
  integer :: aoa_nh_ext_ndx
  integer :: nh_5_ndx
  integer :: nh_50_ndx
  integer :: nh_50w_ndx
  integer :: sad_pbf_ndx
  integer :: cb1_ndx,cb2_ndx,oc1_ndx,oc2_ndx,dst1_ndx,dst2_ndx,sslt1_ndx,sslt2_ndx
  integer :: soa_ndx,soai_ndx,soam_ndx,soat_ndx,soab_ndx,soax_ndx

  character(len=fieldname_len),dimension(rxt_tag_cnt)   :: tag_names
  character(len=fieldname_len),dimension(extcnt)        :: extfrc_name

  logical :: pm25_srf_diag
  logical :: pm25_srf_diag_soa

  logical :: convproc_do_aer
  integer :: ele_temp_ndx, ion_temp_ndx

contains

  subroutine gas_phase_chemdr_inti()

    use mo_chem_utls,      only : get_spc_ndx, get_extfrc_ndx, get_rxt_ndx
    use cam_history,       only : addfld,add_default,horiz_only
    use mo_chm_diags,      only : chm_diags_inti
    use constituents,      only : cnst_get_ind
    use physics_buffer,    only : pbuf_get_index
    use rate_diags,        only : rate_diags_init
    use cam_abortutils,    only : endrun

    implicit none

    character(len=3) :: string
    integer          :: n, m, err, ii
    logical :: history_cesm_forcing
    character(len=16) :: unitstr
    !-----------------------------------------------------------------------
    logical :: history_scwaccm_forcing

    call phys_getopts( history_scwaccm_forcing_out = history_scwaccm_forcing )

    call phys_getopts( convproc_do_aer_out = convproc_do_aer, history_cesm_forcing_out=history_cesm_forcing )
   
    ndx_h2so4 = get_spc_ndx('H2SO4')
!
! CCMI
!
    st80_25_ndx     = get_spc_ndx   ('ST80_25')
    st80_25_tau_ndx = get_rxt_ndx   ('ST80_25_tau')
    aoa_nh_ndx      = get_spc_ndx   ('AOA_NH')
    aoa_nh_ext_ndx  = get_extfrc_ndx('AOA_NH')
    nh_5_ndx        = get_spc_ndx('NH_5')
    nh_50_ndx       = get_spc_ndx('NH_50')
    nh_50w_ndx      = get_spc_ndx('NH_50W')
!
    cb1_ndx         = get_spc_ndx('CB1')
    cb2_ndx         = get_spc_ndx('CB2')
    oc1_ndx         = get_spc_ndx('OC1')
    oc2_ndx         = get_spc_ndx('OC2')
    dst1_ndx        = get_spc_ndx('DST01')
    dst2_ndx        = get_spc_ndx('DST02')
    sslt1_ndx       = get_spc_ndx('SSLT01')
    sslt2_ndx       = get_spc_ndx('SSLT02')
    soa_ndx         = get_spc_ndx('SOA')
    soam_ndx         = get_spc_ndx('SOAM')
    soai_ndx         = get_spc_ndx('SOAI')
    soat_ndx         = get_spc_ndx('SOAT')
    soab_ndx         = get_spc_ndx('SOAB')
    soax_ndx         = get_spc_ndx('SOAX')

    pm25_srf_diag = cb1_ndx>0 .and. cb2_ndx>0 .and. oc1_ndx>0 .and. oc2_ndx>0 &
              .and. dst1_ndx>0 .and. dst2_ndx>0 .and. sslt1_ndx>0 .and. sslt2_ndx>0 &
              .and. soa_ndx>0 

    pm25_srf_diag_soa = cb1_ndx>0 .and. cb2_ndx>0 .and. oc1_ndx>0 .and. oc2_ndx>0 &
              .and. dst1_ndx>0 .and. dst2_ndx>0 .and. sslt1_ndx>0 .and. sslt2_ndx>0 &
              .and. soam_ndx>0 .and. soai_ndx>0 .and. soat_ndx>0 .and. soab_ndx>0 .and. soax_ndx>0
    
    if ( pm25_srf_diag .or. pm25_srf_diag_soa) then
       call addfld('PM25_SRF',horiz_only,'I','kg/kg','bottom layer PM2.5 mixing ratio' )
    endif
    call addfld('U_SRF',horiz_only,'I','m/s','bottom layer wind velocity' )
    call addfld('V_SRF',horiz_only,'I','m/s','bottom layer wind velocity' )
    call addfld('Q_SRF',horiz_only,'I','kg/kg','bottom layer specific humidity' )
!
    het1_ndx= get_rxt_ndx('het1')
    o3_ndx  = get_spc_ndx('O3')
    o3s_ndx = get_spc_ndx('O3S')
    o_ndx   = get_spc_ndx('O')
    o2_ndx  = get_spc_ndx('O2')
    so4_ndx = get_spc_ndx('SO4')
    h2o_ndx = get_spc_ndx('H2O')
    hno3_ndx = get_spc_ndx('HNO3')
    hcl_ndx  = get_spc_ndx('HCL')
    dst_ndx = get_spc_ndx( dust_names(1) )
    synoz_ndx = get_extfrc_ndx( 'SYNOZ' )
    call cnst_get_ind( 'CLDICE', cldice_ndx )
    call cnst_get_ind( 'SNOWQM', snow_ndx, abort=.false. )


    do m = 1,extcnt
       WRITE(UNIT=string, FMT='(I2.2)') m
       extfrc_name(m) = 'extfrc_'// trim(string)
       call addfld( extfrc_name(m), (/ 'lev' /), 'I', ' ', 'ext frcing' )
    end do

    do n = 1,rxt_tag_cnt
       tag_names(n) = trim(rxt_tag_lst(n))
       if (n<=phtcnt) then
          call addfld( tag_names(n), (/ 'lev' /), 'I', '/s', 'photolysis rate constant' )
       else
          ii = n-phtcnt
          select case(num_rnts(ii))
          case(1)
             unitstr='/s'
          case(2)
             unitstr='cm3/molecules/s'
          case(3)
             unitstr='cm6/molecules2/s'
          case default
             call endrun('gas_phase_chemdr_inti: invalid value in num_rnts used to set units in reaction rate constant')
          end select
          call addfld( tag_names(n), (/ 'lev' /), 'I', unitstr, 'reaction rate constant' )
       endif
       if (history_scwaccm_forcing) then
          select case (trim(tag_names(n)))
          case ('jh2o_a', 'jh2o_b', 'jh2o_c' )
             call add_default( tag_names(n), 1, ' ')
          end select
       endif
    enddo

    call addfld( 'DTCBS',   horiz_only, 'I', ' ','photolysis diagnostic black carbon OD' )
    call addfld( 'DTOCS',   horiz_only, 'I', ' ','photolysis diagnostic organic carbon OD' )
    call addfld( 'DTSO4',   horiz_only, 'I', ' ','photolysis diagnostic SO4 OD' )
    call addfld( 'DTSOA',   horiz_only, 'I', ' ','photolysis diagnostic SOA OD' )
    call addfld( 'DTANT',   horiz_only, 'I', ' ','photolysis diagnostic NH4SO4 OD' )
    call addfld( 'DTSAL',   horiz_only, 'I', ' ','photolysis diagnostic salt OD' )
    call addfld( 'DTDUST',  horiz_only, 'I', ' ','photolysis diagnostic dust OD' )
    call addfld( 'DTTOTAL', horiz_only, 'I', ' ','photolysis diagnostic total aerosol OD' )   
    call addfld( 'FRACDAY', horiz_only, 'I', ' ','photolysis diagnostic fraction of day' )

    call addfld( 'QDSAD',      (/ 'lev' /), 'I', '/s',      'water vapor sad delta' )
    call addfld( 'SAD_STRAT',  (/ 'lev' /), 'I', 'cm2/cm3', 'stratospheric aerosol SAD' )
    call addfld( 'SAD_SULFC',  (/ 'lev' /), 'I', 'cm2/cm3', 'chemical sulfate aerosol SAD' )
    call addfld( 'SAD_SAGE',   (/ 'lev' /), 'I', 'cm2/cm3', 'SAGE sulfate aerosol SAD' )
    call addfld( 'SAD_LNAT',   (/ 'lev' /), 'I', 'cm2/cm3', 'large-mode NAT aerosol SAD' )
    call addfld( 'SAD_ICE',    (/ 'lev' /), 'I', 'cm2/cm3', 'water-ice aerosol SAD' )
    call addfld( 'RAD_SULFC',  (/ 'lev' /), 'I', 'cm',      'chemical sad sulfate' )
    call addfld( 'RAD_LNAT',   (/ 'lev' /), 'I', 'cm',      'large nat radius' )
    call addfld( 'RAD_ICE',    (/ 'lev' /), 'I', 'cm',      'sad ice' )
    call addfld( 'SAD_TROP',   (/ 'lev' /), 'I', 'cm2/cm3', 'tropospheric aerosol SAD' )
    call addfld( 'SAD_AERO',   (/ 'lev' /), 'I', 'cm2/cm3', 'aerosol surface area density' )
    if (history_cesm_forcing) then
       call add_default ('SAD_AERO',8,' ')
    endif
    call addfld( 'REFF_AERO',  (/ 'lev' /), 'I', 'cm',      'aerosol effective radius' )
    call addfld( 'SULF_TROP',  (/ 'lev' /), 'I', 'mol/mol', 'tropospheric aerosol SAD' )
    call addfld( 'QDSETT',     (/ 'lev' /), 'I', '/s',      'water vapor settling delta' )
    call addfld( 'QDCHEM',     (/ 'lev' /), 'I', '/s',      'water vapor chemistry delta')
    call addfld( 'HNO3_TOTAL', (/ 'lev' /), 'I', 'mol/mol', 'total HNO3' )
    call addfld( 'HNO3_STS',   (/ 'lev' /), 'I', 'mol/mol', 'STS condensed HNO3' )
    call addfld( 'HNO3_NAT',   (/ 'lev' /), 'I', 'mol/mol', 'NAT condensed HNO3' )
    call addfld( 'HNO3_GAS',   (/ 'lev' /), 'I', 'mol/mol', 'gas-phase hno3' )
    call addfld( 'H2O_GAS',    (/ 'lev' /), 'I', 'mol/mol', 'gas-phase h2o' )
    call addfld( 'HCL_TOTAL',  (/ 'lev' /), 'I', 'mol/mol', 'total hcl' )
    call addfld( 'HCL_GAS',    (/ 'lev' /), 'I', 'mol/mol', 'gas-phase hcl' )
    call addfld( 'HCL_STS',    (/ 'lev' /), 'I', 'mol/mol', 'STS condensed HCL' )

    if (het1_ndx>0) then
       call addfld( 'het1_total', (/ 'lev' /), 'I', '/s', 'total N2O5 + H2O het rate constant' )
    endif
    call addfld( 'SZA', horiz_only, 'I', 'degrees', 'solar zenith angle' )

    call chm_diags_inti()
    call rate_diags_init()

!-----------------------------------------------------------------------
! get pbuf indicies
!-----------------------------------------------------------------------
    ndx_cldfr  = pbuf_get_index('CLD')
    ndx_cmfdqr = pbuf_get_index('RPRDTOT')
    ndx_nevapr = pbuf_get_index('NEVAPR')
    ndx_prain  = pbuf_get_index('PRAIN')
    ndx_cldtop = pbuf_get_index('CLDTOP')

    sad_pbf_ndx= pbuf_get_index('VOLC_SAD',errcode=err) ! prescribed  strat aerosols (volcanic)
    if (.not.sad_pbf_ndx>0) sad_pbf_ndx = pbuf_get_index('SADSULF',errcode=err) ! CARMA's version of strat aerosols

    ele_temp_ndx = pbuf_get_index('TElec',errcode=err)! electron temperature index 
    ion_temp_ndx = pbuf_get_index('TIon',errcode=err) ! ion temperature index

    ! diagnostics for stratospheric heterogeneous reactions
    call addfld( 'GAMMA_HET1', (/ 'lev' /), 'I', '1', 'Reaction Probability' )
    call addfld( 'GAMMA_HET2', (/ 'lev' /), 'I', '1', 'Reaction Probability' )
    call addfld( 'GAMMA_HET3', (/ 'lev' /), 'I', '1', 'Reaction Probability' )
    call addfld( 'GAMMA_HET4', (/ 'lev' /), 'I', '1', 'Reaction Probability' )
    call addfld( 'GAMMA_HET5', (/ 'lev' /), 'I', '1', 'Reaction Probability' )
    call addfld( 'GAMMA_HET6', (/ 'lev' /), 'I', '1', 'Reaction Probability' )
    call addfld( 'WTPER',      (/ 'lev' /), 'I', '%', 'H2SO4 Weight Percent' )

    call addfld( 'O3S_LOSS',      (/ 'lev' /), 'I', '1/sec', 'O3S loss rate const' )

    call chem_prod_loss_diags_init

  end subroutine gas_phase_chemdr_inti


!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
  subroutine gas_phase_chemdr(lchnk, ncol, imozart, q, &
                              phis, zm, zi, calday, &
                              tfld, pmid, pdel, pint,  &
                              cldw, troplev, troplevchem, &
                              ncldwtr, ufld, vfld,  &
                              delt, ps, xactive_prates, &
                              fsds, ts, asdir, ocnfrac, icefrac, &
                              precc, precl, snowhland, ghg_chem, latmapback, &
                              drydepflx, wetdepflx, cflx, fire_sflx, fire_ztop, nhx_nitrogen_flx, noy_nitrogen_flx, qtend, pbuf)

    !-----------------------------------------------------------------------
    !     ... Chem_solver advances the volumetric mixing ratio
    !         forward one time step via a combination of explicit,
    !         ebi, hov, fully implicit, and/or rodas algorithms.
    !-----------------------------------------------------------------------

    use chem_mods,         only : nabscol, nfs, indexm, clscnt4
    use physconst,         only : rga
    use mo_photo,          only : set_ub_col, setcol, table_photo, xactive_photo
    use mo_exp_sol,        only : exp_sol
    use mo_imp_sol,        only : imp_sol
    use mo_setrxt,         only : setrxt
    use mo_adjrxt,         only : adjrxt
    use mo_phtadj,         only : phtadj
    use llnl_O1D_to_2OH_adj,only : O1D_to_2OH_adj
    use mo_usrrxt,         only : usrrxt
    use mo_setinv,         only : setinv
    use mo_negtrc,         only : negtrc
    use mo_sulf,           only : sulf_interp
    use mo_setext,         only : setext
    use fire_emissions,    only : fire_emissions_vrt
    use mo_sethet,         only : sethet
    use mo_drydep,         only : drydep, set_soilw
    use seq_drydep_mod,    only : DD_XLND, DD_XATM, DD_TABL, drydep_method
    use mo_fstrat,         only : set_fstrat_vals, set_fstrat_h2o
    use noy_ubc,           only : noy_ubc_set
    use mo_flbc,           only : flbc_set
    use phys_grid,         only : get_rlat_all_p, get_rlon_all_p, get_lat_all_p, get_lon_all_p
    use mo_mean_mass,      only : set_mean_mass
    use cam_history,       only : outfld
    use wv_saturation,     only : qsat
    use constituents,      only : cnst_mw
    use mo_drydep,         only : has_drydep
    use time_manager,      only : get_ref_date
    use mo_ghg_chem,       only : ghg_chem_set_rates, ghg_chem_set_flbc
    use mo_sad,            only : sad_strat_calc
    use charge_neutrality, only : charge_balance
    use mo_strato_rates,   only : ratecon_sfstrat
    use mo_aero_settling,  only : strat_aer_settling
    use shr_orb_mod,       only : shr_orb_decl
    use cam_control_mod,   only : lambm0, eccen, mvelpp, obliqr
    use mo_strato_rates,   only : has_strato_chem
    use short_lived_species,only: set_short_lived_species,get_short_lived_species
    use mo_chm_diags,      only : chm_diags, het_diags
    use perf_mod,          only : t_startf, t_stopf
    use gas_wetdep_opts,   only : gas_wetdep_method
    use physics_buffer,    only : physics_buffer_desc, pbuf_get_field, pbuf_old_tim_idx
    use infnan,            only : nan, assignment(=)
    use rate_diags,        only : rate_diags_calc, rate_diags_o3s_loss
    use mo_mass_xforms,    only : mmr2vmr, vmr2mmr, h2o_to_vmr, mmr2vmri
    use orbit,             only : zenith
!
! LINOZ
!
    use lin_strat_chem,    only : do_lin_strat_chem, lin_strat_chem_solve
    use linoz_data,        only : has_linoz_data
!
! for aqueous chemistry and aerosol growth
!
    use aero_model,        only : aero_model_gasaerexch

    use aero_model,        only : aero_model_strat_surfarea

    implicit none

    !-----------------------------------------------------------------------
    !        ... Dummy arguments
    !-----------------------------------------------------------------------
    integer,        intent(in)    :: lchnk                          ! chunk index
    integer,        intent(in)    :: ncol                           ! number columns in chunk
    integer,        intent(in)    :: imozart                        ! gas phase start index in q
    real(r8),       intent(in)    :: delt                           ! timestep (s)
    real(r8),       intent(in)    :: calday                         ! day of year
    real(r8),       intent(in)    :: ps(pcols)                      ! surface pressure
    real(r8),       intent(in)    :: phis(pcols)                    ! surface geopotential
    real(r8),target,intent(in)    :: tfld(pcols,pver)               ! midpoint temperature (K)
    real(r8),       intent(in)    :: pmid(pcols,pver)               ! midpoint pressures (Pa)
    real(r8),       intent(in)    :: pdel(pcols,pver)               ! pressure delta about midpoints (Pa)
    real(r8),       intent(in)    :: ufld(pcols,pver)               ! zonal velocity (m/s)
    real(r8),       intent(in)    :: vfld(pcols,pver)               ! meridional velocity (m/s)
    real(r8),       intent(in)    :: cldw(pcols,pver)               ! cloud water (kg/kg)
    real(r8),       intent(in)    :: ncldwtr(pcols,pver)            ! droplet number concentration (#/kg)
    real(r8),       intent(in)    :: zm(pcols,pver)                 ! midpoint geopotential height above the surface (m)
    real(r8),       intent(in)    :: zi(pcols,pver+1)               ! interface geopotential height above the surface (m)
    real(r8),       intent(in)    :: pint(pcols,pver+1)             ! interface pressures (Pa)
    real(r8),       intent(in)    :: q(pcols,pver,pcnst)            ! species concentrations (kg/kg)
    real(r8),pointer, intent(in)  :: fire_sflx(:,:)                 ! fire emssions surface flux (kg/m^2/s)
    real(r8),pointer, intent(in)  :: fire_ztop(:)                   ! top of vertical distribution of fire emssions (m)
    logical,        intent(in)    :: xactive_prates
    real(r8),       intent(in)    :: fsds(pcols)                    ! longwave down at sfc
    real(r8),       intent(in)    :: icefrac(pcols)                 ! sea-ice areal fraction
    real(r8),       intent(in)    :: ocnfrac(pcols)                 ! ocean areal fraction
    real(r8),       intent(in)    :: asdir(pcols)                   ! albedo: shortwave, direct
    real(r8),       intent(in)    :: ts(pcols)                      ! sfc temp (merged w/ocean if coupled)
    real(r8),       intent(in)    :: precc(pcols)                   !
    real(r8),       intent(in)    :: precl(pcols)                   !
    real(r8),       intent(in)    :: snowhland(pcols)               !
    logical,        intent(in)    :: ghg_chem 
    integer,        intent(in)    :: latmapback(pcols)
    integer,        intent(in)    :: troplev(pcols)                 ! trop/strat separation vertical index
    integer,        intent(in)    :: troplevchem(pcols)             ! trop/strat chemistry separation vertical index
    real(r8),       intent(inout) :: qtend(pcols,pver,pcnst)        ! species tendencies (kg/kg/s)
    real(r8),       intent(inout) :: cflx(pcols,pcnst)              ! constituent surface flux (kg/m^2/s)
    real(r8),       intent(out)   :: drydepflx(pcols,pcnst)         ! dry deposition flux (kg/m^2/s)
    real(r8),       intent(in)    :: wetdepflx(pcols,pcnst)         ! wet deposition flux (kg/m^2/s)
    real(r8), intent(out) :: nhx_nitrogen_flx(pcols)
    real(r8), intent(out) :: noy_nitrogen_flx(pcols)

    type(physics_buffer_desc), pointer :: pbuf(:)

    !-----------------------------------------------------------------------
    !     	... Local variables
    !-----------------------------------------------------------------------
    real(r8), parameter :: m2km  = 1.e-3_r8
    real(r8), parameter :: Pa2mb = 1.e-2_r8

    real(r8),       pointer    :: prain(:,:)
    real(r8),       pointer    :: nevapr(:,:)
    real(r8),       pointer    :: cmfdqr(:,:)
    real(r8),       pointer    :: cldfr(:,:)
    real(r8),       pointer    :: cldtop(:)

    integer      ::  i, k, m, n
    integer      ::  tim_ndx
    real(r8)     ::  delt_inverse
    real(r8)     ::  esfact
    integer      ::  latndx(pcols)                         ! chunk lat indicies
    integer      ::  lonndx(pcols)                         ! chunk lon indicies
    real(r8)     ::  invariants(ncol,pver,nfs)
    real(r8)     ::  col_dens(ncol,pver,max(1,nabscol))    ! column densities (molecules/cm^2)
    real(r8)     ::  col_delta(ncol,0:pver,max(1,nabscol)) ! layer column densities (molecules/cm^2)
    real(r8)     ::  extfrc(ncol,pver,max(1,extcnt))
    real(r8)     ::  vmr(ncol,pver,gas_pcnst)         ! xported species (vmr)
    real(r8)     ::  reaction_rates(ncol,pver,max(1,rxntot))      ! reaction rates
    real(r8)     ::  depvel(ncol,gas_pcnst)                ! dry deposition velocity (cm/s)
    real(r8)     ::  het_rates(ncol,pver,max(1,gas_pcnst)) ! washout rate (1/s)
    real(r8), dimension(ncol,pver) :: &
         h2ovmr, &                                         ! water vapor volume mixing ratio
         mbar, &                                           ! mean wet atmospheric mass ( amu )
         zmid, &                                           ! midpoint geopotential in km
         zmidr, &                                          ! midpoint geopotential in km realitive to surf
         sulfate, &                                        ! trop sulfate aerosols
         pmb                                               ! pressure at midpoints ( hPa )
    real(r8), dimension(ncol,pver) :: &
         cwat, &                                           ! cloud water mass mixing ratio (kg/kg)
         wrk
    real(r8), dimension(ncol,pver+1) :: &
         zintr                                              ! interface geopotential in km realitive to surf
    real(r8), dimension(ncol,pver+1) :: &
         zint                                              ! interface geopotential in km
    real(r8), dimension(ncol)  :: &
         zen_angle, &                                      ! solar zenith angles
         zsurf, &                                          ! surface height (m)
         rlats, rlons                                      ! chunk latitudes and longitudes (radians)
    real(r8) :: sza(ncol)                                  ! solar zenith angles (degrees)
    real(r8), parameter :: rad2deg = 180._r8/pi                ! radians to degrees conversion factor
    real(r8) :: relhum(ncol,pver)                          ! relative humidity
    real(r8) :: satv(ncol,pver)                            ! wrk array for relative humidity
    real(r8) :: satq(ncol,pver)                            ! wrk array for relative humidity

    integer                   :: j
    integer                   ::  ltrop_sol(pcols)         ! tropopause vertical index used in chem solvers
    real(r8), pointer         ::  strato_sad(:,:)          ! stratospheric sad (1/cm)

    real(r8)                  ::  sad_trop(pcols,pver)    ! total tropospheric sad (cm^2/cm^3)
    real(r8)                  ::  reff(pcols,pver)        ! aerosol effective radius (cm)
    real(r8)                  ::  reff_strat(pcols,pver)  ! stratospheric aerosol effective radius (cm)

    real(r8) :: tvs(pcols)
    integer  :: ncdate,yr,mon,day,sec
    real(r8) :: wind_speed(pcols)        ! surface wind speed (m/s)
    logical, parameter :: dyn_soilw = .false.
    logical  :: table_soilw
    real(r8) :: soilw(pcols)
    real(r8) :: prect(pcols)
    real(r8) :: sflx(pcols,gas_pcnst)
    real(r8) :: wetdepflx_diag(pcols,gas_pcnst)
    real(r8) :: dust_vmr(ncol,pver,ndust)
    real(r8) :: dt_diag(pcols,8)               ! od diagnostics
    real(r8) :: fracday(pcols)                 ! fraction of day
    real(r8) :: o2mmr(ncol,pver)               ! o2 concentration (kg/kg)
    real(r8) :: ommr(ncol,pver)                ! o concentration (kg/kg)
    real(r8) :: mmr(pcols,pver,gas_pcnst)      ! chem working concentrations (kg/kg)
    real(r8) :: mmr_new(pcols,pver,gas_pcnst)      ! chem working concentrations (kg/kg)
    real(r8) :: hno3_gas(ncol,pver)            ! hno3 gas phase concentration (mol/mol)
    real(r8) :: hno3_cond(ncol,pver,2)         ! hno3 condensed phase concentration (mol/mol)
    real(r8) :: hcl_gas(ncol,pver)             ! hcl gas phase concentration (mol/mol)
    real(r8) :: hcl_cond(ncol,pver)            ! hcl condensed phase concentration (mol/mol)
    real(r8) :: h2o_gas(ncol,pver)             ! h2o gas phase concentration (mol/mol)
    real(r8) :: h2o_cond(ncol,pver)            ! h2o condensed phase concentration (mol/mol)
    real(r8) :: cldice(pcols,pver)             ! cloud water "ice" (kg/kg)
    real(r8) :: radius_strat(ncol,pver,3)      ! radius of sulfate, nat, & ice ( cm )
    real(r8) :: sad_strat(ncol,pver,3)         ! surf area density of sulfate, nat, & ice ( cm^2/cm^3 )
    real(r8) :: mmr_tend(pcols,pver,gas_pcnst) ! chemistry species tendencies (kg/kg/s)
    real(r8) :: qh2o(pcols,pver)               ! specific humidity (kg/kg)
    real(r8) :: delta

  ! for aerosol formation....  
    real(r8) :: del_h2so4_gasprod(ncol,pver)
    real(r8) :: vmr0(ncol,pver,gas_pcnst)

!
! CCMI
!
    real(r8) :: xlat
    real(r8) :: pm25(ncol)

    real(r8) :: dlats(ncol)

    real(r8), dimension(ncol,pver) :: &      ! aerosol reaction diagnostics
        gprob_n2o5, &
        gprob_cnt_hcl, &
        gprob_cnt_h2o, &
        gprob_bnt_h2o, &
        gprob_hocl_hcl, &
        gprob_hobr_hcl, &
        wtper

    real(r8), pointer :: ele_temp_fld(:,:) ! electron temperature pointer
    real(r8), pointer :: ion_temp_fld(:,:) ! ion temperature pointer
    real(r8) :: prod_out(ncol,pver,max(1,clscnt4))
    real(r8) :: loss_out(ncol,pver,max(1,clscnt4))

    real(r8) :: o3s_loss(ncol,pver)

    if ( ele_temp_ndx>0 .and. ion_temp_ndx>0 ) then
       call pbuf_get_field(pbuf, ele_temp_ndx, ele_temp_fld)
       call pbuf_get_field(pbuf, ion_temp_ndx, ion_temp_fld)
    else
       ele_temp_fld => tfld
       ion_temp_fld => tfld
    endif

    ! initialize to NaN to hopefully catch user defined rxts that go unset
    reaction_rates(:,:,:) = nan

    delt_inverse = 1._r8 / delt
    !-----------------------------------------------------------------------      
    !        ... Get chunck latitudes and longitudes
    !-----------------------------------------------------------------------      
    call get_lat_all_p( lchnk, ncol, latndx )
    call get_lon_all_p( lchnk, ncol, lonndx )
    call get_rlat_all_p( lchnk, ncol, rlats )
    call get_rlon_all_p( lchnk, ncol, rlons )
    tim_ndx = pbuf_old_tim_idx()
    call pbuf_get_field(pbuf, ndx_prain,      prain,  start=(/1,1/), kount=(/ncol,pver/))
    call pbuf_get_field(pbuf, ndx_cldfr,        cldfr, start=(/1,1,tim_ndx/), kount=(/ncol,pver,1/) )
    call pbuf_get_field(pbuf, ndx_cmfdqr,     cmfdqr, start=(/1,1/), kount=(/ncol,pver/))
    call pbuf_get_field(pbuf, ndx_nevapr,     nevapr, start=(/1,1/), kount=(/ncol,pver/))
    call pbuf_get_field(pbuf, ndx_cldtop,     cldtop )

    reff_strat(:,:) = 0._r8

    dlats(:) = rlats(:)*rad2deg ! convert to degrees

    !-----------------------------------------------------------------------      
    !        ... Calculate cosine of zenith angle
    !            then cast back to angle (radians)
    !-----------------------------------------------------------------------      
    call zenith( calday, rlats, rlons, zen_angle, ncol )
    zen_angle(:) = acos( zen_angle(:) )

    sza(:) = zen_angle(:) * rad2deg
    call outfld( 'SZA',   sza,    ncol, lchnk )

    !-----------------------------------------------------------------------      
    !        ... Xform geopotential height from m to km 
    !            and pressure from Pa to mb
    !-----------------------------------------------------------------------      
    zsurf(:ncol) = rga * phis(:ncol)
    do k = 1,pver
       zintr(:ncol,k) = m2km * zi(:ncol,k)
       zmidr(:ncol,k) = m2km * zm(:ncol,k)
       zmid(:ncol,k) = m2km * (zm(:ncol,k) + zsurf(:ncol))
       zint(:ncol,k) = m2km * (zi(:ncol,k) + zsurf(:ncol))
       pmb(:ncol,k)  = Pa2mb * pmid(:ncol,k)
    end do
    zint(:ncol,pver+1) = m2km * (zi(:ncol,pver+1) + zsurf(:ncol))
    zintr(:ncol,pver+1)= m2km *  zi(:ncol,pver+1)

    !-----------------------------------------------------------------------      
    !        ... map incoming concentrations to working array
    !-----------------------------------------------------------------------      
    do m = 1,pcnst
       n = map2chm(m)
       if( n > 0 ) then
          mmr(:ncol,:,n) = q(:ncol,:,m)
       end if
    end do

    call get_short_lived_species( mmr, lchnk, ncol, pbuf )

    !-----------------------------------------------------------------------      
    !        ... Set atmosphere mean mass
    !-----------------------------------------------------------------------      
    call set_mean_mass( ncol, mmr, mbar )

    !-----------------------------------------------------------------------      
    !        ... Xform from mmr to vmr
    !-----------------------------------------------------------------------      
    call mmr2vmr( mmr(:ncol,:,:), vmr(:ncol,:,:), mbar(:ncol,:), ncol )
    
!
! CCMI
!
! reset STE tracer to specific vmr of 200 ppbv
!
    if ( st80_25_ndx > 0 ) then 
       where ( pmid(:ncol,:) < 80.e+2_r8 )
          vmr(:ncol,:,st80_25_ndx) = 200.e-9_r8 
       end where
    end if
!
! reset AOA_NH, NH_5, NH_50, NH_50W surface mixing ratios between 30N and 50N
!
    if ( aoa_nh_ndx>0 ) then
      do j=1,ncol
        xlat = dlats(j)
        if ( xlat >= 30._r8 .and. xlat <= 50._r8 ) then
           vmr(j,pver,aoa_nh_ndx) = 0._r8
        end if
      end do
    end if
    if ( nh_5_ndx>0 ) then
      do j=1,ncol
        xlat = dlats(j)
        if ( xlat >= 30._r8 .and. xlat <= 50._r8 ) then
           vmr(j,pver,nh_5_ndx)   = 100.e-9_r8
        end if
      end do
    end if
    if ( nh_50_ndx>0 ) then
      do j=1,ncol
        xlat = dlats(j)
        if ( xlat >= 30._r8 .and. xlat <= 50._r8 ) then
           vmr(j,pver,nh_50_ndx)  = 100.e-9_r8
        end if
      end do
    end if
    if ( nh_50w_ndx>0 ) then
      do j=1,ncol
        xlat = dlats(j)
        if ( xlat >= 30._r8 .and. xlat <= 50._r8 ) then
           vmr(j,pver,nh_50w_ndx) = 100.e-9_r8
        end if
      end do
    end if

    if (h2o_ndx>0) then
       !-----------------------------------------------------------------------      
       !        ... store water vapor in wrk variable
       !-----------------------------------------------------------------------      
       qh2o(:ncol,:) = mmr(:ncol,:,h2o_ndx)
       h2ovmr(:ncol,:) = vmr(:ncol,:,h2o_ndx)
    else
       qh2o(:ncol,:) = q(:ncol,:,1)
       !-----------------------------------------------------------------------      
       !        ... Xform water vapor from mmr to vmr and set upper bndy values
       !-----------------------------------------------------------------------      
       call h2o_to_vmr( q(:ncol,:,1), h2ovmr(:ncol,:), mbar(:ncol,:), ncol )

       call set_fstrat_h2o( h2ovmr, pmid, troplev, calday, ncol, lchnk )

    endif

    !-----------------------------------------------------------------------      
    !        ... force ion/electron balance
    !-----------------------------------------------------------------------      
    call charge_balance( ncol, vmr )

    !-----------------------------------------------------------------------      
    !        ... Set the "invariants"
    !-----------------------------------------------------------------------  
    call setinv( invariants, tfld, h2ovmr, vmr, pmid, ncol, lchnk, pbuf )

    !-----------------------------------------------------------------------      
    !        ... stratosphere aerosol surface area
    !-----------------------------------------------------------------------  
    if (sad_pbf_ndx>0) then
       call pbuf_get_field(pbuf, sad_pbf_ndx, strato_sad)
    else
       allocate(strato_sad(pcols,pver))
       strato_sad(:,:) = 0._r8

       ! Prognostic modal stratospheric sulfate: compute dry strato_sad
       call aero_model_strat_surfarea( ncol, mmr, pmid, tfld, troplevchem, pbuf, strato_sad, reff_strat )

    endif

    stratochem: if ( has_strato_chem ) then
       !-----------------------------------------------------------------------      
       !        ... initialize condensed and gas phases; all hno3 to gas
       !-----------------------------------------------------------------------    
       hcl_cond(:,:)      = 0.0_r8
       hcl_gas (:,:)      = 0.0_r8  
       do k = 1,pver
          hno3_gas(:,k)   = vmr(:,k,hno3_ndx)
          h2o_gas(:,k)    = h2ovmr(:,k)
          hcl_gas(:,k)    = vmr(:,k,hcl_ndx)
          wrk(:,k)        = h2ovmr(:,k)
          if (snow_ndx>0) then
             cldice(:ncol,k) = q(:ncol,k,cldice_ndx) + q(:ncol,k,snow_ndx)
          else
             cldice(:ncol,k) = q(:ncol,k,cldice_ndx)
          endif
       end do
       do m = 1,2
          do k = 1,pver
             hno3_cond(:,k,m) = 0._r8
          end do
       end do

       call mmr2vmri( cldice(:ncol,:), h2o_cond(:ncol,:), mbar(:ncol,:), cnst_mw(cldice_ndx), ncol )

       !-----------------------------------------------------------------------      
       !        ... call SAD routine
       !-----------------------------------------------------------------------      
       call sad_strat_calc( lchnk, invariants(:ncol,:,indexm), pmb, tfld, hno3_gas, &
            hno3_cond, h2o_gas, h2o_cond, hcl_gas, hcl_cond, strato_sad(:ncol,:), radius_strat, &
            sad_strat, ncol, pbuf )

!      NOTE: output of total HNO3 is before vmr is set to gas-phase.
       call outfld( 'HNO3_TOTAL', vmr(:ncol,:,hno3_ndx), ncol ,lchnk )


       do k = 1,pver
          vmr(:,k,hno3_ndx) = hno3_gas(:,k)
          h2ovmr(:,k)       = h2o_gas(:,k)
          vmr(:,k,h2o_ndx)  = h2o_gas(:,k)
          wrk(:,k)          = (h2ovmr(:,k) - wrk(:,k))*delt_inverse
       end do

       call outfld( 'QDSAD', wrk(:,:), ncol, lchnk )
!
       call outfld( 'SAD_STRAT',  strato_sad(:ncol,:), ncol, lchnk )
       call outfld( 'SAD_SULFC',  sad_strat(:,:,1), ncol, lchnk )
       call outfld( 'SAD_LNAT',   sad_strat(:,:,2), ncol, lchnk )
       call outfld( 'SAD_ICE',    sad_strat(:,:,3), ncol, lchnk )
!
       call outfld( 'RAD_SULFC',  radius_strat(:,:,1), ncol, lchnk )
       call outfld( 'RAD_LNAT',   radius_strat(:,:,2), ncol, lchnk )
       call outfld( 'RAD_ICE',    radius_strat(:,:,3), ncol, lchnk )
!
       call outfld( 'HNO3_GAS',   vmr(:ncol,:,hno3_ndx), ncol, lchnk )
       call outfld( 'HNO3_STS',   hno3_cond(:,:,1), ncol, lchnk )
       call outfld( 'HNO3_NAT',   hno3_cond(:,:,2), ncol, lchnk )
!
       call outfld( 'HCL_TOTAL',  vmr(:ncol,:,hcl_ndx), ncol, lchnk )
       call outfld( 'HCL_GAS',    hcl_gas (:,:), ncol ,lchnk )
       call outfld( 'HCL_STS',    hcl_cond(:,:), ncol ,lchnk )

       !-----------------------------------------------------------------------      
       !        ... call aerosol reaction rates
       !-----------------------------------------------------------------------      
       call ratecon_sfstrat( ncol, invariants(:,:,indexm), pmid, tfld, &
            radius_strat(:,:,1), sad_strat(:,:,1), sad_strat(:,:,2), &
            sad_strat(:,:,3), h2ovmr, vmr, reaction_rates, &
            gprob_n2o5, gprob_cnt_hcl, gprob_cnt_h2o, gprob_bnt_h2o, &
            gprob_hocl_hcl, gprob_hobr_hcl, wtper )

       call outfld( 'GAMMA_HET1', gprob_n2o5    (:ncol,:), ncol, lchnk )
       call outfld( 'GAMMA_HET2', gprob_cnt_h2o (:ncol,:), ncol, lchnk )
       call outfld( 'GAMMA_HET3', gprob_bnt_h2o (:ncol,:), ncol, lchnk )
       call outfld( 'GAMMA_HET4', gprob_cnt_hcl (:ncol,:), ncol, lchnk )
       call outfld( 'GAMMA_HET5', gprob_hocl_hcl(:ncol,:), ncol, lchnk )
       call outfld( 'GAMMA_HET6', gprob_hobr_hcl(:ncol,:), ncol, lchnk )
       call outfld( 'WTPER',      wtper         (:ncol,:), ncol, lchnk )

    endif stratochem

!      NOTE: For gas-phase solver only. 
!            ratecon_sfstrat needs total hcl.
    if (hcl_ndx>0) then
       vmr(:,:,hcl_ndx)  = hcl_gas(:,:)
    endif

    !-----------------------------------------------------------------------      
    !        ... Set the column densities at the upper boundary
    !-----------------------------------------------------------------------      
    call set_ub_col( col_delta, vmr, invariants, pint(:,1), pdel, ncol, lchnk)

    !-----------------------------------------------------------------------      
    !       ...  Set rates for "tabular" and user specified reactions
    !-----------------------------------------------------------------------      
    call setrxt( reaction_rates, tfld, invariants(1,1,indexm), ncol )
    
    sulfate(:,:) = 0._r8
    if ( .not. carma_hetchem_feedback ) then
       if( so4_ndx < 1 ) then ! get offline so4 field if not prognostic
          call sulf_interp( ncol, lchnk, sulfate )
       else
          sulfate(:,:) = vmr(:,:,so4_ndx)
       endif
    endif
    
    !-----------------------------------------------------------------
    ! ... zero out sulfate above tropopause
    !-----------------------------------------------------------------
    do k = 1, pver
       do i = 1, ncol
          if (k < troplevchem(i)) then
             sulfate(i,k) = 0.0_r8
          end if
       end do
    end do

    call outfld( 'SULF_TROP', sulfate(:ncol,:), ncol, lchnk )

    !-----------------------------------------------------------------
    !	... compute the relative humidity
    !-----------------------------------------------------------------
    call qsat(tfld(:ncol,:), pmid(:ncol,:), satv, satq)

    do k = 1,pver
       relhum(:,k) = .622_r8 * h2ovmr(:,k) / satq(:,k)
       relhum(:,k) = max( 0._r8,min( 1._r8,relhum(:,k) ) )
    end do
    
    cwat(:ncol,:pver) = cldw(:ncol,:pver)

    call usrrxt( reaction_rates, tfld, ion_temp_fld, ele_temp_fld, invariants, h2ovmr, &
                 pmid, invariants(:,:,indexm), sulfate, mmr, relhum, strato_sad, &
                 troplevchem, dlats, ncol, sad_trop, reff, cwat, mbar, pbuf )

    call outfld( 'SAD_TROP', sad_trop(:ncol,:), ncol, lchnk )

    ! Add trop/strat components of SAD for output
    sad_trop(:ncol,:)=sad_trop(:ncol,:)+strato_sad(:ncol,:)
    call outfld( 'SAD_AERO', sad_trop(:ncol,:), ncol, lchnk )

    ! Add trop/strat components of effective radius for output
    reff(:ncol,:)=reff(:ncol,:)+reff_strat(:ncol,:)
    call outfld( 'REFF_AERO', reff(:ncol,:), ncol, lchnk )
    
    if (het1_ndx>0) then
       call outfld( 'het1_total', reaction_rates(:,:,het1_ndx), ncol, lchnk )
    endif

    if (ghg_chem) then
       call ghg_chem_set_rates( reaction_rates, latmapback, zen_angle, ncol, lchnk )
    endif

    do i = phtcnt+1,rxt_tag_cnt
       call outfld( tag_names(i), reaction_rates(:ncol,:,rxt_tag_map(i)), ncol, lchnk )
    enddo

    call adjrxt( reaction_rates, invariants, invariants(1,1,indexm), ncol,pver )

    !-----------------------------------------------------------------------
    !        ... Compute the photolysis rates at time = t(n+1)
    !-----------------------------------------------------------------------      
    !-----------------------------------------------------------------------      
    !     	... Set the column densities
    !-----------------------------------------------------------------------      
    call setcol( col_delta, col_dens, vmr, pdel,  ncol )

    !-----------------------------------------------------------------------      
    !     	... Calculate the photodissociation rates
    !-----------------------------------------------------------------------      

    esfact = 1._r8
    call shr_orb_decl( calday, eccen, mvelpp, lambm0, obliqr  , &
         delta, esfact )


    if ( xactive_prates ) then
       if ( dst_ndx > 0 ) then
          dust_vmr(:ncol,:,1:ndust) = vmr(:ncol,:,dst_ndx:dst_ndx+ndust-1)
       else 
          dust_vmr(:ncol,:,:) = 0._r8
       endif

       !-----------------------------------------------------------------
       !	... compute the photolysis rates
       !-----------------------------------------------------------------
       call xactive_photo( reaction_rates, vmr, tfld, cwat, cldfr, &
            pmid, zmidr, col_dens, zen_angle, asdir, &
            invariants(1,1,indexm), ps, ts, &
            esfact, relhum, dust_vmr, dt_diag, fracday, ncol, lchnk )

       call outfld('DTCBS',   dt_diag(:ncol,1), ncol, lchnk )
       call outfld('DTOCS',   dt_diag(:ncol,2), ncol, lchnk )
       call outfld('DTSO4',   dt_diag(:ncol,3), ncol, lchnk )
       call outfld('DTANT',   dt_diag(:ncol,4), ncol, lchnk )
       call outfld('DTSAL',   dt_diag(:ncol,5), ncol, lchnk )
       call outfld('DTDUST',  dt_diag(:ncol,6), ncol, lchnk )
       call outfld('DTSOA',   dt_diag(:ncol,7), ncol, lchnk )
       call outfld('DTTOTAL', dt_diag(:ncol,8), ncol, lchnk )
       call outfld('FRACDAY', fracday(:ncol), ncol, lchnk )

    else
       !-----------------------------------------------------------------
       !	... lookup the photolysis rates from table
       !-----------------------------------------------------------------
       call table_photo( reaction_rates, pmid, pdel, tfld, zmid, zint, &
                         col_dens, zen_angle, asdir, cwat, cldfr, &
                         esfact, vmr, invariants, ncol, lchnk, pbuf )
    endif

    do i = 1,phtcnt
       call outfld( tag_names(i), reaction_rates(:ncol,:,rxt_tag_map(i)), ncol, lchnk )
    enddo

    !-----------------------------------------------------------------------      
    !     	... Adjust the photodissociation rates
    !-----------------------------------------------------------------------  
    call O1D_to_2OH_adj( reaction_rates, invariants, invariants(:,:,indexm), ncol, tfld )
    call phtadj( reaction_rates, invariants, invariants(:,:,indexm), ncol,pver )

    !-----------------------------------------------------------------------
    !        ... Compute the extraneous frcing at time = t(n+1)
    !-----------------------------------------------------------------------    
    if ( o2_ndx > 0 .and. o_ndx > 0 ) then
       do k = 1,pver
          o2mmr(:ncol,k) = mmr(:ncol,k,o2_ndx)
          ommr(:ncol,k)  = mmr(:ncol,k,o_ndx)
       end do
    endif
    call setext( extfrc, zint, zintr, cldtop, &
                 zmid, lchnk, tfld, o2mmr, ommr, &
                 pmid, mbar, rlats, calday, ncol, rlons, pbuf )
    ! include forcings from fire emissions ...
    call fire_emissions_vrt( ncol, lchnk, zint, fire_sflx, fire_ztop, extfrc )

    do m = 1,extcnt
       if( m /= synoz_ndx .and. m /= aoa_nh_ext_ndx ) then
          do k = 1,pver
             extfrc(:ncol,k,m) = extfrc(:ncol,k,m) / invariants(:ncol,k,indexm)
          end do
       endif
       call outfld( extfrc_name(m), extfrc(:ncol,:,m), ncol, lchnk )
    end do

    !-----------------------------------------------------------------------
    !        ... Form the washout rates
    !-----------------------------------------------------------------------      
    if ( gas_wetdep_method=='MOZ' ) then
       call sethet( het_rates, pmid, zmid, phis, tfld, &
                    cmfdqr, prain, nevapr, delt, invariants(:,:,indexm), &
                    vmr, ncol, lchnk )
       if (.not. convproc_do_aer) then
          call het_diags( het_rates(:ncol,:,:), mmr(:ncol,:,:), pdel(:ncol,:), lchnk, ncol )
       endif
    else
       het_rates = 0._r8
    end if
!
! CCMI
!
! set loss to below the tropopause only
!
    if ( st80_25_tau_ndx > 0 ) then
       do i = 1,ncol
          reaction_rates(i,1:troplev(i),st80_25_tau_ndx) = 0._r8
       enddo
    end if

    if ( has_linoz_data ) then
       ltrop_sol(:ncol) = troplev(:ncol)
    else
       ltrop_sol(:ncol) = 0 ! apply solver to all levels
    endif

    ! save h2so4 before gas phase chem (for later new particle nucleation)
    if (ndx_h2so4 > 0) then
       del_h2so4_gasprod(1:ncol,:) = vmr(1:ncol,:,ndx_h2so4)
    else
       del_h2so4_gasprod(:,:) = 0.0_r8
    endif

    vmr0(:ncol,:,:) = vmr(:ncol,:,:) ! mixing ratios before chemistry changes

    !=======================================================================
    !        ... Call the class solution algorithms
    !=======================================================================
    !-----------------------------------------------------------------------
    !	... Solve for "Explicit" species
    !-----------------------------------------------------------------------
    call exp_sol( vmr, reaction_rates, het_rates, extfrc, delt, invariants(1,1,indexm), ncol, lchnk, ltrop_sol )

    !-----------------------------------------------------------------------
    !	... Solve for "Implicit" species
    !-----------------------------------------------------------------------
    if ( has_strato_chem ) wrk(:,:) = vmr(:,:,h2o_ndx)
    call t_startf('imp_sol')
    !
    call imp_sol( vmr, reaction_rates, het_rates, extfrc, delt, &
                  ncol,pver, lchnk,  prod_out, loss_out )

    call t_stopf('imp_sol')

    call chem_prod_loss_diags_out( ncol, lchnk, vmr, reaction_rates, prod_out, loss_out, invariants(:ncol,:,indexm) )
    if( h2o_ndx>0) call outfld( 'H2O_GAS',  vmr(1,1,h2o_ndx),  ncol ,lchnk )

    ! reset O3S to O3 in the stratosphere ...
    if ( o3_ndx > 0 .and. o3s_ndx > 0 ) then
       o3s_loss = rate_diags_o3s_loss( reaction_rates, vmr, ncol )
       do i = 1,ncol
          vmr(i,1:troplev(i),o3s_ndx) = vmr(i,1:troplev(i),o3_ndx)
          vmr(i,troplev(i)+1:pver,o3s_ndx) = vmr(i,troplev(i)+1:pver,o3s_ndx) * exp(-delt*o3s_loss(i,troplev(i)+1:pver))
       end do
       call outfld( 'O3S_LOSS',  o3s_loss,  ncol ,lchnk )
    end if

    if (convproc_do_aer) then
       call vmr2mmr( vmr(:ncol,:,:), mmr_new(:ncol,:,:), mbar(:ncol,:), ncol )
       ! mmr_new = average of mmr values before and after imp_sol
       mmr_new(:ncol,:,:) = 0.5_r8*( mmr(:ncol,:,:) + mmr_new(:ncol,:,:) )
       call het_diags( het_rates(:ncol,:,:), mmr_new(:ncol,:,:), pdel(:ncol,:), lchnk, ncol )
    endif

    ! save h2so4 change by gas phase chem (for later new particle nucleation)
    if (ndx_h2so4 > 0) then
       del_h2so4_gasprod(1:ncol,:) = vmr(1:ncol,:,ndx_h2so4) - del_h2so4_gasprod(1:ncol,:)
    endif

!
! Aerosol processes ...
!

    call aero_model_gasaerexch( imozart-1, ncol, lchnk, troplevchem, delt, reaction_rates, &
                                tfld, pmid, pdel, mbar, relhum, &
                                zm,  qh2o, cwat, cldfr, ncldwtr, &
                                invariants(:,:,indexm), invariants, del_h2so4_gasprod,  &
                                vmr0, vmr, pbuf )

    if ( has_strato_chem ) then 

       wrk(:ncol,:) = (vmr(:ncol,:,h2o_ndx) - wrk(:ncol,:))*delt_inverse
       call outfld( 'QDCHEM',   wrk(:ncol,:),         ncol, lchnk )
       call outfld( 'HNO3_GAS', vmr(:ncol,:,hno3_ndx), ncol ,lchnk )

       !-----------------------------------------------------------------------      
       !         ... aerosol settling
       !             first settle hno3(2) using radius ice
       !             secnd settle hno3(3) using radius large nat
       !-----------------------------------------------------------------------      
       wrk(:,:) = vmr(:,:,h2o_ndx)
#ifdef ALT_SETTL
       where( h2o_cond(:,:) > 0._r8 )
          settl_rad(:,:) = radius_strat(:,:,3)
       elsewhere
          settl_rad(:,:) = 0._r8
       endwhere
       call strat_aer_settling( invariants(1,1,indexm), pmid, delt, zmid, tfld, &
            hno3_cond(1,1,2), settl_rad, ncol, lchnk, 1 )

       where( h2o_cond(:,:) == 0._r8 )
          settl_rad(:,:) = radius_strat(:,:,2)
       elsewhere
          settl_rad(:,:) = 0._r8
       endwhere
       call strat_aer_settling( invariants(1,1,indexm), pmid, delt, zmid, tfld, &
            hno3_cond(1,1,2), settl_rad, ncol, lchnk, 2 )
#else
       call strat_aer_settling( invariants(1,1,indexm), pmid, delt, zmid, tfld, &
            hno3_cond(1,1,2), radius_strat(1,1,2), ncol, lchnk, 2 )
#endif

       !-----------------------------------------------------------------------      
       !	... reform total hno3 and hcl = gas + all condensed
       !-----------------------------------------------------------------------      
!      NOTE: vmr for hcl and hno3 is gas-phase at this point.
!            hno3_cond(:,k,1) = STS; hno3_cond(:,k,2) = NAT
   
       do k = 1,pver
          vmr(:,k,hno3_ndx) = vmr(:,k,hno3_ndx) + hno3_cond(:,k,1) &
               + hno3_cond(:,k,2) 
          vmr(:,k,hcl_ndx)  = vmr(:,k,hcl_ndx)  + hcl_cond(:,k) 
              
       end do

       wrk(:,:) = (vmr(:,:,h2o_ndx) - wrk(:,:))*delt_inverse
       call outfld( 'QDSETT', wrk(:,:), ncol, lchnk )

    endif

!
! LINOZ
!
    if ( do_lin_strat_chem ) then
       call lin_strat_chem_solve( ncol, lchnk, vmr(:,:,o3_ndx), col_dens(:,:,1), tfld, zen_angle, pmid, delt, rlats, troplev )
    end if

    !-----------------------------------------------------------------------      
    !         ... Check for negative values and reset to zero
    !-----------------------------------------------------------------------      
    call negtrc( 'After chemistry ', vmr, ncol )

    !-----------------------------------------------------------------------      
    !         ... Set upper boundary mmr values
    !-----------------------------------------------------------------------      
    call set_fstrat_vals( vmr, pmid, pint, troplev, calday, ncol,lchnk )

    !-----------------------------------------------------------------------      
    !         ... Set fixed lower boundary mmr values
    !-----------------------------------------------------------------------      
    call flbc_set( vmr, ncol, lchnk, map2chm )

    !----------------------------------------------------------------------- 
    ! set NOy UBC     
    !-----------------------------------------------------------------------      
    call noy_ubc_set( lchnk, ncol, vmr )

    if ( ghg_chem ) then
       call ghg_chem_set_flbc( vmr, ncol )
    endif

    !-----------------------------------------------------------------------
    ! force ion/electron balance -- ext forcings likely do not conserve charge
    !-----------------------------------------------------------------------      
    call charge_balance( ncol, vmr )

    !-----------------------------------------------------------------------      
    !         ... Xform from vmr to mmr
    !-----------------------------------------------------------------------      
    call vmr2mmr( vmr(:ncol,:,:), mmr_tend(:ncol,:,:), mbar(:ncol,:), ncol )

    call set_short_lived_species( mmr_tend, lchnk, ncol, pbuf )

    !-----------------------------------------------------------------------      
    !         ... Form the tendencies
    !----------------------------------------------------------------------- 
    do m = 1,gas_pcnst 
       mmr_new(:ncol,:,m) = mmr_tend(:ncol,:,m)
       mmr_tend(:ncol,:,m) = (mmr_tend(:ncol,:,m) - mmr(:ncol,:,m))*delt_inverse
    enddo

    do m = 1,pcnst
       n = map2chm(m)
       if( n > 0 ) then
          qtend(:ncol,:,m) = qtend(:ncol,:,m) + mmr_tend(:ncol,:,n) 
       end if
    end do

    tvs(:ncol) = tfld(:ncol,pver) * (1._r8 + qh2o(:ncol,pver))

    sflx(:,:) = 0._r8
    call get_ref_date(yr, mon, day, sec)
    ncdate = yr*10000 + mon*100 + day
    wind_speed(:ncol) = sqrt( ufld(:ncol,pver)*ufld(:ncol,pver) + vfld(:ncol,pver)*vfld(:ncol,pver) )
    prect(:ncol) = precc(:ncol) + precl(:ncol)

    if ( drydep_method == DD_XLND ) then
       soilw = -99
       call drydep( ocnfrac, icefrac, ncdate, ts, ps,  &
            wind_speed, qh2o(:,pver), tfld(:,pver), pmid(:,pver), prect, &
            snowhland, fsds, depvel, sflx, mmr, &
            tvs, soilw, relhum(:,pver:pver), ncol, lonndx, latndx, lchnk )
    else if ( drydep_method == DD_XATM ) then
       table_soilw = has_drydep( 'H2' ) .or. has_drydep( 'CO' )
       if( .not. dyn_soilw .and. table_soilw ) then
          call set_soilw( soilw, lchnk, calday )
       end if
       call drydep( ncdate, ts, ps,  &
            wind_speed, qh2o(:,pver), tfld(:,pver), pmid(:,pver), prect, &
            snowhland, fsds, depvel, sflx, mmr, &
            tvs, soilw, relhum(:,pver:pver), ncol, lonndx, latndx, lchnk )
    else if ( drydep_method == DD_TABL ) then
       call drydep( calday, ts, zen_angle, &
            depvel, sflx, mmr, pmid(:,pver), &
            tvs, ncol, icefrac, ocnfrac, lchnk )
    endif

    drydepflx(:,:) = 0._r8
    wetdepflx_diag(:,:) = 0._r8
    do m = 1,pcnst
       n = map2chm( m )
       if ( n > 0 ) then
         cflx(:ncol,m)      = cflx(:ncol,m) - sflx(:ncol,n)
         drydepflx(:ncol,m) = sflx(:ncol,n)
         wetdepflx_diag(:ncol,n) = wetdepflx(:ncol,m)
       endif
    end do

    call chm_diags( lchnk, ncol, vmr(:ncol,:,:), mmr_new(:ncol,:,:), &
                    reaction_rates(:ncol,:,:), invariants(:ncol,:,:), depvel(:ncol,:),  sflx(:ncol,:), &
                    mmr_tend(:ncol,:,:), pdel(:ncol,:), pmid(:ncol,:), troplev(:ncol), wetdepflx_diag(:ncol,:), &
                    nhx_nitrogen_flx(:ncol), noy_nitrogen_flx(:ncol) )

    call rate_diags_calc( reaction_rates(:,:,:), vmr(:,:,:), invariants(:,:,indexm), ncol, lchnk )
!
! jfl
!
! surface vmr
!
    if ( pm25_srf_diag ) then
       pm25(:ncol) = mmr_new(:ncol,pver,cb1_ndx)   &
            + mmr_new(:ncol,pver,cb2_ndx)   &
            + mmr_new(:ncol,pver,oc1_ndx)   &
            + mmr_new(:ncol,pver,oc2_ndx)   &
            + mmr_new(:ncol,pver,dst1_ndx)  &
            + mmr_new(:ncol,pver,dst2_ndx)  &
            + mmr_new(:ncol,pver,sslt1_ndx) &
            + mmr_new(:ncol,pver,sslt2_ndx) &
            + mmr_new(:ncol,pver,soa_ndx)   &
            + mmr_new(:ncol,pver,so4_ndx)
       call outfld('PM25_SRF',pm25(:ncol) , ncol, lchnk )
    endif
    if ( pm25_srf_diag_soa ) then
       pm25(:ncol) = mmr_new(:ncol,pver,cb1_ndx)   &
            + mmr_new(:ncol,pver,cb2_ndx)   &
            + mmr_new(:ncol,pver,oc1_ndx)   &
            + mmr_new(:ncol,pver,oc2_ndx)   &
            + mmr_new(:ncol,pver,dst1_ndx)  &
            + mmr_new(:ncol,pver,dst2_ndx)  &
            + mmr_new(:ncol,pver,sslt1_ndx) &
            + mmr_new(:ncol,pver,sslt2_ndx) &
            + mmr_new(:ncol,pver,soam_ndx)   &
            + mmr_new(:ncol,pver,soai_ndx)   &
            + mmr_new(:ncol,pver,soat_ndx)   &
            + mmr_new(:ncol,pver,soab_ndx)   &
            + mmr_new(:ncol,pver,soax_ndx)   &
            + mmr_new(:ncol,pver,so4_ndx)
       call outfld('PM25_SRF',pm25(:ncol) , ncol, lchnk )
    endif
!
!
    call outfld('Q_SRF',qh2o(:ncol,pver) , ncol, lchnk )
    call outfld('U_SRF',ufld(:ncol,pver) , ncol, lchnk )
    call outfld('V_SRF',vfld(:ncol,pver) , ncol, lchnk )

!
    if (.not.sad_pbf_ndx>0) then
       deallocate(strato_sad)
    endif

  end subroutine gas_phase_chemdr

end module mo_gas_phase_chemdr
