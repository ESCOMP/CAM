module module_mosaic_init_aerpar
  
  implicit none
  private

  public:: mosaic_init_aer_params

contains
  subroutine mosaic_init_aer_params
    !BSINGH - All initialzations for Mosiac model

    call load_mosaic_parameters
    
  end subroutine mosaic_init_aer_params

  !---------------------------------------------------------------------------------------!
  !BSINGH: load_mosaic_parameters subroutine is directly copied form the mosaic_box.25.f90 
  !        code
  !---------------------------------------------------------------------------------------!

  !---------------------------------------------------------------------------------------!
  ! Called only once per entire simulation to load gas and aerosol
  ! indices, parameters, physico-chemical constants, polynomial coeffs, etc.
  !
  ! author: Rahul A. Zaveri
  ! update: jan 2005
  !---------------------------------------------------------------------------------------!
  subroutine load_mosaic_parameters

    !      include 'v33com2'
    use module_data_mosaic_aero, only: ipmcmos_aero, no_aerosol, all_solid, all_liquid,   &
         mixed, nelectrolyte, naercomp, naer, Ncation, Nanion,                            &
         ngas_aerchtot, ngas_volatile, nsalt,                                             &
         jsulf_poor_NUM, jsulf_rich_NUM, MDRH_T_NUM, d_mdrh_DIM2, phasestate, aer_name,   &
         gas_name, ename, jnh4so4, jlvcite, jnh4hso4, jnh4msa, jnh4no3, jnh4cl, jna2so4,  &
         jna3hso4, jnahso4, jnamsa, jnano3, jnacl, jcano3, jcacl2, jcamsa2, jh2so4, jmsa, &
         jhno3, jhcl, jhhso4, jcaso4, jcaco3, joc, jbc, join, jaro1, jaro2, jalk1, jole1, &
         japi1, japi2, jlim1, jlim2, jh2o, jc_h, jc_nh4, jc_na, jc_ca, ja_hso4, ja_so4,   &
         ja_no3, ja_cl, ja_msa, ih2so4_g, ihno3_g, ihcl_g, inh3_g, imsa_g, iaro1_g,       &
         iaro2_g, ialk1_g, iole1_g, iapi1_g, iapi2_g, ilim1_g, ilim2_g, iso4_a, ino3_a,   &
         icl_a, inh4_a, imsa_a, iaro1_a, iaro2_a, ialk1_a, iole1_a, iapi1_a, iapi2_a,     &
         ilim1_a, ilim2_a, ico3_a, ina_a, ica_a, ioin_a, ioc_a, ibc_a,                    &
         isoa_first, jsoa_first, nmax_ASTEM, b_mtem,                                      &
         zc, za, b_zsr, a_zsr, aw_min, mw_electrolyte, dens_electrolyte,                  &
         partial_molar_vol, MW_c, MW_a, mw_aer_mac,dens_aer_mac, kappa_aer_mac,           &
         dens_comp_a, mw_comp_a, ref_index_a, mw_gas, v_molar_gas,                        &
         rtol_mesa, jsalt_index, jsulf_poor, jsulf_rich, Nmax_mesa, d_mdrh,  &
         use_cam5mam_soa_params

    use module_data_mosaic_kind, only: r8

    implicit none

    ! local variables
    integer iaer, igas, je, ja, j_index, ibin
    logical use_mos31e_rz1_densities, use_uniform_densities !BSINGH - 05/28/2013(RCE updates)
    real(r8), dimension(nelectrolyte) :: G_MX,K_MX

    !BSINGH - 05/28/2013(RCE updates)
    use_mos31e_rz1_densities = .true.
    if ( use_mos31e_rz1_densities ) then
       use_uniform_densities = .false.
    else
       use_uniform_densities = .true.
       if (ipmcmos_aero > 0) use_uniform_densities = .false.
    end if
    !BSINGH - 05/28/2013(RCE updates ENDS)

    ! rce 2013-07-31 -
    !   using a local saved variable like "first" no longer works
    !   the calling routine needs to determine if/when this routine is needed
    !   if(first)then
    !      first=.false.

    !----------------------------------------------------------------
    ! control settings
    ! *** do not change mSIZE_FRAMEWORK here ***
    !      mSIZE_FRAMEWORK = mSECTIONAL        ! mMODAL or mSECTIONAL
    !      mDYNAMIC_SOLVER = mASTEM            ! mASTEM, mLSODES
    !      mGAS_AER_XFER   = mON               ! mON, mOFF

    ! ASTEM parameters
    nmax_ASTEM      = 301              ! max number of time steps in ASTEM
    !      alpha_ASTEM     = 1.0               ! choose a value between 0.01 and 1.0
    !      rtol_eqb_ASTEM = 0.01               ! equilibrium tolerance in ASTEM
    !      ptol_mol_ASTEM = 0.01               ! mol percent tolerance in ASTEM

    ! MESA parameters
    Nmax_MESA       = 80               ! max number of iterations in MESA_PTC
    rtol_mesa       = 0.01_r8             ! MESA equilibrium tolerance
    !----------------------------------------------------------------

    !
    ! set species indices
    !    ixxx_a are for aerosol     species - they apply to aer and total arrays
    !    ixxx_g are for gas         species - they apply to gas array
    !    jyyy   are for electrolyte species - they apply to electrolyte and comp_a arrays
    !
    ! *** CRITICAL "RULES" ABOUT THESE INDICES ***
    !
    ! *** THE CODE WILL FAIL IF THEY ARE NOT FOLLOWED  ***
    !
    ! for the volatile inorganic and organic species (1 thru ngas_volatile)
    !    the ixxx_a and ixxx_g must match (be indentical) for each aerosol-gas pair
    !
    ! for the soa species, the ordering of the ixxx_a, ixxx_g, and jyyy
    !    must be identical
    !
    
    ! electrolyte indices (used for water content calculations)
    ! these indices are order sensitive
    ! inorganic species first
    jnh4so4    =  1    ! soluble
    jlvcite    =  2    ! soluble
    jnh4hso4   =  3    ! soluble
    jnh4msa    =  4    ! soluble: new
    jnh4no3    =  5    ! soluble
    jnh4cl     =  6    ! soluble
    jna2so4    =  7    ! soluble
    jna3hso4   =  8    ! soluble
    jnahso4    =  9    ! soluble
    jnamsa     = 10    ! soluble: new
    jnano3     = 11    ! soluble
    jnacl      = 12    ! soluble
    jcano3     = 13    ! soluble
    jcacl2     = 14    ! soluble
    jcamsa2    = 15    ! soluble     nsalt
    jh2so4     = 16    ! soluble
    jmsa       = 17    ! soluble
    jhno3      = 18    ! soluble
    jhcl       = 19    ! soluble
    jhhso4     = 20    ! soluble
    jcaso4     = 21    ! insoluble
    jcaco3     = 22    ! insoluble
    joc        = 23    ! insoluble - part of naercomp
    jbc        = 24    ! insoluble - part of naercomp
    join       = 25    ! insoluble - part of naercomp

    ! aerosol and gas indices for inorganic species
    iso4_a     =  1 ;  ih2so4_g   =  1     
    ino3_a     =  2 ;  ihno3_g    =  2     
    icl_a      =  3 ;  ihcl_g     =  3     
    inh4_a     =  4 ;  inh3_g     =  4     
    imsa_a     =  5 ;  imsa_g     =  5     

    ! aerosol, gas, and electrolyte indices for secondary organic species
    iaro1_a    =  6 ;  iaro1_g    =  6 ;  jaro1      = 26     
    iaro2_a    =  7 ;  iaro2_g    =  7 ;  jaro2      = 27     
    ialk1_a    =  8 ;  ialk1_g    =  8 ;  jalk1      = 28     
    iole1_a    =  9 ;  iole1_g    =  9 ;  jole1      = 29     
    iapi1_a    = 10 ;  iapi1_g    = 10 ;  japi1      = 30     
    iapi2_a    = 11 ;  iapi2_g    = 11 ;  japi2      = 31     
    ilim1_a    = 12 ;  ilim1_g    = 12 ;  jlim1      = 32     
    ilim2_a    = 13 ;  ilim2_g    = 13 ;  jlim2      = 33     

    isoa_first = iaro1_g
    jsoa_first = jaro1

    ! electrolyte indices for other species
    jh2o       = 34    ! water - part of naercomp

    ! aerosol indices for other species
    ico3_a     = 14
    ina_a      = 15
    ica_a      = 16
    ioin_a     = 17
    ioc_a      = 18
    ibc_a      = 19

    ! gas indices for other species
    ! ico2_g   = 14  ! *** currently not used
                     ! *** if co3_a was treated as a volatile inorganic, 
                     !     then ico3_a and ico2_g would have to be 6 

    ! local aerosol ions
    ! cations
    jc_h       =  1
    jc_nh4     =  2
    jc_na      =  3
    jc_ca      =  4
    !
    ! anions
    ja_hso4    =  1
    ja_so4     =  2
    ja_no3     =  3
    ja_cl      =  4
    ja_msa     =  5
    ! ja_co3   =  6   ! *** currently not used


    !--------------------------------------------------------------------
    ! phase state names
    phasestate(no_aerosol) = "NOAERO"
    phasestate(all_solid)  = "SOLID "
    phasestate(all_liquid) = "LIQUID"
    phasestate(mixed)      = "MIXED "

    ! names of aer species
    do iaer = 1, naer
       write( aer_name(iaer), '(a,i4.4)' ) 'aer', iaer  ! default
    end do

    aer_name(iso4_a) = "SO4"
    aer_name(ino3_a) = "NO3"
    aer_name(icl_a)  = "Cl "
    aer_name(inh4_a) = "NH4"
    aer_name(ioc_a)  = "OC "
    aer_name(imsa_a) = "MSA"
    aer_name(ico3_a) = "CO3"
    aer_name(ina_a)  = "Na "
    aer_name(ica_a)  = "Ca "
    aer_name(ibc_a)  = "BC "
    aer_name(ioin_a) = "OIN"
    aer_name(iaro1_a)= "ARO1"
    aer_name(iaro2_a)= "ARO2"
    aer_name(ialk1_a)= "ALK1"
    aer_name(iole1_a)= "OLE1"
    aer_name(iapi1_a)= "API1"
    aer_name(iapi2_a)= "API2"
    aer_name(ilim1_a)= "LIM1"
    aer_name(ilim2_a)= "LIM2"

    ! names of gas species
    do igas = 1, ngas_aerchtot
       write( gas_name(igas), '(a,i4.4)' ) 'gas', igas  ! default
    end do

    gas_name(ih2so4_g) = "H2SO4"
    gas_name(ihno3_g)  = "HNO3 "
    gas_name(ihcl_g)   = "HCl  "
    gas_name(inh3_g)   = "NH3  "
    gas_name(imsa_g)   = "MSA  "
    gas_name(iaro1_g)   = "ARO1 "
    gas_name(iaro2_g)   = "ARO2 "
    gas_name(ialk1_g)   = "ALK1 "
    gas_name(iole1_g)   = "OLE1 "
    gas_name(iapi1_g)   = "API1 "
    gas_name(iapi2_g)   = "API2 "
    gas_name(ilim1_g)   = "LIM1 "
    gas_name(ilim2_g)   = "LIM2 "

    ! names of electrolytes
    ename(jnh4so4) = "AmSO4"
    ename(jlvcite) = "(NH4)3H(SO4)2"
    ename(jnh4hso4)= "NH4HSO4"
    ename(jnh4msa) = "CH3SO3NH4"
    ename(jnh4no3) = "NH4NO3"
    ename(jnh4cl)  = "NH4Cl"
    ename(jnacl)   = "NaCl"
    ename(jnano3)  = "NaNO3"
    ename(jna2so4) = "Na2SO4"
    ename(jna3hso4)= "Na3H(SO4)2"
    ename(jnamsa)  = "CH3SO3Na"
    ename(jnahso4) = "NaHSO4"
    ename(jcaso4)  = "CaSO4"
    ename(jcamsa2) = "(CH3SO3)2Ca"
    ename(jcano3)  = "Ca(NO3)2"
    ename(jcacl2)  = "CaCl2"
    ename(jcaco3)  = "CaCO3"
    ename(jh2so4)  = "H2SO4"
    ename(jhhso4)  = "HHSO4"
    ename(jhno3)   = "HNO3"
    ename(jhcl)    = "HCl"
    ename(jmsa)    = "CH3SO3H"

    ! molecular weights of electrolytes
    mw_electrolyte(jnh4so4) = 132.0_r8
    mw_electrolyte(jlvcite) = 247.0_r8
    mw_electrolyte(jnh4hso4)= 115.0_r8
    mw_electrolyte(jnh4msa) = 113.0_r8
    mw_electrolyte(jnh4no3) = 80.0_r8
    mw_electrolyte(jnh4cl)  = 53.5_r8
    mw_electrolyte(jnacl)   = 58.5_r8
    mw_electrolyte(jnano3)  = 85.0_r8
    mw_electrolyte(jna2so4) = 142.0_r8
    mw_electrolyte(jna3hso4)= 262.0_r8
    mw_electrolyte(jnahso4) = 120.0_r8
    mw_electrolyte(jnamsa)  = 118.0_r8
    mw_electrolyte(jcaso4)  = 136.0_r8
    mw_electrolyte(jcamsa2) = 230.0_r8
    mw_electrolyte(jcano3)  = 164.0_r8
    mw_electrolyte(jcacl2)  = 111.0_r8
    mw_electrolyte(jcaco3)  = 100.0_r8
    mw_electrolyte(jh2so4)  = 98.0_r8
    mw_electrolyte(jhno3)   = 63.0_r8
    mw_electrolyte(jhcl)    = 36.5_r8
    mw_electrolyte(jmsa)    = 96.0_r8


    ! molecular weights of ions [g/mol]
    MW_c(jc_h)  =  1.0_r8
    MW_c(jc_nh4)= 18.0_r8
    MW_c(jc_na) = 23.0_r8
    MW_c(jc_ca) = 40.0_r8

    MW_a(ja_so4) = 96.0_r8
    MW_a(ja_hso4)= 97.0_r8
    MW_a(ja_no3) = 62.0_r8
    MW_a(ja_cl)  = 35.5_r8
    MW_a(ja_msa) = 95.0_r8


    ! magnitude of the charges on ions
    zc(jc_h)   = 1
    zc(jc_nh4) = 1
    zc(jc_na)  = 1
    zc(jc_ca)  = 2

    za(ja_hso4)= 1
    za(ja_so4) = 2
    za(ja_no3) = 1
    za(ja_cl)  = 1
    za(ja_msa) = 1


    ! densities of pure electrolytes in g/cc
    dens_electrolyte(jnh4so4)  = 1.8_r8
    dens_electrolyte(jlvcite)  = 1.8_r8
    dens_electrolyte(jnh4hso4) = 1.8_r8
    dens_electrolyte(jnh4msa)  = 1.8_r8 ! assumed same as nh4hso4
    dens_electrolyte(jnh4no3)  = 1.8_r8
    dens_electrolyte(jnh4cl)   = 1.8_r8
    dens_electrolyte(jnacl)    = 2.2_r8
    dens_electrolyte(jnano3)   = 2.2_r8
    dens_electrolyte(jna2so4)  = 2.2_r8
    dens_electrolyte(jna3hso4) = 2.2_r8
    dens_electrolyte(jnahso4)  = 2.2_r8
    dens_electrolyte(jnamsa)   = 2.2_r8 ! assumed same as nahso4
    dens_electrolyte(jcaso4)   = 2.6_r8
    dens_electrolyte(jcamsa2)  = 2.6_r8   ! assumed same as caso4
    dens_electrolyte(jcano3)   = 2.6_r8
    dens_electrolyte(jcacl2)   = 2.6_r8
    dens_electrolyte(jcaco3)   = 2.6_r8
    dens_electrolyte(jh2so4)   = 1.8_r8
    dens_electrolyte(jhhso4)   = 1.8_r8
    dens_electrolyte(jhno3)    = 1.8_r8
    dens_electrolyte(jhcl)     = 1.8_r8
    dens_electrolyte(jmsa)     = 1.8_r8 ! assumed same as h2so4
    if ( use_uniform_densities ) then!BSINGH - 05/28/2013(RCE updates)
       do je = 1, nelectrolyte
          dens_electrolyte(je) = 1.6_r8
       enddo
    endif!BSINGH - 05/28/2013(RCE updates)

    ! densities of compounds in g/cc
    dens_comp_a(jnh4so4)  = 1.8_r8
    dens_comp_a(jlvcite)  = 1.8_r8
    dens_comp_a(jnh4hso4) = 1.8_r8
    dens_comp_a(jnh4msa)  = 1.8_r8        ! assumed same as nh4hso4
    dens_comp_a(jnh4no3)  = 1.7_r8
    dens_comp_a(jnh4cl)   = 1.5_r8
    dens_comp_a(jnacl)    = 2.2_r8
    dens_comp_a(jnano3)   = 2.2_r8
    dens_comp_a(jna2so4)  = 2.2_r8
    dens_comp_a(jna3hso4) = 2.2_r8
    dens_comp_a(jnahso4)  = 2.2_r8
    dens_comp_a(jnamsa)   = 2.2_r8        ! assumed same as nahso4
    dens_comp_a(jcaso4)   = 2.6_r8
    dens_comp_a(jcamsa2)  = 2.6_r8        ! assumed same as caso4
    dens_comp_a(jcano3)   = 2.6_r8
    dens_comp_a(jcacl2)   = 2.6_r8
    dens_comp_a(jcaco3)   = 2.6_r8
    dens_comp_a(jh2so4)   = 1.8_r8
    dens_comp_a(jhhso4)   = 1.8_r8
    dens_comp_a(jhno3)    = 1.8_r8
    dens_comp_a(jhcl)     = 1.8_r8
    dens_comp_a(jmsa)     = 1.8_r8        ! assumed same as h2so4
    dens_comp_a(joc)      = 1.0_r8
    dens_comp_a(jbc)      = 1.8_r8
    dens_comp_a(join)     = 2.6_r8
    dens_comp_a(jaro1)    = 1.0_r8
    dens_comp_a(jaro2)    = 1.0_r8
    dens_comp_a(jalk1)    = 1.0_r8
    dens_comp_a(jole1)    = 1.0_r8
    dens_comp_a(japi1)    = 1.0_r8
    dens_comp_a(japi2)    = 1.0_r8
    dens_comp_a(jlim1)    = 1.0_r8
    dens_comp_a(jlim2)    = 1.0_r8
    dens_comp_a(jh2o)     = 1.0_r8
    !BSINGH - 05/28/2013(RCE updates)
    ! following for comparison with mos31d_bs2 and. mos31e_rz1
    if ( use_mos31e_rz1_densities ) then
       dens_comp_a(joc)      = 1.4_r8
       dens_comp_a(jaro1)    = 1.4_r8
       dens_comp_a(jaro2)    = 1.4_r8
       dens_comp_a(jalk1)    = 1.4_r8
       dens_comp_a(jole1)    = 1.4_r8
       dens_comp_a(japi1)    = 1.4_r8
       dens_comp_a(japi2)    = 1.4_r8
       dens_comp_a(jlim1)    = 1.4_r8
       dens_comp_a(jlim2)    = 1.4_r8
    end if

    if ( use_uniform_densities ) then
       !BSINGH - 05/28/2013(RCE updates ENDS)
       do je = 1, naercomp
          dens_comp_a(je) = 1.6_r8
       enddo
       !BSINGH - 05/28/2013(RCE updates)
    endif

    if (ipmcmos_aero > 0) then
       dens_comp_a(jnh4no3)  = 1.8_r8
       dens_comp_a(jnh4cl)   = 1.8_r8
       dens_comp_a(jaro1)    = 1.4_r8
       dens_comp_a(jaro2)    = 1.4_r8
       dens_comp_a(jalk1)    = 1.4_r8
       dens_comp_a(jole1)    = 1.4_r8
       dens_comp_a(japi1)    = 1.4_r8
       dens_comp_a(japi2)    = 1.4_r8
       dens_comp_a(jlim1)    = 1.4_r8
       dens_comp_a(jlim2)    = 1.4_r8
    endif
    !BSINGH - 05/28/2013(RCE updates ENDS)

    ! molecular weights of generic aerosol species
    mw_aer_mac(1:naer) = 200.0_r8  ! default

    mw_aer_mac(iso4_a) = 96.0_r8
    mw_aer_mac(ino3_a) = 62.0_r8
    mw_aer_mac(icl_a)  = 35.5_r8
    mw_aer_mac(imsa_a) = 95.0_r8  ! CH3SO3
    mw_aer_mac(ico3_a) = 60.0_r8
    mw_aer_mac(inh4_a) = 18.0_r8
    mw_aer_mac(ina_a)  = 23.0_r8
    mw_aer_mac(ica_a)  = 40.0_r8
    mw_aer_mac(ioin_a) = 1.0_r8           ! not used
    mw_aer_mac(ibc_a)  = 1.0_r8           ! not used
    mw_aer_mac(ioc_a)  = 1.0_r8   ! 200 assumed for primary organics
    mw_aer_mac(iaro1_a)= 250.0_r8
    mw_aer_mac(iaro2_a)= 250.0_r8
    mw_aer_mac(ialk1_a)= 250.0_r8
    mw_aer_mac(iole1_a)= 250.0_r8
    mw_aer_mac(iapi1_a)= 250.0_r8
    mw_aer_mac(iapi2_a)= 250.0_r8
    mw_aer_mac(ilim1_a)= 250.0_r8
    mw_aer_mac(ilim2_a)= 250.0_r8

    ! molecular weights of compounds
    mw_comp_a(jnh4so4) = 132.0_r8
    mw_comp_a(jlvcite) = 247.0_r8
    mw_comp_a(jnh4hso4)= 115.0_r8
    mw_comp_a(jnh4msa) = 113.0_r8
    mw_comp_a(jnh4no3) = 80.0_r8
    mw_comp_a(jnh4cl)  = 53.5_r8
    mw_comp_a(jnacl)   = 58.5_r8
    mw_comp_a(jnano3)  = 85.0_r8
    mw_comp_a(jna2so4) = 142.0_r8
    mw_comp_a(jna3hso4)= 262.0_r8
    mw_comp_a(jnahso4) = 120.0_r8
    mw_comp_a(jnamsa)  = 118.0_r8
    mw_comp_a(jcaso4)  = 136.0_r8
    mw_comp_a(jcamsa2) = 230.0_r8
    mw_comp_a(jcano3)  = 164.0_r8
    mw_comp_a(jcacl2)  = 111.0_r8
    mw_comp_a(jcaco3)  = 100.0_r8
    mw_comp_a(jh2so4)  = 98.0_r8
    mw_comp_a(jhhso4)  = 98.0_r8
    mw_comp_a(jhno3)   = 63.0_r8
    mw_comp_a(jhcl)    = 36.5_r8
    mw_comp_a(jmsa)    = 96.0_r8
    mw_comp_a(joc)      = 1.0_r8
    mw_comp_a(jbc)      = 1.0_r8
    mw_comp_a(join)    = 1.0_r8
    mw_comp_a(jaro1)    = 250.0_r8
    mw_comp_a(jaro2)    = 250.0_r8
    mw_comp_a(jalk1)    = 250.0_r8
    mw_comp_a(jole1)    = 250.0_r8
    mw_comp_a(japi1)    = 250.0_r8
    mw_comp_a(japi2)    = 250.0_r8
    mw_comp_a(jlim1)    = 250.0_r8
    mw_comp_a(jlim2)    = 250.0_r8
    mw_comp_a(jh2o)    = 18.0_r8
    !BSINGH - 05/28/2013(RCE updates)
    ! partmc-2.2.1 jun-2012
    !#     dens (kg/m^3)   ions in soln (1)    molec wght (kg/mole)   kappa (1)
    ! SO4            1800                0                   96d-3      0.65
    ! NO3            1800                0                   62d-3      0.65
    ! Cl             2200                0                   35.5d-3    0.53
    ! NH4            1800                0                   18d-3      0.65
    ! MSA            1800                0                   95d-3      0.53
    ! ARO1           1400                0                   150d-3     0.1
    ! ARO2           1400                0                   150d-3     0.1
    ! ALK1           1400                0                   140d-3     0.1
    ! OLE1           1400                0                   140d-3     0.1
    ! API1           1400                0                   184d-3     0.1
    ! API2           1400                0                   184d-3     0.1
    ! LIM1           1400                0                   200d-3     0.1
    ! LIM2           1400                0                   200d-3     0.1
    ! CO3            2600                0                   60d-3      0.53
    ! Na             2200                0                   23d-3      0.53
    ! Ca             2600                0                   40d-3      0.53
    ! OIN            2600                0                   1d-3       0.1
    ! OC             1000                0                   1d-3       0.001
    ! BC             1800                0                   1d-3       0
    ! H2O            1000                0                   18d-3      0
    !BSINGH - 05/28/2013(RCE updates ENDS)

    ! densities of generic aerosol species
    dens_aer_mac(1:naer) = 1.0_r8   ! default

    dens_aer_mac(iso4_a) = 1.8_r8 ! used
    dens_aer_mac(ino3_a) = 1.8_r8 ! used
    dens_aer_mac(icl_a)  = 2.2_r8 ! used
    dens_aer_mac(imsa_a) = 1.8_r8 ! used
    dens_aer_mac(ico3_a) = 2.6_r8 ! used
    dens_aer_mac(inh4_a) = 1.8_r8 ! used
    dens_aer_mac(ina_a)  = 2.2_r8 ! used
    dens_aer_mac(ica_a)  = 2.6_r8 ! used
    dens_aer_mac(ioin_a) = 2.6_r8 ! used
    dens_aer_mac(ioc_a)  = 1.0_r8 ! used
    dens_aer_mac(ibc_a)  = 1.8_r8 ! used
    dens_aer_mac(iaro1_a)= 1.0_r8
    dens_aer_mac(iaro2_a)= 1.0_r8
    dens_aer_mac(ialk1_a)= 1.0_r8
    dens_aer_mac(iole1_a)= 1.0_r8
    dens_aer_mac(iapi1_a)= 1.0_r8
    dens_aer_mac(iapi2_a)= 1.0_r8
    dens_aer_mac(ilim1_a)= 1.0_r8
    dens_aer_mac(ilim2_a)= 1.0_r8
    !BSINGH - 05/28/2013(RCE updates)
    ! following for comparison with mos31d_bs2 and. mos31e_rz1
    if ( use_mos31e_rz1_densities ) then
       dens_aer_mac(ioc_a)  = 1.4_r8
       dens_aer_mac(iaro1_a)= 1.4_r8
       dens_aer_mac(iaro2_a)= 1.4_r8
       dens_aer_mac(ialk1_a)= 1.4_r8
       dens_aer_mac(iole1_a)= 1.4_r8
       dens_aer_mac(iapi1_a)= 1.4_r8
       dens_aer_mac(iapi2_a)= 1.4_r8
       dens_aer_mac(ilim1_a)= 1.4_r8
       dens_aer_mac(ilim2_a)= 1.4_r8
    end if

    if ( use_uniform_densities ) then
       !BSINGH - 05/28/2013(RCE updates ENDS)

       do iaer = 1, naer
          dens_aer_mac(iaer) = 1.6_r8
       enddo
    endif!BSINGH - 05/28/2013(RCE updates)

    if (ipmcmos_aero > 0) then
       ! use partmc-mosaic densities
       dens_aer_mac(1:19) = (/ &
            1.80_r8, 1.80_r8, 2.20_r8, 1.80_r8, 1.80_r8, 1.40_r8, 1.40_r8, 1.40_r8, 1.40_r8, 1.40_r8, &
            1.40_r8, 1.40_r8, 1.40_r8, 2.60_r8, 2.20_r8, 2.60_r8, 2.60_r8, 1.00_r8, 1.80_r8 /)!BSINGH - 05/28/2013(RCE updates)
       !           so4   no3   cl    nh4   msa   aro1  aro2  alk1  ole1  api1
       !           api2  lim1  lim2  co3   na    ca    oin   oc    bc
    end if

    if ( use_cam5mam_soa_params > 0 ) then
       dens_aer_mac(ioc_a)   = 1.0_r8
       dens_aer_mac(ilim2_a) = 1.0_r8
       ! for oc, leave mw=1 because some of the mosaic code requires this
       mw_aer_mac(ilim2_a)   = 150.0_r8
       dens_comp_a(joc)      = 1.0_r8
       dens_comp_a(jlim2)    = 1.0_r8
       mw_comp_a(jlim2)      = 250.0_r8
  
!zlu soa
dens_aer_mac(ilim1_a) = 1.0_r8
mw_aer_mac(ilim1_a)   = 250.0_r8
dens_comp_a(jlim1)    = 1.0_r8
mw_comp_a(jlim1)      = 250.0_r8


dens_aer_mac(iaro1_a) = 1.0_r8
mw_aer_mac(iaro1_a)   = 250.0_r8
dens_comp_a(jaro1)    = 1.0_r8
mw_comp_a(jaro1)      = 250.0_r8


dens_aer_mac(iapi1_a) = 1.0_r8
mw_aer_mac(iapi1_a)   = 250.0_r8
dens_comp_a(japi1)    = 1.0_r8
mw_comp_a(japi1)      = 250.0_r8

dens_aer_mac(iapi2_a) = 1.0_r8
mw_aer_mac(iapi2_a)   = 250.0_r8
dens_comp_a(japi2)    = 1.0_r8
mw_comp_a(japi2)      = 250.0_r8



     end if

    ! kappa values (hygroscopicities) of generic aerosol species
    !
    ! for calculation of ccn properties, kappa of electrolytes
    !    should be used
    ! the multi-dimensional sectional code needs a "fixed" kappa
    !    for each generic aerosol species, just as the older
    !    1d sectional code needs a "fixed" dry density
    kappa_aer_mac(1:naer)  = 0.1_r8  ! default

    kappa_aer_mac(iso4_a)  = 0.65_r8
    kappa_aer_mac(ino3_a)  = 0.65_r8
    kappa_aer_mac(imsa_a)  = 0.65_r8
    kappa_aer_mac(inh4_a)  = 0.65_r8
    kappa_aer_mac(icl_a)   = 0.65_r8
    kappa_aer_mac(ina_a)   = 0.65_r8
    kappa_aer_mac(ico3_a)  = 0.001_r8  ! ??
    kappa_aer_mac(ica_a)   = 0.001_r8  ! ??
    kappa_aer_mac(ioin_a)  = 0.001_r8
    kappa_aer_mac(ioc_a)   = 0.001_r8
    kappa_aer_mac(ibc_a)   = 0.001_r8
    kappa_aer_mac(iaro1_a) = 0.1_r8
    kappa_aer_mac(iaro2_a) = 0.1_r8
    kappa_aer_mac(ialk1_a) = 0.1_r8
    kappa_aer_mac(iole1_a) = 0.1_r8
    kappa_aer_mac(iapi1_a) = 0.1_r8
    kappa_aer_mac(iapi2_a) = 0.1_r8
    kappa_aer_mac(ilim1_a) = 0.1_r8
    kappa_aer_mac(ilim2_a) = 0.1_r8
    !BSINGH - 05/28/2013(RCE updates)
    if (ipmcmos_aero > 0) then
       ! use partmc-mosaic kappas
       kappa_aer_mac(1:19) = (/ &
            0.65_r8, 0.65_r8, 0.53_r8, 0.65_r8, 0.53_r8, 0.10_r8, 0.10_r8, 0.10_r8, 0.10_r8, 0.10_r8, &
            0.10_r8, 0.10_r8, 0.10_r8, 0.53_r8, 0.53_r8, 0.53_r8, 0.10_r8, 0.001_r8, 0.0_r8 /)
       !           so4   no3   cl    nh4   msa   aro1  aro2  alk1  ole1  api1
       !           api2  lim1  lim2  co3   na    ca    oin   oc    bc
    end if
    !BSINGH - 05/28/2013(RCE updates ENDS)

    ! partial molar volumes of condensing gases
    partial_molar_vol(1:ngas_aerchtot)   = 200.0_r8  ! default

    partial_molar_vol(ih2so4_g) = 51.83_r8
    partial_molar_vol(ihno3_g)  = 31.45_r8
    partial_molar_vol(ihcl_g)   = 20.96_r8
    partial_molar_vol(inh3_g)   = 24.03_r8
    partial_molar_vol(imsa_g)   = 53.33_r8
    partial_molar_vol(iaro1_g)  = 250.0_r8
    partial_molar_vol(iaro2_g)  = 250.0_r8
    partial_molar_vol(ialk1_g)  = 250.0_r8
    partial_molar_vol(iole1_g)  = 250.0_r8
    partial_molar_vol(iapi1_g)  = 250.0_r8
    partial_molar_vol(iapi2_g)  = 250.0_r8
    partial_molar_vol(ilim1_g)  = 250.0_r8
    partial_molar_vol(ilim2_g)  = 250.0_r8

    ! molecular weights of condensing gases
    mw_gas(1:ngas_aerchtot) = 200.0_r8  ! default
    mw_gas(ih2so4_g) = 98.0_r8
    mw_gas(ihno3_g)  = 63.0_r8
    mw_gas(ihcl_g)   = 36.5_r8
    mw_gas(inh3_g)   = 17.0_r8
    mw_gas(imsa_g)   = 96.0_r8
    mw_gas(iaro1_g)  = 250.0_r8
    mw_gas(iaro2_g)  = 250.0_r8
    mw_gas(ialk1_g)  = 250.0_r8
    mw_gas(iole1_g)  = 250.0_r8
    mw_gas(iapi1_g)  = 250.0_r8
    mw_gas(iapi2_g)  = 250.0_r8
    mw_gas(ilim1_g)  = 250.0_r8
    mw_gas(ilim2_g)  = 250.0_r8

    if ( use_cam5mam_soa_params > 0 ) then
       partial_molar_vol(ilim2_g)  = 250.0_r8
       mw_gas(ilim2_g)  = 250.0_r8
 !zlu soa
       partial_molar_vol(ilim1_g)  = 250.0_r8
       mw_gas(ilim1_g)  = 250.0_r8

       partial_molar_vol(iaro1_g)  = 250.0_r8
       mw_gas(iaro1_g)  = 250.0_r8

       partial_molar_vol(iapi1_g)  = 250.0_r8
       mw_gas(iapi1_g)  = 250.0_r8

       partial_molar_vol(iapi2_g)  = 250.0_r8
       mw_gas(iapi2_g)  = 250.0_r8



    endif

    ! used to calculate diffusivities of condensing gases
    v_molar_gas(1:ngas_aerchtot) = 60.0_r8  ! default
    v_molar_gas(ih2so4_g)= 42.88_r8
    v_molar_gas(ihno3_g) = 24.11_r8
    v_molar_gas(ihcl_g)  = 21.48_r8
    v_molar_gas(inh3_g)  = 14.90_r8
    v_molar_gas(imsa_g)  = 58.00_r8

    ! refractive index
    ref_index_a(jnh4so4) = cmplx(1.52_r8,0._r8)
    ref_index_a(jlvcite) = cmplx(1.50_r8,0._r8)
    ref_index_a(jnh4hso4)= cmplx(1.47_r8,0._r8)
    ref_index_a(jnh4msa) = cmplx(1.50_r8,0._r8)      ! assumed
    ref_index_a(jnh4no3) = cmplx(1.50_r8,0._r8)
    ref_index_a(jnh4cl)  = cmplx(1.50_r8,0._r8)
    ref_index_a(jnacl)   = cmplx(1.45_r8,0._r8)
    ref_index_a(jnano3)  = cmplx(1.50_r8,0._r8)
    ref_index_a(jna2so4) = cmplx(1.50_r8,0._r8)
    ref_index_a(jna3hso4)= cmplx(1.50_r8,0._r8)
    ref_index_a(jnahso4) = cmplx(1.50_r8,0._r8)
    ref_index_a(jnamsa)  = cmplx(1.50_r8,0._r8)      ! assumed
    ref_index_a(jcaso4)  = cmplx(1.56_r8,0.006_r8)
    ref_index_a(jcamsa2) = cmplx(1.56_r8,0.006_r8)   ! assumed
    ref_index_a(jcano3)  = cmplx(1.56_r8,0.006_r8)
    ref_index_a(jcacl2)  = cmplx(1.52_r8,0.006_r8)
    ref_index_a(jcaco3)  = cmplx(1.68_r8,0.006_r8)
    ref_index_a(jh2so4)  = cmplx(1.43_r8,0._r8)
    ref_index_a(jhhso4)  = cmplx(1.43_r8,0._r8)
    ref_index_a(jhno3)   = cmplx(1.50_r8,0._r8)
    ref_index_a(jhcl)    = cmplx(1.50_r8,0._r8)
    ref_index_a(jmsa)    = cmplx(1.43_r8,0._r8)      ! assumed
    ref_index_a(joc)      = cmplx(1.45_r8,0._r8)
    ref_index_a(jbc)      = cmplx(1.82_r8,0.74_r8)
    ref_index_a(join)    = cmplx(1.55_r8,0.006_r8)
    ref_index_a(jaro1)   = cmplx(1.45_r8,0._r8)
    ref_index_a(jaro2)   = cmplx(1.45_r8,0._r8)
    ref_index_a(jalk1)   = cmplx(1.45_r8,0._r8)
    ref_index_a(jole1)   = cmplx(1.45_r8,0._r8)
    ref_index_a(japi1)   = cmplx(1.45_r8,0._r8)
    ref_index_a(japi2)   = cmplx(1.45_r8,0._r8)
    ref_index_a(jlim1)   = cmplx(1.45_r8,0._r8)
    ref_index_a(jlim2)   = cmplx(1.45_r8,0._r8)
    ref_index_a(jh2o)    = cmplx(1.33_r8,0._r8)

    ! jsalt_index
    jsalt_index(jnh4so4) = 5           ! AS
    jsalt_index(jlvcite) = 2           ! LV
    jsalt_index(jnh4hso4)= 1           ! AB
    jsalt_index(jnh4no3) = 2           ! AN
    jsalt_index(jnh4cl)  = 1           ! AC
    jsalt_index(jna2so4) = 60          ! SS
    jsalt_index(jnahso4) = 10          ! SB
    jsalt_index(jnano3)  = 40          ! SN
    jsalt_index(jnacl)   = 10          ! SC
    jsalt_index(jcano3)  = 120 ! CN
    jsalt_index(jcacl2)  = 80          ! CC
    jsalt_index(jnh4msa) = 0           ! AM    zero for now
    jsalt_index(jnamsa)  = 0           ! SM    zero for now
    jsalt_index(jcamsa2) = 0           ! CM    zero for now

    ! Aerosol Indices
    !  AC = 1, AN = 2, AS = 5, SC = 10, SN = 40, SS = 60, CC = 80, CN = 120,
    !  AB = 1, LV = 2, SB = 10
    !
    ! SULFATE-POOR DOMAIN
    jsulf_poor(1)   =  1       !       AC
    jsulf_poor(2)   =  2       !       AN
    jsulf_poor(5)   =  3       !       AS
    jsulf_poor(10)  =  4       !       SC
    jsulf_poor(40)  =  5       !       SN
    jsulf_poor(60)  =  6       !       SS
    jsulf_poor(80)  =  7       !       CC
    jsulf_poor(120) =  8       !       CN
    jsulf_poor(3)   =  9       !       AN + AC
    jsulf_poor(6)   =  10      !       AS + AC
    jsulf_poor(7)   =  11      !       AS + AN
    jsulf_poor(8)   =          12      !       AS + AN + AC
    jsulf_poor(11)  =  13      !       SC + AC
    jsulf_poor(41)  =  14      !       SN + AC
    jsulf_poor(42)  =  15      !       SN + AN
    jsulf_poor(43)  =  16      !       SN + AN + AC
    jsulf_poor(50)  =  17      !       SN + SC
    jsulf_poor(51)  =  18      !       SN + SC + AC
    jsulf_poor(61)  =  19      !       SS + AC
    jsulf_poor(62)  =  20      !       SS + AN
    jsulf_poor(63)  =  21      !       SS + AN + AC
    jsulf_poor(65)  =  22      !       SS + AS
    jsulf_poor(66)  =  23      !       SS + AS + AC
    jsulf_poor(67)  =  24      !       SS + AS + AN
    jsulf_poor(68)  =  25      !       SS + AS + AN + AC
    jsulf_poor(70)  =  26      !       SS + SC
    jsulf_poor(71)  =  27      !       SS + SC + AC
    jsulf_poor(100) =  28      !       SS + SN
    jsulf_poor(101) =  29      !       SS + SN + AC
    jsulf_poor(102) =  30      !       SS + SN + AN
    jsulf_poor(103) =  31      !       SS + SN + AN + AC
    jsulf_poor(110) =  32      !       SS + SN + SC
    jsulf_poor(111) =  33      !       SS + SN + SC + AC
    jsulf_poor(81)  =  34      !       CC + AC
    jsulf_poor(90)  =  35      !       CC + SC
    jsulf_poor(91)  =  36      !       CC + SC + AC
    jsulf_poor(121) =  37      !       CN + AC
    jsulf_poor(122) =  38      !       CN + AN
    jsulf_poor(123) =  39      !       CN + AN + AC
    jsulf_poor(130) =  40      !       CN + SC
    jsulf_poor(131) =  41      !       CN + SC + AC
    jsulf_poor(160) =  42      !       CN + SN
    jsulf_poor(161) =  43      !       CN + SN + AC
    jsulf_poor(162) =  44      !       CN + SN + AN
    jsulf_poor(163) =  45      !       CN + SN + AN + AC
    jsulf_poor(170) =  46      !       CN + SN + SC
    jsulf_poor(171) =  47      !       CN + SN + SC + AC
    jsulf_poor(200) =  48      !       CN + CC
    jsulf_poor(201) =  49      !       CN + CC + AC
    jsulf_poor(210) =  50      !       CN + CC + SC
    jsulf_poor(211) =  51      !       CN + CC + SC + AC
    !
    ! SULFATE-RICH DOMAIN
    jsulf_rich(1)   =  52      !       AB
    jsulf_rich(2)   =  53      !       LV
    jsulf_rich(10)  =  54      !       SB
    jsulf_rich(3)   =  55      !       AB + LV
    jsulf_rich(7)   =  56      !       AS + LV
    jsulf_rich(70)  =  57      !       SS + SB
    jsulf_rich(62)  =  58      !       SS + LV
    jsulf_rich(67)  =  59      !       SS + AS + LV
    jsulf_rich(61)  =  60      !       SS + AB
    jsulf_rich(63)  =  61      !       SS + LV + AB
    jsulf_rich(11)  =  62      !       SB + AB
    jsulf_rich(71)  =  63      !       SS + SB + AB
    jsulf_rich(5)   =  3       !       AS
    jsulf_rich(60)  =  6       !       SS
    jsulf_rich(65)  =  22      !       SS + AS



    !
    ! polynomial coefficients for binary molality (used in ZSR equation)
    !
    !
    ! a_zsr for aw < 0.97
    !
    ! (NH4)2SO4
    je = jnh4so4
    a_zsr(1,je)  =  1.30894_r8
    a_zsr(2,je)  = -7.09922_r8
    a_zsr(3,je)  =  20.62831_r8
    a_zsr(4,je)  = -32.19965_r8
    a_zsr(5,je)  =  25.17026_r8
    a_zsr(6,je)  = -7.81632_r8
    aw_min(je)   = 0.1_r8
    !
    ! (NH4)3H(SO4)2
    je = jlvcite
    a_zsr(1,je)  =  1.10725_r8
    a_zsr(2,je)  = -5.17978_r8
    a_zsr(3,je)  =  12.29534_r8
    a_zsr(4,je)  = -16.32545_r8
    a_zsr(5,je)  =  11.29274_r8
    a_zsr(6,je)  = -3.19164_r8
    aw_min(je)   = 0.1_r8
    !
    ! NH4HSO4
    je = jnh4hso4
    a_zsr(1,je)  =  1.15510_r8
    a_zsr(2,je)  = -3.20815_r8
    a_zsr(3,je)  =  2.71141_r8
    a_zsr(4,je)  =  2.01155_r8
    a_zsr(5,je)  = -4.71014_r8
    a_zsr(6,je)  =  2.04616_r8
    aw_min(je)   = 0.1_r8
    !
    ! NH4MSA (assumed same as NH4HSO4)
    je = jnh4msa
    a_zsr(1,je)  =  1.15510_r8
    a_zsr(2,je)  = -3.20815_r8
    a_zsr(3,je)  =  2.71141_r8
    a_zsr(4,je)  =  2.01155_r8
    a_zsr(5,je)  = -4.71014_r8
    a_zsr(6,je)  =  2.04616_r8
    aw_min(je)   = 0.1_r8
    !
    ! NH4NO3
    je = jnh4no3
    a_zsr(1,je)  =  0.43507_r8
    a_zsr(2,je)  =  6.38220_r8
    a_zsr(3,je)  = -30.19797_r8
    a_zsr(4,je)  =  53.36470_r8
    a_zsr(5,je)  = -43.44203_r8
    a_zsr(6,je)  =  13.46158_r8
    aw_min(je)   = 0.1_r8
    !
    ! NH4Cl: revised on Nov 13, 2003. based on Chan and Ha (1999) JGR.
    je = jnh4cl
    a_zsr(1,je)  =  0.45309_r8
    a_zsr(2,je)  =  2.65606_r8
    a_zsr(3,je)  = -14.7730_r8
    a_zsr(4,je)  =  26.2936_r8
    a_zsr(5,je)  = -20.5735_r8
    a_zsr(6,je)  =  5.94255_r8
    aw_min(je)   = 0.1_r8
    !
    ! NaCl
    je = jnacl
    a_zsr(1,je)  =  0.42922_r8
    a_zsr(2,je)  = -1.17718_r8
    a_zsr(3,je)  =  2.80208_r8
    a_zsr(4,je)  = -4.51097_r8
    a_zsr(5,je)  =  3.76963_r8
    a_zsr(6,je)  = -1.31359_r8
    aw_min(je)   = 0.1_r8
    !
    ! NaNO3
    je = jnano3
    a_zsr(1,je)  =  1.34966_r8
    a_zsr(2,je)  = -5.20116_r8
    a_zsr(3,je)  =  11.49011_r8
    a_zsr(4,je)  = -14.41380_r8
    a_zsr(5,je)  =  9.07037_r8
    a_zsr(6,je)  = -2.29769_r8
    aw_min(je)   = 0.1_r8
    !
    ! Na2SO4
    je = jna2so4
    a_zsr(1,je)  =  0.39888_r8
    a_zsr(2,je)  = -1.27150_r8
    a_zsr(3,je)  =  3.42792_r8
    a_zsr(4,je)  = -5.92632_r8
    a_zsr(5,je)  =  5.33351_r8
    a_zsr(6,je)  = -1.96541_r8
    aw_min(je)   = 0.1_r8
    !
    ! Na3H(SO4)2  added on 1/14/2004
    je = jna3hso4
    a_zsr(1,je)  =  0.31480_r8
    a_zsr(2,je)  = -1.01087_r8
    a_zsr(3,je)  =  2.44029_r8
    a_zsr(4,je)  = -3.66095_r8
    a_zsr(5,je)  =  2.77632_r8
    a_zsr(6,je)  = -0.86058_r8
    aw_min(je)   = 0.1_r8
    !
    ! NaHSO4
    je = jnahso4
    a_zsr(1,je)  =  0.62764_r8
    a_zsr(2,je)  = -1.63520_r8
    a_zsr(3,je)  =  4.62531_r8
    a_zsr(4,je)  = -10.06925_r8
    a_zsr(5,je)  =  10.33547_r8
    a_zsr(6,je)  = -3.88729_r8
    aw_min(je)   = 0.1_r8
    !
    ! NaMSA (assumed same as NaHSO4)
    je = jnamsa
    a_zsr(1,je)  =  0.62764_r8
    a_zsr(2,je)  = -1.63520_r8
    a_zsr(3,je)  =  4.62531_r8
    a_zsr(4,je)  = -10.06925_r8
    a_zsr(5,je)  =  10.33547_r8
    a_zsr(6,je)  = -3.88729_r8
    aw_min(je)   = 0.1_r8
    !
    ! Ca(NO3)2
    je = jcano3
    a_zsr(1,je)  =  0.38895_r8
    a_zsr(2,je)  = -1.16013_r8
    a_zsr(3,je)  =  2.16819_r8
    a_zsr(4,je)  = -2.23079_r8
    a_zsr(5,je)  =  1.00268_r8
    a_zsr(6,je)  = -0.16923_r8
    aw_min(je)   = 0.1_r8
    !
    ! CaCl2: Kim and Seinfeld
    je = jcacl2
    a_zsr(1,je)  =  0.29891_r8
    a_zsr(2,je)  = -1.31104_r8
    a_zsr(3,je)  =  3.68759_r8
    a_zsr(4,je)  = -5.81708_r8
    a_zsr(5,je)  =  4.67520_r8
    a_zsr(6,je)  = -1.53223_r8
    aw_min(je)   = 0.1_r8
    !
    ! H2SO4
    je = jh2so4
    a_zsr(1,je) =  0.32751_r8
    a_zsr(2,je) = -1.00692_r8
    a_zsr(3,je) =  2.59750_r8
    a_zsr(4,je) = -4.40014_r8
    a_zsr(5,je) =  3.88212_r8
    a_zsr(6,je) = -1.39916_r8
    aw_min(je)  = 0.1_r8
    !
    ! MSA (assumed same as H2SO4)
    je = jmsa
    a_zsr(1,je) =  0.32751_r8
    a_zsr(2,je) = -1.00692_r8
    a_zsr(3,je) =  2.59750_r8
    a_zsr(4,je) = -4.40014_r8
    a_zsr(5,je) =  3.88212_r8
    a_zsr(6,je) = -1.39916_r8
    aw_min(je)  = 0.1_r8
    !
    ! HHSO4
    je = jhhso4
    a_zsr(1,je) =  0.32751_r8
    a_zsr(2,je) = -1.00692_r8
    a_zsr(3,je) =  2.59750_r8
    a_zsr(4,je) = -4.40014_r8
    a_zsr(5,je) =  3.88212_r8
    a_zsr(6,je) = -1.39916_r8
    aw_min(je)  = 1.0_r8
    !
    ! HNO3
    je = jhno3
    a_zsr(1,je) =  0.75876_r8
    a_zsr(2,je) = -3.31529_r8
    a_zsr(3,je) =  9.26392_r8
    a_zsr(4,je) = -14.89799_r8
    a_zsr(5,je) =  12.08781_r8
    a_zsr(6,je) = -3.89958_r8
    aw_min(je)  = 0.1_r8
    !
    ! HCl
    je = jhcl
    a_zsr(1,je) =  0.31133_r8
    a_zsr(2,je) = -0.79688_r8
    a_zsr(3,je) =  1.93995_r8
    a_zsr(4,je) = -3.31582_r8
    a_zsr(5,je) =  2.93513_r8
    a_zsr(6,je) = -1.07268_r8
    aw_min(je)  = 0.1_r8
    !
    ! CaSO4
    je = jcaso4
    a_zsr(1,je)  =  0.0_r8
    a_zsr(2,je)  =  0.0_r8
    a_zsr(3,je)  =  0.0_r8
    a_zsr(4,je)  =  0.0_r8
    a_zsr(5,je)  =  0.0_r8
    a_zsr(6,je)  =  0.0_r8
    aw_min(je)   = 1.0_r8
    !
    ! Ca(MSA)2 (assumed same as Ca(NO3)2)
    je = jcamsa2
    a_zsr(1,je)  =  0.38895_r8
    a_zsr(2,je)  = -1.16013_r8
    a_zsr(3,je)  =  2.16819_r8
    a_zsr(4,je)  = -2.23079_r8
    a_zsr(5,je)  =  1.00268_r8
    a_zsr(6,je)  = -0.16923_r8
    aw_min(je)   = 0.1_r8
    !
    ! CaCO3
    je = jcaco3
    a_zsr(1,je)  =  0.0_r8
    a_zsr(2,je)  =  0.0_r8
    a_zsr(3,je)  =  0.0_r8
    a_zsr(4,je)  =  0.0_r8
    a_zsr(5,je)  =  0.0_r8
    a_zsr(6,je)  =  0.0_r8
    aw_min(je)   = 1.0_r8



    !-------------------------------------------
    ! b_zsr for aw => 0.97 to 0.99999
    !
    ! (NH4)2SO4
    b_zsr(jnh4so4)  = 28.0811_r8
    !
    ! (NH4)3H(SO4)2
    b_zsr(jlvcite)  = 14.7178_r8
    !
    ! NH4HSO4
    b_zsr(jnh4hso4) = 29.4779_r8
    !
    ! NH4MSA
    b_zsr(jnh4msa)  = 29.4779_r8 ! assumed same as NH4HSO4
    !
    ! NH4NO3
    b_zsr(jnh4no3)  = 33.4049_r8
    !
    ! NH4Cl
    b_zsr(jnh4cl)   = 30.8888_r8
    !
    ! NaCl
    b_zsr(jnacl)    = 29.8375_r8
    !
    ! NaNO3
    b_zsr(jnano3)   = 32.2756_r8
    !
    ! Na2SO4
    b_zsr(jna2so4)  = 27.6889_r8
    !
    ! Na3H(SO4)2
    b_zsr(jna3hso4) = 14.2184_r8
    !
    ! NaHSO4
    b_zsr(jnahso4)  = 28.3367_r8
    !
    ! NaMSA
    b_zsr(jnamsa)   = 28.3367_r8 ! assumed same as NaHSO4
    !
    ! Ca(NO3)2
    b_zsr(jcano3)   = 18.3661_r8
    !
    ! CaCl2
    b_zsr(jcacl2)   = 20.8792_r8
    !
    ! H2SO4
    b_zsr(jh2so4)   = 26.7347_r8
    !
    ! HHSO4
    b_zsr(jhhso4)   = 26.7347_r8
    !
    ! HNO3
    b_zsr(jhno3)    = 28.8257_r8
    !
    ! HCl
    b_zsr(jhcl)     = 27.7108_r8
    !
    ! MSA
    b_zsr(jmsa)     = 26.7347_r8 ! assumed same as H2SO4
    !
    ! CaSO4
    b_zsr(jcaso4)   = 0.0_r8
    !
    ! Ca(MSA)2
    b_zsr(jcamsa2)  = 18.3661_r8 ! assumed same as Ca(NO3)2
    !
    ! CaCO3
    b_zsr(jcaco3)   = 0.0_r8









    !-------------------------------------------
    ! Li and Lu (2001) Surface tension model
    ! G_MX [mol/cm^2]; K_MX [-]
    !
    ! (NH4)2SO4
    G_MX(jnh4so4)  = -8.79e-7_r8*1.e-4_r8
    K_MX(jnh4so4)  =  3.84e+1_r8
    !
    ! (NH4)3H(SO4)2
    G_MX(jlvcite)  = -8.79e-7_r8*1.e-4_r8    ! assumed same as (NH4)2SO4
    K_MX(jlvcite)  =  3.84e+1_r8          ! assumed same as (NH4)2SO4
    !
    ! NH4HSO4
    G_MX(jnh4hso4) = -8.79e-7_r8*1.e-4_r8    ! assumed same as (NH4)2SO4
    K_MX(jnh4hso4) =  3.84e+1_r8          ! assumed same as (NH4)2SO4
    !
    ! NH4MSA
    G_MX(jnh4msa)  = -8.79e-7_r8*1.e-4_r8    ! assumed same as (NH4)2SO4
    K_MX(jnh4msa)  =  3.84e+1_r8          ! assumed same as (NH4)2SO4
    !
    ! NH4NO3
    G_MX(jnh4no3)  = -3.08e-6_r8*1.e-4_r8
    K_MX(jnh4no3)  =  4.89e-1_r8
    !
    ! NH4Cl
    G_MX(jnh4cl)   = -1.01e-6_r8*1.e-4_r8
    K_MX(jnh4cl)   =  1.3_r8
    !
    ! NaCl
    G_MX(jnacl)    = -1.05e-6_r8*1.e-4_r8
    K_MX(jnacl)    =  1.2_r8
    !
    ! NaNO3
    G_MX(jnano3)   = -1.66e-6_r8*1.e-4_r8
    K_MX(jnano3)   =  1.25_r8
    !
    ! Na2SO4
    G_MX(jna2so4)  = -8.37e-7_r8*1.e-4_r8
    K_MX(jna2so4)  =  7.57e+1_r8
    !
    ! Na3H(SO4)2
    G_MX(jna3hso4) = -8.37e-7_r8*1.e-4_r8    ! assumed same as Na2SO4
    K_MX(jna3hso4) =  7.57e+1_r8          ! assumed same as Na2SO4
    !
    ! NaHSO4
    G_MX(jnahso4)  = -8.37e-7_r8*1.e-4_r8    ! assumed same as Na2SO4
    K_MX(jnahso4)  =  7.57e+1_r8          ! assumed same as Na2SO4
    !
    ! NaMSA
    G_MX(jnamsa)   = -8.37e-7_r8*1.e-4_r8
    K_MX(jnamsa)   =  7.57e+1_r8
    !
    ! Ca(NO3)2
    G_MX(jcano3)   = -4.88e-7_r8*1.e-4_r8    ! assumed same as CaCl2
    K_MX(jcano3)   =  1.50e+1_r8          ! assumed same as CaCl2
    !
    ! CaCl2
    G_MX(jcacl2)   = -4.88e-7_r8*1.e-4_r8
    K_MX(jcacl2)   =  1.50e+1_r8
    !
    ! H2SO4
    G_MX(jh2so4)   = -6.75e-8_r8*1.e-4_r8
    K_MX(jh2so4)   =  1.65e+3_r8
    !
    ! HHSO4
    G_MX(jh2so4)   = -6.75e-8_r8*1.e-4_r8    ! assumed same as H2SO4
    K_MX(jh2so4)   =  1.65e+3_r8          ! assumed same as H2SO4
    !
    ! HNO3
    G_MX(jhno3)    =  8.05e-7_r8*1.e-4_r8
    K_MX(jhno3)    =  1.06e-1_r8
    !
    ! HCl
    G_MX(jhcl)     =  4.12e-7_r8*1.e-4_r8
    K_MX(jhcl)     =  4.68e-3_r8
    !&

    ! MSA
    G_MX(jmsa)     =  8.05e-7_r8*1.e-4_r8    ! assumed same as HNO3
    K_MX(jmsa)     =  1.06e-1_r8          ! assumed same as HNO3
    !
    ! CaSO4
    G_MX(jmsa)     =  0.0_r8*1.e-4_r8        ! assumed
    K_MX(jmsa)     =  0.0_r8              ! assumed
    !
    ! Ca(MSA)2
    G_MX(jcamsa2)  =  0.0_r8*1.e-4_r8        ! assumed
    K_MX(jcamsa2)  =  0.0_r8              ! assumed
    !
    ! CaCO3
    G_MX(jcaco3)   =  0.0_r8*1.e-4_r8        ! assumed
    K_MX(jcaco3)   =  0.0_r8              ! assumed







    !----------------------------------------------------------------
    ! parameters for MTEM mixing rule (Zaveri, Easter, and Wexler, 2005)
    ! log_gamZ(jA,jE)   A in E
    !----------------------------------------------------------------
    !
    b_mtem(:,:,:) = 0.0_r8 !BSINGH - Temporarily initialized, please modify if required *Ask dick about it* the code blows up if i initialize it with nan
    ! (NH4)2SO4 in E
    jA = jnh4so4

    ! in (NH4)2SO4
    jE = jnh4so4
    b_mtem(1,jA,jE) = -2.94685_r8
    b_mtem(2,jA,jE) = 17.3328_r8
    b_mtem(3,jA,jE) = -64.8441_r8
    b_mtem(4,jA,jE) = 122.7070_r8
    b_mtem(5,jA,jE) = -114.4373_r8
    b_mtem(6,jA,jE) = 41.6811_r8

    ! in NH4NO3
    jE = jnh4no3
    b_mtem(1,jA,jE) = -2.7503_r8
    b_mtem(2,jA,jE) = 4.3806_r8
    b_mtem(3,jA,jE) = -1.1110_r8
    b_mtem(4,jA,jE) = -1.7005_r8
    b_mtem(5,jA,jE) = -4.4207_r8
    b_mtem(6,jA,jE) = 5.1990_r8

    ! in NH4Cl (revised on 11/15/2003)
    jE = jnh4cl
    b_mtem(1,jA,jE) = -2.06952_r8
    b_mtem(2,jA,jE) = 7.1240_r8
    b_mtem(3,jA,jE) = -24.4274_r8
    b_mtem(4,jA,jE) = 51.1458_r8
    b_mtem(5,jA,jE) = -54.2056_r8
    b_mtem(6,jA,jE) = 22.0606_r8

    ! in Na2SO4
    jE = jna2so4
    b_mtem(1,jA,jE) = -2.17361_r8
    b_mtem(2,jA,jE) = 15.9919_r8
    b_mtem(3,jA,jE) = -69.0952_r8
    b_mtem(4,jA,jE) = 139.8860_r8
    b_mtem(5,jA,jE) = -134.9890_r8
    b_mtem(6,jA,jE) = 49.8877_r8

    ! in NaNO3
    jE = jnano3
    b_mtem(1,jA,jE) = -4.4370_r8
    b_mtem(2,jA,jE) = 24.0243_r8
    b_mtem(3,jA,jE) = -76.2437_r8
    b_mtem(4,jA,jE) = 128.6660_r8
    b_mtem(5,jA,jE) = -110.0900_r8
    b_mtem(6,jA,jE) = 37.7414_r8

    ! in NaCl
    jE = jnacl
    b_mtem(1,jA,jE) = -1.5394_r8
    b_mtem(2,jA,jE) = 5.8671_r8
    b_mtem(3,jA,jE) = -22.7726_r8
    b_mtem(4,jA,jE) = 47.0547_r8
    b_mtem(5,jA,jE) = -47.8266_r8
    b_mtem(6,jA,jE) = 18.8489_r8

    ! in HNO3
    jE = jhno3
    b_mtem(1,jA,jE) = -0.35750_r8
    b_mtem(2,jA,jE) = -3.82466_r8
    b_mtem(3,jA,jE) = 4.55462_r8
    b_mtem(4,jA,jE) = 5.05402_r8
    b_mtem(5,jA,jE) = -14.7476_r8
    b_mtem(6,jA,jE) = 8.8009_r8

    ! in HCl
    jE = jhcl
    b_mtem(1,jA,jE) = -2.15146_r8
    b_mtem(2,jA,jE) = 5.50205_r8
    b_mtem(3,jA,jE) = -19.1476_r8
    b_mtem(4,jA,jE) = 39.1880_r8
    b_mtem(5,jA,jE) = -39.9460_r8
    b_mtem(6,jA,jE) = 16.0700_r8

    ! in H2SO4
    jE = jh2so4
    b_mtem(1,jA,jE) = -2.52604_r8
    b_mtem(2,jA,jE) = 9.76022_r8
    b_mtem(3,jA,jE) = -35.2540_r8
    b_mtem(4,jA,jE) = 71.2981_r8
    b_mtem(5,jA,jE) = -71.8207_r8
    b_mtem(6,jA,jE) = 28.0758_r8

    !
    ! in NH4HSO4
    jE = jnh4hso4
    b_mtem(1,jA,jE) = -4.13219_r8
    b_mtem(2,jA,jE) = 13.8863_r8
    b_mtem(3,jA,jE) = -34.5387_r8
    b_mtem(4,jA,jE) = 56.5012_r8
    b_mtem(5,jA,jE) = -51.8702_r8
    b_mtem(6,jA,jE) = 19.6232_r8

    !
    ! in (NH4)3H(SO4)2
    jE = jlvcite
    b_mtem(1,jA,jE) = -2.53482_r8
    b_mtem(2,jA,jE) = 12.3333_r8
    b_mtem(3,jA,jE) = -46.1020_r8
    b_mtem(4,jA,jE) = 90.4775_r8
    b_mtem(5,jA,jE) = -88.1254_r8
    b_mtem(6,jA,jE) = 33.4715_r8

    !
    ! in NaHSO4
    jE = jnahso4
    b_mtem(1,jA,jE) = -3.23425_r8
    b_mtem(2,jA,jE) = 18.7842_r8
    b_mtem(3,jA,jE) = -78.7807_r8
    b_mtem(4,jA,jE) = 161.517_r8
    b_mtem(5,jA,jE) = -154.940_r8
    b_mtem(6,jA,jE) = 56.2252_r8

    !
    ! in Na3H(SO4)2
    jE = jna3hso4
    b_mtem(1,jA,jE) = -1.25316_r8
    b_mtem(2,jA,jE) = 7.40960_r8
    b_mtem(3,jA,jE) = -34.8929_r8
    b_mtem(4,jA,jE) = 72.8853_r8
    b_mtem(5,jA,jE) = -72.4503_r8
    b_mtem(6,jA,jE) = 27.7706_r8


    !-----------------
    ! NH4NO3 in E
    jA = jnh4no3

    ! in (NH4)2SO4
    jE = jnh4so4
    b_mtem(1,jA,jE) = -3.5201_r8
    b_mtem(2,jA,jE) = 21.6584_r8
    b_mtem(3,jA,jE) = -72.1499_r8
    b_mtem(4,jA,jE) = 126.7000_r8
    b_mtem(5,jA,jE) = -111.4550_r8
    b_mtem(6,jA,jE) = 38.5677_r8

    ! in NH4NO3
    jE = jnh4no3
    b_mtem(1,jA,jE) = -2.2630_r8
    b_mtem(2,jA,jE) = -0.1518_r8
    b_mtem(3,jA,jE) = 17.0898_r8
    b_mtem(4,jA,jE) = -36.7832_r8
    b_mtem(5,jA,jE) = 29.8407_r8
    b_mtem(6,jA,jE) = -7.9314_r8

    ! in NH4Cl (revised on 11/15/2003)
    jE = jnh4cl
    b_mtem(1,jA,jE) = -1.3851_r8
    b_mtem(2,jA,jE) = -0.4462_r8
    b_mtem(3,jA,jE) = 8.4567_r8
    b_mtem(4,jA,jE) = -11.5988_r8
    b_mtem(5,jA,jE) = 2.9802_r8
    b_mtem(6,jA,jE) = 1.8132_r8

    ! in Na2SO4
    jE = jna2so4
    b_mtem(1,jA,jE) = -1.7602_r8
    b_mtem(2,jA,jE) = 10.4044_r8
    b_mtem(3,jA,jE) = -35.5894_r8
    b_mtem(4,jA,jE) = 64.3584_r8
    b_mtem(5,jA,jE) = -57.8931_r8
    b_mtem(6,jA,jE) = 20.2141_r8

    ! in NaNO3
    jE = jnano3
    b_mtem(1,jA,jE) = -3.24346_r8
    b_mtem(2,jA,jE) = 16.2794_r8
    b_mtem(3,jA,jE) = -48.7601_r8
    b_mtem(4,jA,jE) = 79.2246_r8
    b_mtem(5,jA,jE) = -65.8169_r8
    b_mtem(6,jA,jE) = 22.1500_r8

    ! in NaCl
    jE = jnacl
    b_mtem(1,jA,jE) = -1.75658_r8
    b_mtem(2,jA,jE) = 7.71384_r8
    b_mtem(3,jA,jE) = -22.7984_r8
    b_mtem(4,jA,jE) = 39.1532_r8
    b_mtem(5,jA,jE) = -34.6165_r8
    b_mtem(6,jA,jE) = 12.1283_r8

    ! in Ca(NO3)2
    jE = jcano3
    b_mtem(1,jA,jE) = -0.97178_r8
    b_mtem(2,jA,jE) = 6.61964_r8
    b_mtem(3,jA,jE) = -26.2353_r8
    b_mtem(4,jA,jE) = 50.5259_r8
    b_mtem(5,jA,jE) = -47.6586_r8
    b_mtem(6,jA,jE) = 17.5074_r8

    ! in CaCl2 added on 12/22/2003
    jE = jcacl2
    b_mtem(1,jA,jE) = -0.41515_r8
    b_mtem(2,jA,jE) = 6.44101_r8
    b_mtem(3,jA,jE) = -26.4473_r8
    b_mtem(4,jA,jE) = 49.0718_r8
    b_mtem(5,jA,jE) = -44.2631_r8
    b_mtem(6,jA,jE) = 15.3771_r8

    ! in HNO3
    jE = jhno3
    b_mtem(1,jA,jE) = -1.20644_r8
    b_mtem(2,jA,jE) = 5.70117_r8
    b_mtem(3,jA,jE) = -18.2783_r8
    b_mtem(4,jA,jE) = 31.7199_r8
    b_mtem(5,jA,jE) = -27.8703_r8
    b_mtem(6,jA,jE) = 9.7299_r8

    ! in HCl
    jE = jhcl
    b_mtem(1,jA,jE) = -0.680862_r8
    b_mtem(2,jA,jE) = 3.59456_r8
    b_mtem(3,jA,jE) = -10.7969_r8
    b_mtem(4,jA,jE) = 17.8434_r8
    b_mtem(5,jA,jE) = -15.3165_r8
    b_mtem(6,jA,jE) = 5.17123_r8


    !----------
    ! NH4Cl in E
    jA = jnh4cl

    ! in (NH4)2SO4
    jE = jnh4so4
    b_mtem(1,jA,jE) = -2.8850_r8
    b_mtem(2,jA,jE) = 20.6970_r8
    b_mtem(3,jA,jE) = -70.6810_r8
    b_mtem(4,jA,jE) = 124.3690_r8
    b_mtem(5,jA,jE) = -109.2880_r8
    b_mtem(6,jA,jE) = 37.5831_r8

    ! in NH4NO3
    jE = jnh4no3
    b_mtem(1,jA,jE) = -1.9386_r8
    b_mtem(2,jA,jE) = 1.3238_r8
    b_mtem(3,jA,jE) = 11.8500_r8
    b_mtem(4,jA,jE) = -28.1168_r8
    b_mtem(5,jA,jE) = 21.8543_r8
    b_mtem(6,jA,jE) = -5.1671_r8

    ! in NH4Cl (revised on 11/15/2003)
    jE = jnh4cl
    b_mtem(1,jA,jE) = -0.9559_r8
    b_mtem(2,jA,jE) = 0.8121_r8
    b_mtem(3,jA,jE) = 4.3644_r8
    b_mtem(4,jA,jE) = -8.9258_r8
    b_mtem(5,jA,jE) = 4.2362_r8
    b_mtem(6,jA,jE) = 0.2891_r8

    ! in Na2SO4
    jE = jna2so4
    b_mtem(1,jA,jE) = 0.0377_r8
    b_mtem(2,jA,jE) = 6.0752_r8
    b_mtem(3,jA,jE) = -30.8641_r8
    b_mtem(4,jA,jE) = 63.3095_r8
    b_mtem(5,jA,jE) = -61.0070_r8
    b_mtem(6,jA,jE) = 22.1734_r8

    ! in NaNO3
    jE = jnano3
    b_mtem(1,jA,jE) = -1.8336_r8
    b_mtem(2,jA,jE) = 12.8160_r8
    b_mtem(3,jA,jE) = -42.3388_r8
    b_mtem(4,jA,jE) = 71.1816_r8
    b_mtem(5,jA,jE) = -60.5708_r8
    b_mtem(6,jA,jE) = 20.5853_r8

    ! in NaCl
    jE = jnacl
    b_mtem(1,jA,jE) = -0.1429_r8
    b_mtem(2,jA,jE) = 2.3561_r8
    b_mtem(3,jA,jE) = -10.4425_r8
    b_mtem(4,jA,jE) = 20.8951_r8
    b_mtem(5,jA,jE) = -20.7739_r8
    b_mtem(6,jA,jE) = 7.9355_r8

    ! in Ca(NO3)2
    jE = jcano3
    b_mtem(1,jA,jE) = 0.76235_r8
    b_mtem(2,jA,jE) = 3.08323_r8
    b_mtem(3,jA,jE) = -23.6772_r8
    b_mtem(4,jA,jE) = 53.7415_r8
    b_mtem(5,jA,jE) = -55.4043_r8
    b_mtem(6,jA,jE) = 21.2944_r8

    ! in CaCl2 (revised on 11/27/2003)
    jE = jcacl2
    b_mtem(1,jA,jE) = 1.13864_r8
    b_mtem(2,jA,jE) = -0.340539_r8
    b_mtem(3,jA,jE) = -8.67025_r8
    b_mtem(4,jA,jE) = 22.8008_r8
    b_mtem(5,jA,jE) = -24.5181_r8
    b_mtem(6,jA,jE) = 9.3663_r8

    ! in HNO3
    jE = jhno3
    b_mtem(1,jA,jE) = 2.42532_r8
    b_mtem(2,jA,jE) = -14.1755_r8
    b_mtem(3,jA,jE) = 38.804_r8
    b_mtem(4,jA,jE) = -58.2437_r8
    b_mtem(5,jA,jE) = 43.5431_r8
    b_mtem(6,jA,jE) = -12.5824_r8

    ! in HCl
    jE = jhcl
    b_mtem(1,jA,jE) = 0.330337_r8
    b_mtem(2,jA,jE) = 0.0778934_r8
    b_mtem(3,jA,jE) = -2.30492_r8
    b_mtem(4,jA,jE) = 4.73003_r8
    b_mtem(5,jA,jE) = -4.80849_r8
    b_mtem(6,jA,jE) = 1.78866_r8



    !----------
    ! Na2SO4 in E
    jA = jna2so4

    ! in (NH4)2SO4
    jE = jnh4so4
    b_mtem(1,jA,jE) = -2.6982_r8
    b_mtem(2,jA,jE) = 22.9875_r8
    b_mtem(3,jA,jE) = -98.9840_r8
    b_mtem(4,jA,jE) = 198.0180_r8
    b_mtem(5,jA,jE) = -188.7270_r8
    b_mtem(6,jA,jE) = 69.0548_r8

    ! in NH4NO3
    jE = jnh4no3
    b_mtem(1,jA,jE) = -2.4844_r8
    b_mtem(2,jA,jE) = 6.5420_r8
    b_mtem(3,jA,jE) = -9.8998_r8
    b_mtem(4,jA,jE) = 11.3884_r8
    b_mtem(5,jA,jE) = -13.6842_r8
    b_mtem(6,jA,jE) = 7.7411_r8

    ! in NH4Cl (revised on 11/15/2003)
    jE = jnh4cl
    b_mtem(1,jA,jE) = -1.3325_r8
    b_mtem(2,jA,jE) = 13.0406_r8
    b_mtem(3,jA,jE) = -56.1935_r8
    b_mtem(4,jA,jE) = 107.1170_r8
    b_mtem(5,jA,jE) = -97.3721_r8
    b_mtem(6,jA,jE) = 34.3763_r8

    ! in Na2SO4
    jE = jna2so4
    b_mtem(1,jA,jE) = -1.2832_r8
    b_mtem(2,jA,jE) = 12.8526_r8
    b_mtem(3,jA,jE) = -62.2087_r8
    b_mtem(4,jA,jE) = 130.3876_r8
    b_mtem(5,jA,jE) = -128.2627_r8
    b_mtem(6,jA,jE) = 48.0340_r8

    ! in NaNO3
    jE = jnano3
    b_mtem(1,jA,jE) = -3.5384_r8
    b_mtem(2,jA,jE) = 21.3758_r8
    b_mtem(3,jA,jE) = -70.7638_r8
    b_mtem(4,jA,jE) = 121.1580_r8
    b_mtem(5,jA,jE) = -104.6230_r8
    b_mtem(6,jA,jE) = 36.0557_r8


    ! in NaCl
    jE = jnacl
    b_mtem(1,jA,jE) = 0.2175_r8
    b_mtem(2,jA,jE) = -0.5648_r8
    b_mtem(3,jA,jE) = -8.0288_r8
    b_mtem(4,jA,jE) = 25.9734_r8
    b_mtem(5,jA,jE) = -32.3577_r8
    b_mtem(6,jA,jE) = 14.3924_r8

    ! in HNO3
    jE = jhno3
    b_mtem(1,jA,jE) = -0.309617_r8
    b_mtem(2,jA,jE) = -1.82899_r8
    b_mtem(3,jA,jE) = -1.5505_r8
    b_mtem(4,jA,jE) = 13.3847_r8
    b_mtem(5,jA,jE) = -20.1284_r8
    b_mtem(6,jA,jE) = 9.93163_r8

    ! in HCl
    jE = jhcl
    b_mtem(1,jA,jE) = -0.259455_r8
    b_mtem(2,jA,jE) = -0.819366_r8
    b_mtem(3,jA,jE) = -4.28964_r8
    b_mtem(4,jA,jE) = 16.4305_r8
    b_mtem(5,jA,jE) = -21.8546_r8
    b_mtem(6,jA,jE) = 10.3044_r8

    ! in H2SO4
    jE = jh2so4
    b_mtem(1,jA,jE) = -1.84257_r8
    b_mtem(2,jA,jE) = 7.85788_r8
    b_mtem(3,jA,jE) = -29.9275_r8
    b_mtem(4,jA,jE) = 61.7515_r8
    b_mtem(5,jA,jE) = -63.2308_r8
    b_mtem(6,jA,jE) = 24.9542_r8

    ! in NH4HSO4
    jE = jnh4hso4
    b_mtem(1,jA,jE) = -1.05891_r8
    b_mtem(2,jA,jE) = 2.84831_r8
    b_mtem(3,jA,jE) = -21.1827_r8
    b_mtem(4,jA,jE) = 57.5175_r8
    b_mtem(5,jA,jE) = -64.8120_r8
    b_mtem(6,jA,jE) = 26.1986_r8

    ! in (NH4)3H(SO4)2
    jE = jlvcite
    b_mtem(1,jA,jE) = -1.16584_r8
    b_mtem(2,jA,jE) = 8.50075_r8
    b_mtem(3,jA,jE) = -44.3420_r8
    b_mtem(4,jA,jE) = 97.3974_r8
    b_mtem(5,jA,jE) = -98.4549_r8
    b_mtem(6,jA,jE) = 37.6104_r8

    ! in NaHSO4
    jE = jnahso4
    b_mtem(1,jA,jE) = -1.95805_r8
    b_mtem(2,jA,jE) = 6.62417_r8
    b_mtem(3,jA,jE) = -31.8072_r8
    b_mtem(4,jA,jE) = 77.8603_r8
    b_mtem(5,jA,jE) = -84.6458_r8
    b_mtem(6,jA,jE) = 33.4963_r8

    ! in Na3H(SO4)2
    jE = jna3hso4
    b_mtem(1,jA,jE) = -0.36045_r8
    b_mtem(2,jA,jE) = 3.55223_r8
    b_mtem(3,jA,jE) = -24.0327_r8
    b_mtem(4,jA,jE) = 54.4879_r8
    b_mtem(5,jA,jE) = -56.6531_r8
    b_mtem(6,jA,jE) = 22.4956_r8


    !----------
    ! NaNO3 in E
    jA = jnano3

    ! in (NH4)2SO4
    jE = jnh4so4
    b_mtem(1,jA,jE) = -2.5888_r8
    b_mtem(2,jA,jE) = 17.6192_r8
    b_mtem(3,jA,jE) = -63.2183_r8
    b_mtem(4,jA,jE) = 115.3520_r8
    b_mtem(5,jA,jE) = -104.0860_r8
    b_mtem(6,jA,jE) = 36.7390_r8

    ! in NH4NO3
    jE = jnh4no3

    b_mtem(1,jA,jE) = -2.0669_r8
    b_mtem(2,jA,jE) = 1.4792_r8
    b_mtem(3,jA,jE) = 10.5261_r8
    b_mtem(4,jA,jE) = -27.0987_r8
    b_mtem(5,jA,jE) = 23.0591_r8
    b_mtem(6,jA,jE) = -6.0938_r8

    ! in NH4Cl (revised on 11/15/2003)
    jE = jnh4cl
    b_mtem(1,jA,jE) = -0.8325_r8
    b_mtem(2,jA,jE) = 3.9933_r8
    b_mtem(3,jA,jE) = -15.3789_r8
    b_mtem(4,jA,jE) = 30.4050_r8
    b_mtem(5,jA,jE) = -29.4204_r8
    b_mtem(6,jA,jE) = 11.0597_r8

    ! in Na2SO4
    jE = jna2so4
    b_mtem(1,jA,jE) = -1.1233_r8
    b_mtem(2,jA,jE) = 8.3998_r8
    b_mtem(3,jA,jE) = -31.9002_r8
    b_mtem(4,jA,jE) = 60.1450_r8
    b_mtem(5,jA,jE) = -55.5503_r8
    b_mtem(6,jA,jE) = 19.7757_r8

    ! in NaNO3
    jE = jnano3
    b_mtem(1,jA,jE) = -2.5386_r8
    b_mtem(2,jA,jE) = 13.9039_r8
    b_mtem(3,jA,jE) = -42.8467_r8
    b_mtem(4,jA,jE) = 69.7442_r8
    b_mtem(5,jA,jE) = -57.8988_r8
    b_mtem(6,jA,jE) = 19.4635_r8

    ! in NaCl
    jE = jnacl
    b_mtem(1,jA,jE) = -0.4351_r8
    b_mtem(2,jA,jE) = 2.8311_r8
    b_mtem(3,jA,jE) = -11.4485_r8
    b_mtem(4,jA,jE) = 22.7201_r8
    b_mtem(5,jA,jE) = -22.4228_r8
    b_mtem(6,jA,jE) = 8.5792_r8

    ! in Ca(NO3)2
    jE = jcano3
    b_mtem(1,jA,jE) = -0.72060_r8
    b_mtem(2,jA,jE) = 5.64915_r8
    b_mtem(3,jA,jE) = -23.5020_r8
    b_mtem(4,jA,jE) = 46.0078_r8
    b_mtem(5,jA,jE) = -43.8075_r8
    b_mtem(6,jA,jE) = 16.1652_r8

    ! in CaCl2
    jE = jcacl2

    b_mtem(1,jA,jE) = 0.003928_r8
    b_mtem(2,jA,jE) = 3.54724_r8
    b_mtem(3,jA,jE) = -18.6057_r8
    b_mtem(4,jA,jE) = 38.1445_r8
    b_mtem(5,jA,jE) = -36.7745_r8
    b_mtem(6,jA,jE) = 13.4529_r8

    ! in HNO3
    jE = jhno3
    b_mtem(1,jA,jE) = -1.1712_r8
    b_mtem(2,jA,jE) = 7.20907_r8
    b_mtem(3,jA,jE) = -22.9215_r8
    b_mtem(4,jA,jE) = 38.1257_r8
    b_mtem(5,jA,jE) = -32.0759_r8
    b_mtem(6,jA,jE) = 10.6443_r8

    ! in HCl
    jE = jhcl
    b_mtem(1,jA,jE) = 0.738022_r8
    b_mtem(2,jA,jE) = -1.14313_r8
    b_mtem(3,jA,jE) = 0.32251_r8
    b_mtem(4,jA,jE) = 0.838679_r8
    b_mtem(5,jA,jE) = -1.81747_r8
    b_mtem(6,jA,jE) = 0.873986_r8


    !----------
    ! NaCl in E
    jA = jnacl

    ! in (NH4)2SO4
    jE = jnh4so4
    b_mtem(1,jA,jE) = -1.9525_r8
    b_mtem(2,jA,jE) = 16.6433_r8
    b_mtem(3,jA,jE) = -61.7090_r8
    b_mtem(4,jA,jE) = 112.9910_r8
    b_mtem(5,jA,jE) = -101.9370_r8
    b_mtem(6,jA,jE) = 35.7760_r8

    ! in NH4NO3
    jE = jnh4no3
    b_mtem(1,jA,jE) = -1.7525_r8
    b_mtem(2,jA,jE) = 3.0713_r8
    b_mtem(3,jA,jE) = 4.8063_r8
    b_mtem(4,jA,jE) = -17.5334_r8
    b_mtem(5,jA,jE) = 14.2872_r8
    b_mtem(6,jA,jE) = -3.0690_r8

    ! in NH4Cl (revised on 11/15/2003)
    jE = jnh4cl
    b_mtem(1,jA,jE) = -0.4021_r8
    b_mtem(2,jA,jE) = 5.2399_r8
    b_mtem(3,jA,jE) = -19.4278_r8
    b_mtem(4,jA,jE) = 33.0027_r8
    b_mtem(5,jA,jE) = -28.1020_r8
    b_mtem(6,jA,jE) = 9.5159_r8

    ! in Na2SO4
    jE = jna2so4
    b_mtem(1,jA,jE) = 0.6692_r8
    b_mtem(2,jA,jE) = 4.1207_r8
    b_mtem(3,jA,jE) = -27.3314_r8
    b_mtem(4,jA,jE) = 59.3112_r8
    b_mtem(5,jA,jE) = -58.7998_r8
    b_mtem(6,jA,jE) = 21.7674_r8

    ! in NaNO3
    jE = jnano3
    b_mtem(1,jA,jE) = -1.17444_r8
    b_mtem(2,jA,jE) = 10.9927_r8
    b_mtem(3,jA,jE) = -38.9013_r8
    b_mtem(4,jA,jE) = 66.8521_r8
    b_mtem(5,jA,jE) = -57.6564_r8
    b_mtem(6,jA,jE) = 19.7296_r8

    ! in NaCl
    jE = jnacl
    b_mtem(1,jA,jE) = 1.17679_r8
    b_mtem(2,jA,jE) = -2.5061_r8
    b_mtem(3,jA,jE) = 0.8508_r8
    b_mtem(4,jA,jE) = 4.4802_r8
    b_mtem(5,jA,jE) = -8.4945_r8
    b_mtem(6,jA,jE) = 4.3182_r8

    ! in Ca(NO3)2
    jE = jcano3
    b_mtem(1,jA,jE) = 1.01450_r8
    b_mtem(2,jA,jE) = 2.10260_r8
    b_mtem(3,jA,jE) = -20.9036_r8
    b_mtem(4,jA,jE) = 49.1481_r8
    b_mtem(5,jA,jE) = -51.4867_r8
    b_mtem(6,jA,jE) = 19.9301_r8

    ! in CaCl2 (PSC92: revised on 11/27/2003)
    jE = jcacl2
    b_mtem(1,jA,jE) = 1.55463_r8
    b_mtem(2,jA,jE) = -3.20122_r8
    b_mtem(3,jA,jE) = -0.957075_r8
    b_mtem(4,jA,jE) = 12.103_r8
    b_mtem(5,jA,jE) = -17.221_r8
    b_mtem(6,jA,jE) = 7.50264_r8

    ! in HNO3
    jE = jhno3
    b_mtem(1,jA,jE) = 2.46187_r8
    b_mtem(2,jA,jE) = -12.6845_r8
    b_mtem(3,jA,jE) = 34.2383_r8
    b_mtem(4,jA,jE) = -51.9992_r8
    b_mtem(5,jA,jE) = 39.4934_r8
    b_mtem(6,jA,jE) = -11.7247_r8

    ! in HCl
    jE = jhcl
    b_mtem(1,jA,jE) = 1.74915_r8
    b_mtem(2,jA,jE) = -4.65768_r8
    b_mtem(3,jA,jE) = 8.80287_r8
    b_mtem(4,jA,jE) = -12.2503_r8
    b_mtem(5,jA,jE) = 8.668751_r8
    b_mtem(6,jA,jE) = -2.50158_r8


    !----------
    ! Ca(NO3)2 in E
    jA = jcano3

    ! in NH4NO3
    jE = jnh4no3
    b_mtem(1,jA,jE) = -1.86260_r8
    b_mtem(2,jA,jE) = 11.6178_r8
    b_mtem(3,jA,jE) = -30.9069_r8
    b_mtem(4,jA,jE) = 41.7578_r8
    b_mtem(5,jA,jE) = -33.7338_r8
    b_mtem(6,jA,jE) = 12.7541_r8

    ! in NH4Cl (revised on 11/15/2003)
    jE = jnh4cl
    b_mtem(1,jA,jE) = -1.1798_r8
    b_mtem(2,jA,jE) = 25.9608_r8
    b_mtem(3,jA,jE) = -98.9373_r8
    b_mtem(4,jA,jE) = 160.2300_r8
    b_mtem(5,jA,jE) = -125.9540_r8
    b_mtem(6,jA,jE) = 39.5130_r8

    ! in NaNO3
    jE = jnano3
    b_mtem(1,jA,jE) = -1.44384_r8
    b_mtem(2,jA,jE) = 13.6044_r8
    b_mtem(3,jA,jE) = -54.4300_r8
    b_mtem(4,jA,jE) = 100.582_r8
    b_mtem(5,jA,jE) = -91.2364_r8
    b_mtem(6,jA,jE) = 32.5970_r8

    ! in NaCl
    jE = jnacl
    b_mtem(1,jA,jE) = -0.099114_r8
    b_mtem(2,jA,jE) = 2.84091_r8
    b_mtem(3,jA,jE) = -16.9229_r8
    b_mtem(4,jA,jE) = 37.4839_r8
    b_mtem(5,jA,jE) = -39.5132_r8
    b_mtem(6,jA,jE) = 15.8564_r8

    ! in Ca(NO3)2
    jE = jcano3
    b_mtem(1,jA,jE) = 0.055116_r8
    b_mtem(2,jA,jE) = 4.58610_r8
    b_mtem(3,jA,jE) = -27.6629_r8
    b_mtem(4,jA,jE) = 60.8288_r8
    b_mtem(5,jA,jE) = -61.4988_r8
    b_mtem(6,jA,jE) = 23.3136_r8

    ! in CaCl2 (PSC92: revised on 11/27/2003)
    jE = jcacl2
    b_mtem(1,jA,jE) = 1.57155_r8
    b_mtem(2,jA,jE) = -3.18486_r8
    b_mtem(3,jA,jE) = -3.35758_r8
    b_mtem(4,jA,jE) = 18.7501_r8
    b_mtem(5,jA,jE) = -24.5604_r8
    b_mtem(6,jA,jE) = 10.3798_r8

    ! in HNO3
    jE = jhno3
    b_mtem(1,jA,jE) = 1.04446_r8
    b_mtem(2,jA,jE) = -3.19066_r8
    b_mtem(3,jA,jE) = 2.44714_r8
    b_mtem(4,jA,jE) = 2.07218_r8
    b_mtem(5,jA,jE) = -6.43949_r8
    b_mtem(6,jA,jE) = 3.66471_r8

    ! in HCl
    jE = jhcl
    b_mtem(1,jA,jE) = 1.05723_r8
    b_mtem(2,jA,jE) = -1.46826_r8
    b_mtem(3,jA,jE) = -1.0713_r8
    b_mtem(4,jA,jE) = 4.64439_r8
    b_mtem(5,jA,jE) = -6.32402_r8
    b_mtem(6,jA,jE) = 2.78202_r8


    !----------
    ! CaCl2 in E
    jA = jcacl2

    ! in NH4NO3 (PSC92: revised on 12/22/2003)
    jE = jnh4no3
    b_mtem(1,jA,jE) = -1.43626_r8
    b_mtem(2,jA,jE) = 13.6598_r8
    b_mtem(3,jA,jE) = -38.2068_r8
    b_mtem(4,jA,jE) = 53.9057_r8
    b_mtem(5,jA,jE) = -44.9018_r8
    b_mtem(6,jA,jE) = 16.6120_r8

    ! in NH4Cl (PSC92: revised on 11/27/2003)
    jE = jnh4cl
    b_mtem(1,jA,jE) = -0.603965_r8
    b_mtem(2,jA,jE) = 27.6027_r8
    b_mtem(3,jA,jE) = -104.258_r8
    b_mtem(4,jA,jE) = 163.553_r8
    b_mtem(5,jA,jE) = -124.076_r8
    b_mtem(6,jA,jE) = 37.4153_r8

    ! in NaNO3 (PSC92: revised on 12/22/2003)
    jE = jnano3
    b_mtem(1,jA,jE) = 0.44648_r8
    b_mtem(2,jA,jE) = 8.8850_r8
    b_mtem(3,jA,jE) = -45.5232_r8
    b_mtem(4,jA,jE) = 89.3263_r8
    b_mtem(5,jA,jE) = -83.8604_r8
    b_mtem(6,jA,jE) = 30.4069_r8

    ! in NaCl (PSC92: revised on 11/27/2003)
    jE = jnacl
    b_mtem(1,jA,jE) = 1.61927_r8
    b_mtem(2,jA,jE) = 0.247547_r8
    b_mtem(3,jA,jE) = -18.1252_r8
    b_mtem(4,jA,jE) = 45.2479_r8
    b_mtem(5,jA,jE) = -48.6072_r8
    b_mtem(6,jA,jE) = 19.2784_r8

    ! in Ca(NO3)2 (PSC92: revised on 11/27/2003)
    jE = jcano3
    b_mtem(1,jA,jE) = 2.36667_r8
    b_mtem(2,jA,jE) = -0.123309_r8
    b_mtem(3,jA,jE) = -24.2723_r8
    b_mtem(4,jA,jE) = 65.1486_r8
    b_mtem(5,jA,jE) = -71.8504_r8
    b_mtem(6,jA,jE) = 28.3696_r8

    ! in CaCl2 (PSC92: revised on 11/27/2003)
    jE = jcacl2
    b_mtem(1,jA,jE) = 3.64023_r8
    b_mtem(2,jA,jE) = -12.1926_r8
    b_mtem(3,jA,jE) = 20.2028_r8
    b_mtem(4,jA,jE) = -16.0056_r8
    b_mtem(5,jA,jE) = 1.52355_r8
    b_mtem(6,jA,jE) = 2.44709_r8

    ! in HNO3
    jE = jhno3
    b_mtem(1,jA,jE) = 5.88794_r8
    b_mtem(2,jA,jE) = -29.7083_r8
    b_mtem(3,jA,jE) = 78.6309_r8
    b_mtem(4,jA,jE) = -118.037_r8
    b_mtem(5,jA,jE) = 88.932_r8
    b_mtem(6,jA,jE) = -26.1407_r8

    ! in HCl
    jE = jhcl
    b_mtem(1,jA,jE) = 2.40628_r8
    b_mtem(2,jA,jE) = -6.16566_r8
    b_mtem(3,jA,jE) = 10.2851_r8
    b_mtem(4,jA,jE) = -12.9035_r8
    b_mtem(5,jA,jE) = 7.7441_r8
    b_mtem(6,jA,jE) = -1.74821_r8


    !----------
    ! HNO3 in E
    jA = jhno3

    ! in (NH4)2SO4
    jE = jnh4so4
    b_mtem(1,jA,jE) = -3.57598_r8
    b_mtem(2,jA,jE) = 21.5469_r8
    b_mtem(3,jA,jE) = -77.4111_r8
    b_mtem(4,jA,jE) = 144.136_r8
    b_mtem(5,jA,jE) = -132.849_r8
    b_mtem(6,jA,jE) = 47.9412_r8

    ! in NH4NO3
    jE = jnh4no3
    b_mtem(1,jA,jE) = -2.00209_r8
    b_mtem(2,jA,jE) = -3.48399_r8
    b_mtem(3,jA,jE) = 34.9906_r8
    b_mtem(4,jA,jE) = -68.6653_r8
    b_mtem(5,jA,jE) = 54.0992_r8
    b_mtem(6,jA,jE) = -15.1343_r8

    ! in NH4Cl revised on 12/22/2003
    jE = jnh4cl
    b_mtem(1,jA,jE) = -0.63790_r8
    b_mtem(2,jA,jE) = -1.67730_r8
    b_mtem(3,jA,jE) = 10.1727_r8
    b_mtem(4,jA,jE) = -14.9097_r8
    b_mtem(5,jA,jE) = 7.67410_r8
    b_mtem(6,jA,jE) = -0.79586_r8

    ! in NaCl
    jE = jnacl
    b_mtem(1,jA,jE) = 1.3446_r8
    b_mtem(2,jA,jE) = -2.5578_r8
    b_mtem(3,jA,jE) = 1.3464_r8
    b_mtem(4,jA,jE) = 2.90537_r8
    b_mtem(5,jA,jE) = -6.53014_r8
    b_mtem(6,jA,jE) = 3.31339_r8

    ! in NaNO3
    jE = jnano3
    b_mtem(1,jA,jE) = -0.546636_r8
    b_mtem(2,jA,jE) = 10.3127_r8
    b_mtem(3,jA,jE) = -39.9603_r8
    b_mtem(4,jA,jE) = 71.4609_r8
    b_mtem(5,jA,jE) = -63.4958_r8
    b_mtem(6,jA,jE) = 22.0679_r8

    ! in Na2SO4
    jE = jna2so4
    b_mtem(1,jA,jE) = 1.35059_r8
    b_mtem(2,jA,jE) = 4.34557_r8
    b_mtem(3,jA,jE) = -35.8425_r8
    b_mtem(4,jA,jE) = 80.9868_r8
    b_mtem(5,jA,jE) = -81.6544_r8
    b_mtem(6,jA,jE) = 30.4841_r8

    ! in Ca(NO3)2
    jE = jcano3
    b_mtem(1,jA,jE) = 0.869414_r8
    b_mtem(2,jA,jE) = 2.98486_r8
    b_mtem(3,jA,jE) = -22.255_r8
    b_mtem(4,jA,jE) = 50.1863_r8
    b_mtem(5,jA,jE) = -51.214_r8
    b_mtem(6,jA,jE) = 19.2235_r8

    ! in CaCl2 (KM) revised on 12/22/2003
    jE = jcacl2
    b_mtem(1,jA,jE) = 1.42800_r8
    b_mtem(2,jA,jE) = -1.78959_r8
    b_mtem(3,jA,jE) = -2.49075_r8
    b_mtem(4,jA,jE) = 10.1877_r8
    b_mtem(5,jA,jE) = -12.1948_r8
    b_mtem(6,jA,jE) = 4.64475_r8

    ! in HNO3 (added on 12/06/2004)
    jE = jhno3
    b_mtem(1,jA,jE) = 0.22035_r8
    b_mtem(2,jA,jE) = 2.94973_r8
    b_mtem(3,jA,jE) = -12.1469_r8
    b_mtem(4,jA,jE) = 20.4905_r8
    b_mtem(5,jA,jE) = -17.3966_r8
    b_mtem(6,jA,jE) = 5.70779_r8

    ! in HCl (added on 12/06/2004)
    jE = jhcl
    b_mtem(1,jA,jE) = 1.55503_r8
    b_mtem(2,jA,jE) = -3.61226_r8
    b_mtem(3,jA,jE) = 6.28265_r8
    b_mtem(4,jA,jE) = -8.69575_r8
    b_mtem(5,jA,jE) = 6.09372_r8
    b_mtem(6,jA,jE) = -1.80898_r8

    ! in H2SO4
    jE = jh2so4
    b_mtem(1,jA,jE) = 1.10783_r8
    b_mtem(2,jA,jE) = -1.3363_r8
    b_mtem(3,jA,jE) = -1.83525_r8
    b_mtem(4,jA,jE) = 7.47373_r8
    b_mtem(5,jA,jE) = -9.72954_r8
    b_mtem(6,jA,jE) = 4.12248_r8

    ! in NH4HSO4
    jE = jnh4hso4
    b_mtem(1,jA,jE) = -0.851026_r8
    b_mtem(2,jA,jE) = 12.2515_r8
    b_mtem(3,jA,jE) = -49.788_r8
    b_mtem(4,jA,jE) = 91.6215_r8
    b_mtem(5,jA,jE) = -81.4877_r8
    b_mtem(6,jA,jE) = 28.0002_r8

    ! in (NH4)3H(SO4)2
    jE = jlvcite
    b_mtem(1,jA,jE) = -3.09464_r8
    b_mtem(2,jA,jE) = 14.9303_r8
    b_mtem(3,jA,jE) = -43.0454_r8
    b_mtem(4,jA,jE) = 72.6695_r8
    b_mtem(5,jA,jE) = -65.2140_r8
    b_mtem(6,jA,jE) = 23.4814_r8

    ! in NaHSO4
    jE = jnahso4
    b_mtem(1,jA,jE) = 1.22973_r8
    b_mtem(2,jA,jE) = 2.82702_r8
    b_mtem(3,jA,jE) = -17.5869_r8
    b_mtem(4,jA,jE) = 28.9564_r8
    b_mtem(5,jA,jE) = -23.5814_r8
    b_mtem(6,jA,jE) = 7.91153_r8

    ! in Na3H(SO4)2
    jE = jna3hso4
    b_mtem(1,jA,jE) = 1.64773_r8
    b_mtem(2,jA,jE) = 0.94188_r8
    b_mtem(3,jA,jE) = -19.1242_r8
    b_mtem(4,jA,jE) = 46.9887_r8
    b_mtem(5,jA,jE) = -50.9494_r8
    b_mtem(6,jA,jE) = 20.2169_r8


    !----------
    ! HCl in E
    jA = jhcl

    ! in (NH4)2SO4
    jE = jnh4so4
    b_mtem(1,jA,jE) = -2.93783_r8
    b_mtem(2,jA,jE) = 20.5546_r8
    b_mtem(3,jA,jE) = -75.8548_r8
    b_mtem(4,jA,jE) = 141.729_r8
    b_mtem(5,jA,jE) = -130.697_r8
    b_mtem(6,jA,jE) = 46.9905_r8

    ! in NH4NO3
    jE = jnh4no3
    b_mtem(1,jA,jE) = -1.69063_r8
    b_mtem(2,jA,jE) = -1.85303_r8
    b_mtem(3,jA,jE) = 29.0927_r8
    b_mtem(4,jA,jE) = -58.7401_r8
    b_mtem(5,jA,jE) = 44.999_r8
    b_mtem(6,jA,jE) = -11.9988_r8

    ! in NH4Cl (revised on 11/15/2003)
    jE = jnh4cl
    b_mtem(1,jA,jE) = -0.2073_r8
    b_mtem(2,jA,jE) = -0.4322_r8
    b_mtem(3,jA,jE) = 6.1271_r8
    b_mtem(4,jA,jE) = -12.3146_r8
    b_mtem(5,jA,jE) = 8.9919_r8
    b_mtem(6,jA,jE) = -2.3388_r8

    ! in NaCl
    jE = jnacl
    b_mtem(1,jA,jE) = 2.95913_r8
    b_mtem(2,jA,jE) = -7.92254_r8
    b_mtem(3,jA,jE) = 13.736_r8
    b_mtem(4,jA,jE) = -15.433_r8
    b_mtem(5,jA,jE) = 7.40386_r8
    b_mtem(6,jA,jE) = -0.918641_r8

    ! in NaNO3
    jE = jnano3
    b_mtem(1,jA,jE) = 0.893272_r8
    b_mtem(2,jA,jE) = 6.53768_r8
    b_mtem(3,jA,jE) = -32.3458_r8
    b_mtem(4,jA,jE) = 61.2834_r8
    b_mtem(5,jA,jE) = -56.4446_r8
    b_mtem(6,jA,jE) = 19.9202_r8

    ! in Na2SO4
    jE = jna2so4
    b_mtem(1,jA,jE) = 3.14484_r8
    b_mtem(2,jA,jE) = 0.077019_r8
    b_mtem(3,jA,jE) = -31.4199_r8
    b_mtem(4,jA,jE) = 80.5865_r8
    b_mtem(5,jA,jE) = -85.392_r8
    b_mtem(6,jA,jE) = 32.6644_r8

    ! in Ca(NO3)2
    jE = jcano3
    b_mtem(1,jA,jE) = 2.60432_r8
    b_mtem(2,jA,jE) = -0.55909_r8
    b_mtem(3,jA,jE) = -19.6671_r8
    b_mtem(4,jA,jE) = 53.3446_r8
    b_mtem(5,jA,jE) = -58.9076_r8
    b_mtem(6,jA,jE) = 22.9927_r8

    ! in CaCl2 (KM) revised on 3/13/2003 and again on 11/27/2003
    jE = jcacl2
    b_mtem(1,jA,jE) = 2.98036_r8
    b_mtem(2,jA,jE) = -8.55365_r8
    b_mtem(3,jA,jE) = 15.2108_r8
    b_mtem(4,jA,jE) = -15.9359_r8
    b_mtem(5,jA,jE) = 7.41772_r8
    b_mtem(6,jA,jE) = -1.32143_r8

    ! in HNO3 (added on 12/06/2004)
    jE = jhno3
    b_mtem(1,jA,jE) = 3.8533_r8
    b_mtem(2,jA,jE) = -16.9427_r8
    b_mtem(3,jA,jE) = 45.0056_r8
    b_mtem(4,jA,jE) = -69.6145_r8
    b_mtem(5,jA,jE) = 54.1491_r8
    b_mtem(6,jA,jE) = -16.6513_r8

    ! in HCl (added on 12/06/2004)
    jE = jhcl
    b_mtem(1,jA,jE) = 2.56665_r8
    b_mtem(2,jA,jE) = -7.13585_r8
    b_mtem(3,jA,jE) = 14.8103_r8
    b_mtem(4,jA,jE) = -21.8881_r8
    b_mtem(5,jA,jE) = 16.6808_r8
    b_mtem(6,jA,jE) = -5.22091_r8

    ! in H2SO4
    jE = jh2so4
    b_mtem(1,jA,jE) = 2.50179_r8
    b_mtem(2,jA,jE) = -6.69364_r8
    b_mtem(3,jA,jE) = 11.6551_r8
    b_mtem(4,jA,jE) = -13.6897_r8
    b_mtem(5,jA,jE) = 7.36796_r8
    b_mtem(6,jA,jE) = -1.33245_r8

    ! in NH4HSO4
    jE = jnh4hso4
    b_mtem(1,jA,jE) = 0.149955_r8
    b_mtem(2,jA,jE) = 11.8213_r8
    b_mtem(3,jA,jE) = -53.9164_r8
    b_mtem(4,jA,jE) = 101.574_r8
    b_mtem(5,jA,jE) = -91.4123_r8
    b_mtem(6,jA,jE) = 31.5487_r8

    ! in (NH4)3H(SO4)2
    jE = jlvcite
    b_mtem(1,jA,jE) = -2.36927_r8
    b_mtem(2,jA,jE) = 14.8359_r8
    b_mtem(3,jA,jE) = -44.3443_r8
    b_mtem(4,jA,jE) = 73.6229_r8
    b_mtem(5,jA,jE) = -65.3366_r8
    b_mtem(6,jA,jE) = 23.3250_r8

    ! in NaHSO4
    jE = jnahso4
    b_mtem(1,jA,jE) = 2.72993_r8
    b_mtem(2,jA,jE) = -0.23406_r8
    b_mtem(3,jA,jE) = -10.4103_r8
    b_mtem(4,jA,jE) = 13.1586_r8
    b_mtem(5,jA,jE) = -7.79925_r8
    b_mtem(6,jA,jE) = 2.30843_r8

    ! in Na3H(SO4)2
    jE = jna3hso4
    b_mtem(1,jA,jE) = 3.51258_r8
    b_mtem(2,jA,jE) = -3.95107_r8
    b_mtem(3,jA,jE) = -11.0175_r8
    b_mtem(4,jA,jE) = 38.8617_r8
    b_mtem(5,jA,jE) = -48.1575_r8
    b_mtem(6,jA,jE) = 20.4717_r8


    !----------
    ! 2H.SO4 in E
    jA = jh2so4

    ! in H2SO4
    jE = jh2so4
    b_mtem(1,jA,jE) = 0.76734_r8
    b_mtem(2,jA,jE) = -1.12263_r8
    b_mtem(3,jA,jE) = -9.08728_r8
    b_mtem(4,jA,jE) = 30.3836_r8
    b_mtem(5,jA,jE) = -38.4133_r8
    b_mtem(6,jA,jE) = 17.0106_r8

    ! in NH4HSO4
    jE = jnh4hso4
    b_mtem(1,jA,jE) = -2.03879_r8
    b_mtem(2,jA,jE) = 15.7033_r8
    b_mtem(3,jA,jE) = -58.7363_r8
    b_mtem(4,jA,jE) = 109.242_r8
    b_mtem(5,jA,jE) = -102.237_r8
    b_mtem(6,jA,jE) = 37.5350_r8

    ! in (NH4)3H(SO4)2
    jE = jlvcite
    b_mtem(1,jA,jE) = -3.10228_r8
    b_mtem(2,jA,jE) = 16.6920_r8
    b_mtem(3,jA,jE) = -59.1522_r8
    b_mtem(4,jA,jE) = 113.487_r8
    b_mtem(5,jA,jE) = -110.890_r8
    b_mtem(6,jA,jE) = 42.4578_r8

    ! in (NH4)2SO4
    jE = jnh4so4
    b_mtem(1,jA,jE) = -3.43885_r8
    b_mtem(2,jA,jE) = 21.0372_r8
    b_mtem(3,jA,jE) = -84.7026_r8
    b_mtem(4,jA,jE) = 165.324_r8
    b_mtem(5,jA,jE) = -156.101_r8
    b_mtem(6,jA,jE) = 57.3101_r8

    ! in NaHSO4
    jE = jnahso4
    b_mtem(1,jA,jE) = 0.33164_r8
    b_mtem(2,jA,jE) = 6.55864_r8
    b_mtem(3,jA,jE) = -33.5876_r8
    b_mtem(4,jA,jE) = 65.1798_r8
    b_mtem(5,jA,jE) = -63.2046_r8
    b_mtem(6,jA,jE) = 24.1783_r8

    ! in Na3H(SO4)2
    jE = jna3hso4
    b_mtem(1,jA,jE) = 3.06830_r8
    b_mtem(2,jA,jE) = -3.18408_r8
    b_mtem(3,jA,jE) = -19.6332_r8
    b_mtem(4,jA,jE) = 61.3657_r8
    b_mtem(5,jA,jE) = -73.4438_r8
    b_mtem(6,jA,jE) = 31.2334_r8

    ! in Na2SO4
    jE = jna2so4
    b_mtem(1,jA,jE) = 2.58649_r8
    b_mtem(2,jA,jE) = 0.87921_r8
    b_mtem(3,jA,jE) = -39.3023_r8
    b_mtem(4,jA,jE) = 101.603_r8
    b_mtem(5,jA,jE) = -109.469_r8
    b_mtem(6,jA,jE) = 43.0188_r8

    ! in HNO3
    jE = jhno3
    b_mtem(1,jA,jE) = 1.54587_r8
    b_mtem(2,jA,jE) = -7.50976_r8
    b_mtem(3,jA,jE) = 12.8237_r8
    b_mtem(4,jA,jE) = -10.1452_r8
    b_mtem(5,jA,jE) = -0.541956_r8
    b_mtem(6,jA,jE) = 3.34536_r8

    ! in HCl
    jE = jhcl
    b_mtem(1,jA,jE) = 0.829757_r8
    b_mtem(2,jA,jE) = -4.11316_r8
    b_mtem(3,jA,jE) = 3.67111_r8
    b_mtem(4,jA,jE) = 3.6833_r8
    b_mtem(5,jA,jE) = -11.2711_r8
    b_mtem(6,jA,jE) = 6.71421_r8


    !----------
    ! H.HSO4 in E
    jA = jhhso4

    ! in H2SO4
    jE = jh2so4
    b_mtem(1,jA,jE) = 2.63953_r8
    b_mtem(2,jA,jE) = -6.01532_r8
    b_mtem(3,jA,jE) = 10.0204_r8
    b_mtem(4,jA,jE) = -12.4840_r8
    b_mtem(5,jA,jE) = 7.78853_r8
    b_mtem(6,jA,jE) = -2.12638_r8

    ! in NH4HSO4
    jE = jnh4hso4
    b_mtem(1,jA,jE) = -0.77412_r8
    b_mtem(2,jA,jE) = 14.1656_r8
    b_mtem(3,jA,jE) = -53.4087_r8
    b_mtem(4,jA,jE) = 93.2013_r8
    b_mtem(5,jA,jE) = -80.5723_r8
    b_mtem(6,jA,jE) = 27.1577_r8

    ! in (NH4)3H(SO4)2
    jE = jlvcite
    b_mtem(1,jA,jE) = -2.98882_r8
    b_mtem(2,jA,jE) = 14.4436_r8
    b_mtem(3,jA,jE) = -40.1774_r8
    b_mtem(4,jA,jE) = 67.5937_r8
    b_mtem(5,jA,jE) = -61.5040_r8
    b_mtem(6,jA,jE) = 22.3695_r8

    ! in (NH4)2SO4
    jE = jnh4so4
    b_mtem(1,jA,jE) = -1.15502_r8
    b_mtem(2,jA,jE) = 8.12309_r8
    b_mtem(3,jA,jE) = -38.4726_r8
    b_mtem(4,jA,jE) = 80.8861_r8
    b_mtem(5,jA,jE) = -80.1644_r8
    b_mtem(6,jA,jE) = 30.4717_r8

    ! in NaHSO4
    jE = jnahso4
    b_mtem(1,jA,jE) = 1.99641_r8
    b_mtem(2,jA,jE) = -2.96061_r8
    b_mtem(3,jA,jE) = 5.54778_r8
    b_mtem(4,jA,jE) = -14.5488_r8
    b_mtem(5,jA,jE) = 14.8492_r8
    b_mtem(6,jA,jE) = -5.1389_r8

    ! in Na3H(SO4)2
    jE = jna3hso4
    b_mtem(1,jA,jE) = 2.23816_r8
    b_mtem(2,jA,jE) = -3.20847_r8
    b_mtem(3,jA,jE) = -4.82853_r8
    b_mtem(4,jA,jE) = 20.9192_r8
    b_mtem(5,jA,jE) = -27.2819_r8
    b_mtem(6,jA,jE) = 11.8655_r8

    ! in Na2SO4
    jE = jna2so4
    b_mtem(1,jA,jE) = 2.56907_r8
    b_mtem(2,jA,jE) = 1.13444_r8
    b_mtem(3,jA,jE) = -34.6853_r8
    b_mtem(4,jA,jE) = 87.9775_r8
    b_mtem(5,jA,jE) = -93.2330_r8
    b_mtem(6,jA,jE) = 35.9260_r8

    ! in HNO3
    jE = jhno3
    b_mtem(1,jA,jE) = 2.00024_r8
    b_mtem(2,jA,jE) = -4.80868_r8
    b_mtem(3,jA,jE) = 8.29222_r8
    b_mtem(4,jA,jE) = -11.0849_r8
    b_mtem(5,jA,jE) = 7.51262_r8
    b_mtem(6,jA,jE) = -2.07654_r8

    ! in HCl
    jE = jhcl
    b_mtem(1,jA,jE) = 2.8009_r8
    b_mtem(2,jA,jE) = -6.98416_r8
    b_mtem(3,jA,jE) = 14.3146_r8
    b_mtem(4,jA,jE) = -22.0068_r8
    b_mtem(5,jA,jE) = 17.5557_r8
    b_mtem(6,jA,jE) = -5.84917_r8


    !----------
    ! NH4HSO4 in E
    jA = jnh4hso4

    ! in H2SO4
    jE = jh2so4
    b_mtem(1,jA,jE) = 0.169160_r8
    b_mtem(2,jA,jE) = 2.15094_r8
    b_mtem(3,jA,jE) = -9.62904_r8
    b_mtem(4,jA,jE) = 18.2631_r8
    b_mtem(5,jA,jE) = -17.3333_r8
    b_mtem(6,jA,jE) = 6.19835_r8

    ! in NH4HSO4
    jE = jnh4hso4
    b_mtem(1,jA,jE) = -2.34457_r8
    b_mtem(2,jA,jE) = 12.8035_r8
    b_mtem(3,jA,jE) = -35.2513_r8
    b_mtem(4,jA,jE) = 53.6153_r8
    b_mtem(5,jA,jE) = -42.7655_r8
    b_mtem(6,jA,jE) = 13.7129_r8

    ! in (NH4)3H(SO4)2
    jE = jlvcite
    b_mtem(1,jA,jE) = -2.56109_r8
    b_mtem(2,jA,jE) = 11.1414_r8
    b_mtem(3,jA,jE) = -30.2361_r8
    b_mtem(4,jA,jE) = 50.0320_r8
    b_mtem(5,jA,jE) = -44.1586_r8
    b_mtem(6,jA,jE) = 15.5393_r8

    ! in (NH4)2SO4
    jE = jnh4so4
    b_mtem(1,jA,jE) = -0.97315_r8
    b_mtem(2,jA,jE) = 7.06295_r8
    b_mtem(3,jA,jE) = -29.3032_r8
    b_mtem(4,jA,jE) = 57.6101_r8
    b_mtem(5,jA,jE) = -54.9020_r8
    b_mtem(6,jA,jE) = 20.2222_r8

    ! in NaHSO4
    jE = jnahso4
    b_mtem(1,jA,jE) = -0.44450_r8
    b_mtem(2,jA,jE) = 3.33451_r8
    b_mtem(3,jA,jE) = -15.2791_r8
    b_mtem(4,jA,jE) = 30.1413_r8
    b_mtem(5,jA,jE) = -26.7710_r8
    b_mtem(6,jA,jE) = 8.78462_r8

    ! in Na3H(SO4)2
    jE = jna3hso4
    b_mtem(1,jA,jE) = -0.99780_r8
    b_mtem(2,jA,jE) = 4.69200_r8
    b_mtem(3,jA,jE) = -16.1219_r8
    b_mtem(4,jA,jE) = 29.3100_r8
    b_mtem(5,jA,jE) = -26.3383_r8
    b_mtem(6,jA,jE) = 9.20695_r8

    ! in Na2SO4
    jE = jna2so4
    b_mtem(1,jA,jE) = -0.52694_r8
    b_mtem(2,jA,jE) = 7.02684_r8
    b_mtem(3,jA,jE) = -33.7508_r8
    b_mtem(4,jA,jE) = 70.0565_r8
    b_mtem(5,jA,jE) = -68.3226_r8
    b_mtem(6,jA,jE) = 25.2692_r8

    ! in HNO3
    jE = jhno3
    b_mtem(1,jA,jE) = 0.572926_r8
    b_mtem(2,jA,jE) = -2.04791_r8
    b_mtem(3,jA,jE) = 2.1134_r8
    b_mtem(4,jA,jE) = 0.246654_r8
    b_mtem(5,jA,jE) = -3.06019_r8
    b_mtem(6,jA,jE) = 1.98126_r8

    ! in HCl
    jE = jhcl
    b_mtem(1,jA,jE) = 0.56514_r8
    b_mtem(2,jA,jE) = 0.22287_r8
    b_mtem(3,jA,jE) = -2.76973_r8
    b_mtem(4,jA,jE) = 4.54444_r8
    b_mtem(5,jA,jE) = -3.86549_r8
    b_mtem(6,jA,jE) = 1.13441_r8


    !----------
    ! (NH4)3H(SO4)2 in E
    jA = jlvcite

    ! in H2SO4
    jE = jh2so4
    b_mtem(1,jA,jE) = -1.44811_r8
    b_mtem(2,jA,jE) = 6.71815_r8
    b_mtem(3,jA,jE) = -25.0141_r8
    b_mtem(4,jA,jE) = 50.1109_r8
    b_mtem(5,jA,jE) = -50.0561_r8
    b_mtem(6,jA,jE) = 19.3370_r8

    ! in NH4HSO4
    jE = jnh4hso4
    b_mtem(1,jA,jE) = -3.41707_r8
    b_mtem(2,jA,jE) = 13.4496_r8
    b_mtem(3,jA,jE) = -34.8018_r8
    b_mtem(4,jA,jE) = 55.2987_r8
    b_mtem(5,jA,jE) = -48.1839_r8
    b_mtem(6,jA,jE) = 17.2444_r8

    ! in (NH4)3H(SO4)2
    jE = jlvcite
    b_mtem(1,jA,jE) = -2.54479_r8
    b_mtem(2,jA,jE) = 11.8501_r8
    b_mtem(3,jA,jE) = -39.7286_r8
    b_mtem(4,jA,jE) = 74.2479_r8
    b_mtem(5,jA,jE) = -70.4934_r8
    b_mtem(6,jA,jE) = 26.2836_r8

    ! in (NH4)2SO4
    jE = jnh4so4
    b_mtem(1,jA,jE) = -2.30561_r8
    b_mtem(2,jA,jE) = 14.5806_r8
    b_mtem(3,jA,jE) = -55.1238_r8
    b_mtem(4,jA,jE) = 103.451_r8
    b_mtem(5,jA,jE) = -95.2571_r8
    b_mtem(6,jA,jE) = 34.2218_r8

    ! in NaHSO4
    jE = jnahso4
    b_mtem(1,jA,jE) = -2.20809_r8
    b_mtem(2,jA,jE) = 13.6391_r8
    b_mtem(3,jA,jE) = -57.8246_r8
    b_mtem(4,jA,jE) = 117.907_r8
    b_mtem(5,jA,jE) = -112.154_r8
    b_mtem(6,jA,jE) = 40.3058_r8

    ! in Na3H(SO4)2
    jE = jna3hso4
    b_mtem(1,jA,jE) = -1.15099_r8
    b_mtem(2,jA,jE) = 6.32269_r8
    b_mtem(3,jA,jE) = -27.3860_r8
    b_mtem(4,jA,jE) = 55.4592_r8
    b_mtem(5,jA,jE) = -54.0100_r8
    b_mtem(6,jA,jE) = 20.3469_r8

    ! in Na2SO4
    jE = jna2so4
    b_mtem(1,jA,jE) = -1.15678_r8
    b_mtem(2,jA,jE) = 8.28718_r8
    b_mtem(3,jA,jE) = -37.3231_r8
    b_mtem(4,jA,jE) = 76.6124_r8
    b_mtem(5,jA,jE) = -74.9307_r8
    b_mtem(6,jA,jE) = 28.0559_r8

    ! in HNO3
    jE = jhno3
    b_mtem(1,jA,jE) = 0.01502_r8
    b_mtem(2,jA,jE) = -3.1197_r8
    b_mtem(3,jA,jE) = 3.61104_r8
    b_mtem(4,jA,jE) = 3.05196_r8
    b_mtem(5,jA,jE) = -9.98957_r8
    b_mtem(6,jA,jE) = 6.04155_r8

    ! in HCl
    jE = jhcl
    b_mtem(1,jA,jE) = -1.06477_r8
    b_mtem(2,jA,jE) = 3.38801_r8
    b_mtem(3,jA,jE) = -12.5784_r8
    b_mtem(4,jA,jE) = 25.2823_r8
    b_mtem(5,jA,jE) = -25.4611_r8
    b_mtem(6,jA,jE) = 10.0754_r8


    !----------
    ! NaHSO4 in E
    jA = jnahso4

    ! in H2SO4
    jE = jh2so4
    b_mtem(1,jA,jE) = 0.68259_r8
    b_mtem(2,jA,jE) = 0.71468_r8
    b_mtem(3,jA,jE) = -5.59003_r8
    b_mtem(4,jA,jE) = 11.0089_r8
    b_mtem(5,jA,jE) = -10.7983_r8
    b_mtem(6,jA,jE) = 3.82335_r8

    ! in NH4HSO4
    jE = jnh4hso4
    b_mtem(1,jA,jE) = -0.03956_r8
    b_mtem(2,jA,jE) = 4.52828_r8
    b_mtem(3,jA,jE) = -25.2557_r8
    b_mtem(4,jA,jE) = 54.4225_r8
    b_mtem(5,jA,jE) = -52.5105_r8
    b_mtem(6,jA,jE) = 18.6562_r8

    ! in (NH4)3H(SO4)2
    jE = jlvcite
    b_mtem(1,jA,jE) = -1.53503_r8
    b_mtem(2,jA,jE) = 8.27608_r8
    b_mtem(3,jA,jE) = -28.9539_r8
    b_mtem(4,jA,jE) = 55.2876_r8
    b_mtem(5,jA,jE) = -51.9563_r8
    b_mtem(6,jA,jE) = 18.6576_r8

    ! in (NH4)2SO4
    jE = jnh4so4
    b_mtem(1,jA,jE) = -0.38793_r8
    b_mtem(2,jA,jE) = 7.14680_r8
    b_mtem(3,jA,jE) = -38.7201_r8
    b_mtem(4,jA,jE) = 84.3965_r8
    b_mtem(5,jA,jE) = -84.7453_r8
    b_mtem(6,jA,jE) = 32.1283_r8

    ! in NaHSO4
    jE = jnahso4
    b_mtem(1,jA,jE) = -0.41982_r8
    b_mtem(2,jA,jE) = 4.26491_r8
    b_mtem(3,jA,jE) = -20.2351_r8
    b_mtem(4,jA,jE) = 42.6764_r8
    b_mtem(5,jA,jE) = -40.7503_r8
    b_mtem(6,jA,jE) = 14.2868_r8

    ! in Na3H(SO4)2
    jE = jna3hso4
    b_mtem(1,jA,jE) = -0.32912_r8
    b_mtem(2,jA,jE) = 1.80808_r8
    b_mtem(3,jA,jE) = -8.01286_r8
    b_mtem(4,jA,jE) = 15.5791_r8
    b_mtem(5,jA,jE) = -14.5494_r8
    b_mtem(6,jA,jE) = 5.27052_r8

    ! in Na2SO4
    jE = jna2so4
    b_mtem(1,jA,jE) = 0.10271_r8
    b_mtem(2,jA,jE) = 5.09559_r8
    b_mtem(3,jA,jE) = -30.3295_r8
    b_mtem(4,jA,jE) = 66.2975_r8
    b_mtem(5,jA,jE) = -66.3458_r8
    b_mtem(6,jA,jE) = 24.9443_r8

    ! in HNO3
    jE = jhno3
    b_mtem(1,jA,jE) = 0.608309_r8
    b_mtem(2,jA,jE) = -0.541905_r8
    b_mtem(3,jA,jE) = -2.52084_r8
    b_mtem(4,jA,jE) = 6.63297_r8
    b_mtem(5,jA,jE) = -7.24599_r8
    b_mtem(6,jA,jE) = 2.88811_r8

    ! in HCl
    jE = jhcl
    b_mtem(1,jA,jE) = 1.98399_r8
    b_mtem(2,jA,jE) = -4.51562_r8
    b_mtem(3,jA,jE) = 8.36059_r8
    b_mtem(4,jA,jE) = -12.4948_r8
    b_mtem(5,jA,jE) = 9.67514_r8
    b_mtem(6,jA,jE) = -3.18004_r8


    !----------
    ! Na3H(SO4)2 in E
    jA = jna3hso4

    ! in H2SO4
    jE = jh2so4
    b_mtem(1,jA,jE) = -0.83214_r8
    b_mtem(2,jA,jE) = 4.99572_r8
    b_mtem(3,jA,jE) = -20.1697_r8
    b_mtem(4,jA,jE) = 41.4066_r8
    b_mtem(5,jA,jE) = -42.2119_r8
    b_mtem(6,jA,jE) = 16.4855_r8

    ! in NH4HSO4
    jE = jnh4hso4
    b_mtem(1,jA,jE) = -0.65139_r8
    b_mtem(2,jA,jE) = 3.52300_r8
    b_mtem(3,jA,jE) = -22.8220_r8
    b_mtem(4,jA,jE) = 56.2956_r8
    b_mtem(5,jA,jE) = -59.9028_r8
    b_mtem(6,jA,jE) = 23.1844_r8

    ! in (NH4)3H(SO4)2
    jE = jlvcite
    b_mtem(1,jA,jE) = -1.31331_r8
    b_mtem(2,jA,jE) = 8.40835_r8
    b_mtem(3,jA,jE) = -38.1757_r8
    b_mtem(4,jA,jE) = 80.5312_r8
    b_mtem(5,jA,jE) = -79.8346_r8
    b_mtem(6,jA,jE) = 30.0219_r8

    ! in (NH4)2SO4
    jE = jnh4so4
    b_mtem(1,jA,jE) = -1.03054_r8
    b_mtem(2,jA,jE) = 8.08155_r8
    b_mtem(3,jA,jE) = -38.1046_r8
    b_mtem(4,jA,jE) = 78.7168_r8
    b_mtem(5,jA,jE) = -77.2263_r8
    b_mtem(6,jA,jE) = 29.1521_r8

    ! in NaHSO4
    jE = jnahso4
    b_mtem(1,jA,jE) = -1.90695_r8
    b_mtem(2,jA,jE) = 11.6241_r8
    b_mtem(3,jA,jE) = -50.3175_r8
    b_mtem(4,jA,jE) = 105.884_r8
    b_mtem(5,jA,jE) = -103.258_r8
    b_mtem(6,jA,jE) = 37.6588_r8

    ! in Na3H(SO4)2
    jE = jna3hso4
    b_mtem(1,jA,jE) = -0.34780_r8
    b_mtem(2,jA,jE) = 2.85363_r8
    b_mtem(3,jA,jE) = -17.6224_r8
    b_mtem(4,jA,jE) = 38.9220_r8
    b_mtem(5,jA,jE) = -39.8106_r8
    b_mtem(6,jA,jE) = 15.6055_r8

    ! in Na2SO4
    jE = jna2so4
    b_mtem(1,jA,jE) = -0.75230_r8
    b_mtem(2,jA,jE) = 10.0140_r8
    b_mtem(3,jA,jE) = -50.5677_r8
    b_mtem(4,jA,jE) = 106.941_r8
    b_mtem(5,jA,jE) = -105.534_r8
    b_mtem(6,jA,jE) = 39.5196_r8

    ! in HNO3
    jE = jhno3
    b_mtem(1,jA,jE) = 0.057456_r8
    b_mtem(2,jA,jE) = -1.31264_r8
    b_mtem(3,jA,jE) = -1.94662_r8
    b_mtem(4,jA,jE) = 10.7024_r8
    b_mtem(5,jA,jE) = -14.9946_r8
    b_mtem(6,jA,jE) = 7.12161_r8

    ! in HCl
    jE = jhcl
    b_mtem(1,jA,jE) = 0.637894_r8
    b_mtem(2,jA,jE) = -2.29719_r8
    b_mtem(3,jA,jE) = 0.765361_r8
    b_mtem(4,jA,jE) = 4.8748_r8
    b_mtem(5,jA,jE) = -9.25978_r8
    b_mtem(6,jA,jE) = 4.91773_r8
    !
    !
    !
    !----------------------------------------------------------
    ! Coefficients for %MDRH(T) = d1 + d2*T + d3*T^2 + d4*T^3    (T in Kelvin)
    ! valid Temperature Range: 240 - 320 K
    !----------------------------------------------------------
    !
    ! SULFATE-POOR SYSTEMS
    ! AC
    j_index = 1
    d_mdrh(j_index,1) = -58.00268351_r8
    d_mdrh(j_index,2) = 2.031077573_r8
    d_mdrh(j_index,3) = -0.008281218_r8
    d_mdrh(j_index,4) = 1.00447E-05_r8

    ! AN
    j_index = 2
    d_mdrh(j_index,1) = 1039.137773_r8
    d_mdrh(j_index,2) = -11.47847095_r8
    d_mdrh(j_index,3) = 0.047702786_r8
    d_mdrh(j_index,4) = -6.77675E-05_r8

    ! AS
    j_index = 3
    d_mdrh(j_index,1) = 115.8366357_r8
    d_mdrh(j_index,2) = 0.491881663_r8
    d_mdrh(j_index,3) = -0.00422807_r8
    d_mdrh(j_index,4) = 7.29274E-06_r8

    ! SC
    j_index = 4
    d_mdrh(j_index,1) = 253.2424151_r8
    d_mdrh(j_index,2) = -1.429957864_r8
    d_mdrh(j_index,3) = 0.003727554_r8
    d_mdrh(j_index,4) = -3.13037E-06_r8

    ! SN
    j_index = 5
    d_mdrh(j_index,1) = -372.4306506_r8
    d_mdrh(j_index,2) = 5.3955633_r8
    d_mdrh(j_index,3) = -0.019804438_r8
    d_mdrh(j_index,4) = 2.25662E-05_r8

    ! SS
    j_index = 6
    d_mdrh(j_index,1) = 286.1271416_r8
    d_mdrh(j_index,2) = -1.670787758_r8
    d_mdrh(j_index,3) = 0.004431373_r8
    d_mdrh(j_index,4) = -3.57757E-06_r8

    ! CC
    j_index = 7
    d_mdrh(j_index,1) = -1124.07059_r8
    d_mdrh(j_index,2) = 14.26364209_r8
    d_mdrh(j_index,3) = -0.054816822_r8
    d_mdrh(j_index,4) = 6.70107E-05_r8

    ! CN
    j_index = 8
    d_mdrh(j_index,1) = 1855.413934_r8
    d_mdrh(j_index,2) = -20.29219473_r8
    d_mdrh(j_index,3) = 0.07807482_r8
    d_mdrh(j_index,4) = -1.017887858e-4_r8

    ! AN + AC
    j_index = 9
    d_mdrh(j_index,1) = 1761.176886_r8
    d_mdrh(j_index,2) = -19.29811062_r8
    d_mdrh(j_index,3) = 0.075676987_r8
    d_mdrh(j_index,4) = -1.0116959e-4_r8

    ! AS + AC
    j_index = 10
    d_mdrh(j_index,1) = 122.1074303_r8
    d_mdrh(j_index,2) = 0.429692122_r8
    d_mdrh(j_index,3) = -0.003928277_r8
    d_mdrh(j_index,4) = 6.43275E-06_r8

    ! AS + AN
    j_index = 11
    d_mdrh(j_index,1) = 2424.634678_r8
    d_mdrh(j_index,2) = -26.54031307_r8
    d_mdrh(j_index,3) = 0.101625387_r8
    d_mdrh(j_index,4) = -1.31544547798e-4_r8

    ! AS + AN + AC
    j_index = 12
    d_mdrh(j_index,1) = 2912.082599_r8
    d_mdrh(j_index,2) = -31.8894185_r8
    d_mdrh(j_index,3) = 0.121185849_r8
    d_mdrh(j_index,4) = -1.556534623e-4_r8

    ! SC + AC
    j_index = 13
    d_mdrh(j_index,1) = 172.2596493_r8
    d_mdrh(j_index,2) = -0.511006195_r8
    d_mdrh(j_index,3) = 4.27244597e-4_r8
    d_mdrh(j_index,4) = 4.12797E-07_r8

    ! SN + AC
    j_index = 14
    d_mdrh(j_index,1) = 1596.184935_r8
    d_mdrh(j_index,2) = -16.37945565_r8
    d_mdrh(j_index,3) = 0.060281218_r8
    d_mdrh(j_index,4) = -7.6161E-05_r8

    ! SN + AN
    j_index = 15
    d_mdrh(j_index,1) = 1916.072988_r8
    d_mdrh(j_index,2) = -20.85594868_r8
    d_mdrh(j_index,3) = 0.081140141_r8
    d_mdrh(j_index,4) = -1.07954274796e-4_r8

    ! SN + AN + AC
    j_index = 16
    d_mdrh(j_index,1) = 1467.165935_r8
    d_mdrh(j_index,2) = -16.01166196_r8
    d_mdrh(j_index,3) = 0.063505582_r8
    d_mdrh(j_index,4) = -8.66722E-05_r8

    ! SN + SC
    j_index = 17
    d_mdrh(j_index,1) = 158.447059_r8
    d_mdrh(j_index,2) = -0.628167358_r8
    d_mdrh(j_index,3) = 0.002014448_r8
    d_mdrh(j_index,4) = -3.13037E-06_r8

    ! SN + SC + AC
    j_index = 18
    d_mdrh(j_index,1) = 1115.892468_r8
    d_mdrh(j_index,2) = -11.76936534_r8
    d_mdrh(j_index,3) = 0.045577399_r8
    d_mdrh(j_index,4) = -6.05779E-05_r8

    ! SS + AC
    j_index = 19
    d_mdrh(j_index,1) = 269.5432407_r8
    d_mdrh(j_index,2) = -1.319963885_r8
    d_mdrh(j_index,3) = 0.002592363_r8
    d_mdrh(j_index,4) = -1.44479E-06_r8

    ! SS + AN
    j_index = 20
    d_mdrh(j_index,1) = 2841.334784_r8
    d_mdrh(j_index,2) = -31.1889487_r8
    d_mdrh(j_index,3) = 0.118809274_r8
    d_mdrh(j_index,4) = -1.53007e-4_r8

    ! SS + AN + AC
    j_index = 21
    d_mdrh(j_index,1) = 2199.36914_r8
    d_mdrh(j_index,2) = -24.11926569_r8
    d_mdrh(j_index,3) = 0.092932361_r8
    d_mdrh(j_index,4) = -1.21774e-4_r8

    ! SS + AS
    j_index = 22
    d_mdrh(j_index,1) = 395.0051604_r8
    d_mdrh(j_index,2) = -2.521101657_r8
    d_mdrh(j_index,3) = 0.006139319_r8
    d_mdrh(j_index,4) = -4.43756E-06_r8

    ! SS + AS + AC
    j_index = 23
    d_mdrh(j_index,1) = 386.5150675_r8
    d_mdrh(j_index,2) = -2.4632138_r8
    d_mdrh(j_index,3) = 0.006139319_r8
    d_mdrh(j_index,4) = -4.98796E-06_r8

    ! SS + AS + AN
    j_index = 24
    d_mdrh(j_index,1) = 3101.538491_r8
    d_mdrh(j_index,2) = -34.19978105_r8
    d_mdrh(j_index,3) = 0.130118605_r8
    d_mdrh(j_index,4) = -1.66873e-4_r8

    ! SS + AS + AN + AC
    j_index = 25
    d_mdrh(j_index,1) = 2307.579403_r8
    d_mdrh(j_index,2) = -25.43136774_r8
    d_mdrh(j_index,3) = 0.098064728_r8
    d_mdrh(j_index,4) = -1.28301e-4_r8

    ! SS + SC
    j_index = 26
    d_mdrh(j_index,1) = 291.8309602_r8
    d_mdrh(j_index,2) = -1.828912974_r8
    d_mdrh(j_index,3) = 0.005053148_r8
    d_mdrh(j_index,4) = -4.57516E-06_r8

    ! SS + SC + AC
    j_index = 27
    d_mdrh(j_index,1) = 188.3914345_r8
    d_mdrh(j_index,2) = -0.631345031_r8
    d_mdrh(j_index,3) = 0.000622807_r8
    d_mdrh(j_index,4) = 4.47196E-07_r8

    ! SS + SN
    j_index = 28
    d_mdrh(j_index,1) = -167.1252839_r8
    d_mdrh(j_index,2) = 2.969828002_r8
    d_mdrh(j_index,3) = -0.010637255_r8
    d_mdrh(j_index,4) = 1.13175E-05_r8

    ! SS + SN + AC
    j_index = 29
    d_mdrh(j_index,1) = 1516.782768_r8
    d_mdrh(j_index,2) = -15.7922661_r8
    d_mdrh(j_index,3) = 0.058942209_r8
    d_mdrh(j_index,4) = -7.5301E-05_r8

    ! SS + SN + AN
    j_index = 30
    d_mdrh(j_index,1) = 1739.963163_r8
    d_mdrh(j_index,2) = -19.06576022_r8
    d_mdrh(j_index,3) = 0.07454963_r8
    d_mdrh(j_index,4) = -9.94302E-05_r8

    ! SS + SN + AN + AC
    j_index = 31
    d_mdrh(j_index,1) = 2152.104877_r8
    d_mdrh(j_index,2) = -23.74998008_r8
    d_mdrh(j_index,3) = 0.092256654_r8
    d_mdrh(j_index,4) = -1.21953e-4_r8

    ! SS + SN + SC
    j_index = 32
    d_mdrh(j_index,1) = 221.9976265_r8
    d_mdrh(j_index,2) = -1.311331272_r8
    d_mdrh(j_index,3) = 0.004406089_r8
    d_mdrh(j_index,4) = -5.88235E-06_r8

    ! SS + SN + SC + AC
    j_index = 33
    d_mdrh(j_index,1) = 1205.645615_r8
    d_mdrh(j_index,2) = -12.71353459_r8
    d_mdrh(j_index,3) = 0.048803922_r8
    d_mdrh(j_index,4) = -6.41899E-05_r8

    ! CC + AC
    j_index = 34
    d_mdrh(j_index,1) = 506.6737879_r8
    d_mdrh(j_index,2) = -3.723520818_r8
    d_mdrh(j_index,3) = 0.010814242_r8
    d_mdrh(j_index,4) = -1.21087E-05_r8

    ! CC + SC
    j_index = 35
    d_mdrh(j_index,1) = -1123.523841_r8
    d_mdrh(j_index,2) = 14.08345977_r8
    d_mdrh(j_index,3) = -0.053687823_r8
    d_mdrh(j_index,4) = 6.52219E-05_r8

    ! CC + SC + AC
    j_index = 36
    d_mdrh(j_index,1) = -1159.98607_r8
    d_mdrh(j_index,2) = 14.44309169_r8
    d_mdrh(j_index,3) = -0.054841073_r8
    d_mdrh(j_index,4) = 6.64259E-05_r8

    ! CN + AC
    j_index = 37
    d_mdrh(j_index,1) = 756.0747916_r8
    d_mdrh(j_index,2) = -8.546826257_r8
    d_mdrh(j_index,3) = 0.035798677_r8
    d_mdrh(j_index,4) = -5.06629E-05_r8

    ! CN + AN
    j_index = 38
    d_mdrh(j_index,1) = 338.668191_r8
    d_mdrh(j_index,2) = -2.971223403_r8
    d_mdrh(j_index,3) = 0.012294866_r8
    d_mdrh(j_index,4) = -1.87558E-05_r8

    ! CN + AN + AC
    j_index = 39
    d_mdrh(j_index,1) = -53.18033508_r8
    d_mdrh(j_index,2) = 0.663911748_r8
    d_mdrh(j_index,3) = 9.16326e-4_r8
    d_mdrh(j_index,4) = -6.70354E-06_r8

    ! CN + SC
    j_index = 40
    d_mdrh(j_index,1) = 3623.831129_r8
    d_mdrh(j_index,2) = -39.27226457_r8
    d_mdrh(j_index,3) = 0.144559515_r8
    d_mdrh(j_index,4) = -1.78159e-4_r8

    ! CN + SC + AC
    j_index = 41
    d_mdrh(j_index,1) = 3436.656743_r8
    d_mdrh(j_index,2) = -37.16192684_r8
    d_mdrh(j_index,3) = 0.136641377_r8
    d_mdrh(j_index,4) = -1.68262e-4_r8

    ! CN + SN
    j_index = 42
    d_mdrh(j_index,1) = 768.608476_r8
    d_mdrh(j_index,2) = -8.051517149_r8
    d_mdrh(j_index,3) = 0.032342332_r8
    d_mdrh(j_index,4) = -4.52224E-05_r8

    ! CN + SN + AC
    j_index = 43
    d_mdrh(j_index,1) = 33.58027951_r8
    d_mdrh(j_index,2) = -0.308772182_r8
    d_mdrh(j_index,3) = 0.004713639_r8
    d_mdrh(j_index,4) = -1.19658E-05_r8

    ! CN + SN + AN
    j_index = 44
    d_mdrh(j_index,1) = 57.80183041_r8
    d_mdrh(j_index,2) = 0.215264604_r8
    d_mdrh(j_index,3) = 4.11406e-4_r8
    d_mdrh(j_index,4) = -4.30702E-06_r8

    ! CN + SN + AN + AC
    j_index = 45
    d_mdrh(j_index,1) = -234.368984_r8
    d_mdrh(j_index,2) = 2.721045204_r8
    d_mdrh(j_index,3) = -0.006688341_r8
    d_mdrh(j_index,4) = 2.31729E-06_r8

    ! CN + SN + SC
    j_index = 46
    d_mdrh(j_index,1) = 3879.080557_r8
    d_mdrh(j_index,2) = -42.13562874_r8
    d_mdrh(j_index,3) = 0.155235005_r8
    d_mdrh(j_index,4) = -1.91387e-4_r8

    ! CN + SN + SC + AC
    j_index = 47
    d_mdrh(j_index,1) = 3600.576985_r8
    d_mdrh(j_index,2) = -39.0283489_r8
    d_mdrh(j_index,3) = 0.143710316_r8
    d_mdrh(j_index,4) = -1.77167e-4_r8

    ! CN + CC
    j_index = 48
    d_mdrh(j_index,1) = -1009.729826_r8
    d_mdrh(j_index,2) = 12.9145339_r8
    d_mdrh(j_index,3) = -0.049811146_r8
    d_mdrh(j_index,4) = 6.09563E-05_r8

    ! CN + CC + AC
    j_index = 49
    d_mdrh(j_index,1) = -577.0919514_r8
    d_mdrh(j_index,2) = 8.020324227_r8
    d_mdrh(j_index,3) = -0.031469556_r8
    d_mdrh(j_index,4) = 3.82181E-05_r8

    ! CN + CC + SC
    j_index = 50
    d_mdrh(j_index,1) = -728.9983499_r8
    d_mdrh(j_index,2) = 9.849458215_r8
    d_mdrh(j_index,3) = -0.03879257_r8
    d_mdrh(j_index,4) = 4.78844E-05_r8

    ! CN + CC + SC + AC
    j_index = 51
    d_mdrh(j_index,1) = -803.7026845_r8
    d_mdrh(j_index,2) = 10.61881494_r8
    d_mdrh(j_index,3) = -0.041402993_r8
    d_mdrh(j_index,4) = 5.08084E-05_r8

    !
    ! SULFATE-RICH SYSTEMS
    ! AB
    j_index = 52
    d_mdrh(j_index,1) = -493.6190458_r8
    d_mdrh(j_index,2) = 6.747053851_r8
    d_mdrh(j_index,3) = -0.026955267_r8
    d_mdrh(j_index,4) = 3.45118E-05_r8

    ! LV
    j_index = 53
    d_mdrh(j_index,1) = 53.37874093_r8
    d_mdrh(j_index,2) = 1.01368249_r8
    d_mdrh(j_index,3) = -0.005887513_r8
    d_mdrh(j_index,4) = 8.94393E-06_r8

    ! SB
    j_index = 54
    d_mdrh(j_index,1) = 206.619047_r8
    d_mdrh(j_index,2) = -1.342735684_r8
    d_mdrh(j_index,3) = 0.003197691_r8
    d_mdrh(j_index,4) = -1.93603E-06_r8

    ! AB + LV
    j_index = 55
    d_mdrh(j_index,1) = -493.6190458_r8
    d_mdrh(j_index,2) = 6.747053851_r8
    d_mdrh(j_index,3) = -0.026955267_r8
    d_mdrh(j_index,4) = 3.45118E-05_r8

    ! AS + LV
    j_index = 56
    d_mdrh(j_index,1) = 53.37874093_r8
    d_mdrh(j_index,2) = 1.01368249_r8
    d_mdrh(j_index,3) = -0.005887513_r8
    d_mdrh(j_index,4) = 8.94393E-06_r8

    ! SS + SB
    j_index = 57
    d_mdrh(j_index,1) = 206.619047_r8
    d_mdrh(j_index,2) = -1.342735684_r8
    d_mdrh(j_index,3) = 0.003197691_r8
    d_mdrh(j_index,4) = -1.93603E-06_r8

    ! SS + LV
    j_index = 58
    d_mdrh(j_index,1) = 41.7619047_r8
    d_mdrh(j_index,2) = 1.303872053_r8
    d_mdrh(j_index,3) = -0.007647908_r8
    d_mdrh(j_index,4) = 1.17845E-05_r8

    ! SS + AS + LV
    j_index = 59
    d_mdrh(j_index,1) = 41.7619047_r8
    d_mdrh(j_index,2) = 1.303872053_r8
    d_mdrh(j_index,3) = -0.007647908_r8
    d_mdrh(j_index,4) = 1.17845E-05_r8

    ! SS + AB
    j_index = 60
    d_mdrh(j_index,1) = -369.7142842_r8
    d_mdrh(j_index,2) = 5.512878771_r8
    d_mdrh(j_index,3) = -0.02301948_r8
    d_mdrh(j_index,4) = 3.0303E-05_r8

    ! SS + LV + AB
    j_index = 61
    d_mdrh(j_index,1) = -369.7142842_r8
    d_mdrh(j_index,2) = 5.512878771_r8
    d_mdrh(j_index,3) = -0.02301948_r8
    d_mdrh(j_index,4) = 3.0303E-05_r8

    ! SB + AB
    j_index = 62
    d_mdrh(j_index,1) = -162.8095232_r8
    d_mdrh(j_index,2) = 2.399326592_r8
    d_mdrh(j_index,3) = -0.009336219_r8
    d_mdrh(j_index,4) = 1.17845E-05_r8

    ! SS + SB + AB
    j_index = 63
    d_mdrh(j_index,1) = -735.4285689_r8
    d_mdrh(j_index,2) = 8.885521857_r8
    d_mdrh(j_index,3) = -0.033488456_r8
    d_mdrh(j_index,4) = 4.12458E-05_r8


    !   endif ! first

    return
  end subroutine load_mosaic_parameters



end module module_mosaic_init_aerpar
