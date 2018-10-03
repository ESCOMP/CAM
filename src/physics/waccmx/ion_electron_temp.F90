module ion_electron_temp

!---------------------------------------------------------------------------------
! Purpose:
!
! Module to compute the ion/electron temperature and dry static heating
!
! Authors: Joe McInerney/Hanli Liu/Art Richmond
!
!---------------------------------------------------------------------------------
  use shr_kind_mod,   only : r8 => shr_kind_r8            ! Real kind to declare variables
  use ppgrid,         only : pcols, pver, pverp           ! Dimensions and chunk bounds
  use cam_history,    only : outfld, hist_fld_active, write_inithist ! Routine to output fields to history files
  use cam_history,    only : horiz_only, addfld, add_default ! Routines and variables for adding fields to history output
  use physics_types,  only : physics_state, &             ! Structures containing physics state variables
                             physics_ptend, &             ! Structures containing physics tendency variables
                             physics_ptend_init           ! Routine to initialize physics tendency variables
  use physics_buffer, only : pbuf_add_field, &            ! 
                             pbuf_get_index,dtype_r8, &   !
                             physics_buffer_desc, &       !
                             pbuf_get_field, &            ! Needed to access physics buffer
                             pbuf_set_field
  use mo_jeuv,        only : nIonRates                    ! Number of ionization rates in mo_photo
  use shr_const_mod,  only : kboltz => shr_const_boltz, &
                             pi => shr_const_pi           ! Boltzmann constant and pi
  use chem_mods,      only : adv_mass                     ! Array holding mass values for short lived species
  use cam_abortutils, only : endrun
  use mo_chem_utls,   only : get_spc_ndx                  ! Routine to get index of adv_mass array for short lived species
  use constituents,   only : cnst_get_ind, cnst_mw        ! Routines to get molecular weights for constituents
  use solar_parms_data, only : f107=>solar_parms_f107     ! 10.7 cm solar flux
  use steady_state_tei, only : steady_state_tei_init, steady_state_tei_tend
  use perf_mod,         only : t_startf, t_stopf          ! timing utils
  use spmd_utils,       only : masterproc
  use cam_logfile,      only : iulog ! Output unit
  use ionos_state_mod,  only : ionos_state

  implicit none

  save
  
  private   ! Make default type private to the module

  !------------------------
  ! PUBLIC: interfaces 
  !------------------------
  public :: ion_electron_temp_init     ! Initialization
  public :: ion_electron_temp_register ! Registration of ionosphere variables in pbuf physics buffer
  public :: ion_electron_temp_inidat   ! Get fields from initial condition file into physics buffer
  public :: ion_electron_temp_tend     ! Calculate tendencies for extended model ionosphere
  public :: ion_electron_temp_readnl

  !------------------------------------------------------------------------
  ! PRIVATE: Rest of the data and interfaces are private to this module
  !------------------------------------------------------------------------   
  real(r8), parameter :: kboltz_ev = 8.617E-5_r8 ! Boltzmann constant (eV/K)
  real(r8), parameter :: temax = 7.0E3_r8        ! maximum electron temperature (K)
  real(r8), parameter :: dayOPFlux = 2.0E8_r8    ! Daytime O+ flux at upper boundary (
  real(r8), parameter :: nightOPFlux = -2.0E8_r8 ! Nighttime O+ flux at upper boundary (

  real(r8), parameter :: rads2Degs   = 180._r8/pi ! radians to degrees


! private data
  real(r8) :: rMassOp ! O+ molecular weight kg/kmol

  logical :: steady_state_ion_elec_temp = .true.

  integer :: index_te=-1, index_ti=-1  ! Indices to find ion and electron temperature in pbuf

contains

!==============================================================================

  subroutine ion_electron_temp_readnl(nlfile)

    use namelist_utils,  only: find_group_name
    use units,           only: getunit, freeunit
    use spmd_utils, only: mpicom, masterprocid, mpicom, mpi_logical

    character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

    ! Local variables
    integer :: unitn, ierr
    character(len=*), parameter :: subname = 'ion_electron_temp_readnl'

    namelist /ion_electron_temp_nl/ steady_state_ion_elec_temp

    if (masterproc) then
       unitn = getunit()
       open( unitn, file=trim(nlfile), status='old' )
       call find_group_name(unitn, 'ion_electron_temp_nl', status=ierr)
       if (ierr == 0) then
          read(unitn, ion_electron_temp_nl, iostat=ierr)
          if (ierr /= 0) then
             call endrun(subname // ':: ERROR reading namelist')
          end if
       end if
       close(unitn)
       call freeunit(unitn)

    end if

    call mpi_bcast (steady_state_ion_elec_temp, 1, mpi_logical, masterprocid, mpicom, ierr)

    if (masterproc) then
       write(iulog,*) subname//': steady_state_ion_elec_temp = ',steady_state_ion_elec_temp
    endif

  end subroutine ion_electron_temp_readnl

!==============================================================================
!==============================================================================

  subroutine ion_electron_temp_init(pbuf2d)
  
!-----------------------------------------------------------------------
! Time independent initialization for ionosphere simulation.
!-----------------------------------------------------------------------

    use phys_control,     only : phys_getopts !Method used to get flag for waccmx ionosphere output variables     

    type(physics_buffer_desc), pointer :: pbuf2d(:,:)

    logical :: history_waccmx
    integer :: indxOp,sIndxOp                           ! state%q or pbuf index for O+ mixing ratio

    if (steady_state_ion_elec_temp) then
       call steady_state_tei_init(pbuf2d)
    end if
    
    call phys_getopts(history_waccmx_out=history_waccmx)

    !-------------------------------------------------------------------------------
    !  Add history variables for ionosphere 
    !-------------------------------------------------------------------------------
    call addfld ('QIonElec' ,(/ 'lev' /), 'I', 'K/s', 'Electron Ion Thermal Heating Rate')
    call addfld ('TElec&IC'     ,(/ 'lev' /), 'I', 'K',      'Electron Temperature')
    call addfld ('TIon&IC'      ,(/ 'lev' /), 'I', 'K',      'Ion Temperature')
    call addfld ('TElec'        ,(/ 'lev' /), 'I', 'K',      'Electron Temperature')
    call addfld ('TIon'         ,(/ 'lev' /), 'I', 'K',      'Ion Temperature')
    call addfld ('ElecColDens'  ,horiz_only , 'I', 'TECU',   'Electron Column Density')
    if (.not.steady_state_ion_elec_temp) then
       call addfld ('QIN'          ,(/ 'lev' /), 'I', 'J/kg/s', 'Ion-neutral Heating')
       call addfld ('QEN'          ,(/ 'lev' /), 'I', ' ',      'Electron-neutral Heating')
       call addfld ('QEI'          ,(/ 'lev' /), 'I', ' ',      'Electron-ion Heating')
       call addfld ('LOSS_g3'      ,(/ 'lev' /), 'I', ' ',      'Loss Term g3')
       call addfld ('LOSS_EI'      ,(/ 'lev' /), 'I', ' ',      'Loss Term EI')
       call addfld ('LOSS_IN'      ,(/ 'lev' /), 'I', ' ',      'Loss Term IN')
       call addfld ('SOURCER'      ,(/ 'lev' /), 'I', ' ',      'SOURCER')
       call addfld ('SOURCEEff'    ,(/ 'lev' /), 'I', ' ',      'SOURCEEff')
       call addfld ('AURIPRATESUM' ,(/ 'lev' /), 'I', ' ',      'Auroral ionization')
       call addfld ('OpI'          ,(/ 'lev' /), 'I', ' ',      'O+ Ionosphere')
       call addfld ('eI'           ,(/ 'lev' /), 'I', ' ',      'e Ionosphere')
    end if

    call add_default ('TElec&IC'      , 0, ' ')
    call add_default ('TIon&IC'       , 0, ' ')

    !-------------------------------------------------------------------------------
    !  Set default values for ionosphere history variables
    !-------------------------------------------------------------------------------
    if (history_waccmx) then
       call add_default ('TElec'         , 1, ' ')
       call add_default ('TIon'          , 1, ' ')
       if (.not.steady_state_ion_elec_temp) then
          call add_default ('QIN'           , 1, ' ')
          call add_default ('QEN'           , 1, ' ')
          call add_default ('QEI'           , 1, ' ')
          call add_default ('SOURCER'  , 1, ' ')
          call add_default ('SOURCEEff'  , 1, ' ')
          call add_default ('AURIPRATESUM'  , 1, ' ')
       end if
    end if

    if (.not.steady_state_ion_elec_temp) then
       call cnst_get_ind( 'Op',  indxOp, abort=.false. )
       if (indxOp > 0) then
          rMassOp = cnst_mw(indxOP)
       else
          sIndxOp = get_spc_ndx( 'Op' )
          if (sIndxOp > 0) then
             rMassOp  = adv_mass(sIndxOp)
          else
             call endrun('update_teti: Cannot find short-lived index for Op in update_teti')         
          endif
       endif
    endif

  end subroutine ion_electron_temp_init

!==============================================================================     

  subroutine ion_electron_temp_register

    !-----------------------------------------------------------------------
    ! Register ionosphere variables with physics buffer:
    !
    ! Ion production rates pcols,pver,nIonRates,
    !   so firstdim = 1 middledim = pver lastdim = nIonRates.
    ! 
    ! pcols dimension and lchnk assumed here
    !
    !-----------------------------------------------------------------------
  
    !------------------------------------------------------------------------------
    ! Electron temperature in physics buffer (global so can write to history files) 
    !------------------------------------------------------------------------------
    call pbuf_add_field('TElec','global',dtype_r8,(/pcols,pver/), index_te)
    
    !--------------------------------------------------------------------------
    ! Ion temperature in physics buffer (global so can write to history files)
    !--------------------------------------------------------------------------
    call pbuf_add_field('TIon', 'global',dtype_r8,(/pcols,pver/), index_ti)

  end subroutine ion_electron_temp_register

!==============================================================================

  subroutine ion_electron_temp_inidat(ncid_ini, pbuf2d)

    !-----------------------------------------------------------------------
    ! Grab fields from initial condition file and put in physics buffer
    !-----------------------------------------------------------------------

    use pio,              only : file_desc_t
    use cam_grid_support, only : cam_grid_check, cam_grid_id
    use cam_grid_support, only : cam_grid_get_dim_names
    use cam_abortutils,   only : endrun
    use physics_buffer,   only : pbuf_set_field
    use ncdio_atm,        only : infld
    use ppgrid,           only : pcols, pver, begchunk, endchunk

    type(file_desc_t), intent(inout)   :: ncid_ini    ! Initial condition file id
    type(physics_buffer_desc), pointer :: pbuf2d(:,:) ! Physics buffer

    integer          :: grid_id
    character(len=4) :: dim1name, dim2name
    logical          :: found
    real(r8),pointer :: tE(:,:,:)   ! Electron temperature pointer
    real(r8),pointer :: tI(:,:,:)   ! Ion temperature pointer
    integer          :: ierr
    character(len=*), parameter :: subname='ION_ELECTRON_TEMP_INIDAT'
 
    found = .false.

    grid_id = cam_grid_id('physgrid')
    if (.not. cam_grid_check(grid_id)) then
      call endrun(trim(subname)//': Internal error, no "physgrid" grid')
    end if
    call cam_grid_get_dim_names(grid_id, dim1name, dim2name)

    if (index_te>0) then
       !---------------------------------------------------------------------------------
       ! Electron temperature in to physics buffer.  If not found use neutral temperature 
       !---------------------------------------------------------------------------------
       allocate(tE(pcols,pver,begchunk:endchunk))
       call infld( 'TElec',ncid_ini,dim1name, 'lev', dim2name, 1, pcols, 1, pver, begchunk, endchunk, &
            tE, found, gridname='physgrid')

       if (.not.found) then
          if (masterproc) write(iulog,*) 'ion_electron_temp_inidat: Could not find electron temperature in ic file. ' &
                                      // 'Using neutral temperature'
          call infld( 'T',ncid_ini,dim1name, 'lev', dim2name, 1, pcols, 1, pver, begchunk, endchunk, &
               tE, found, gridname='physgrid')
       endif

       call pbuf_set_field(pbuf2d, index_te, tE)

       deallocate(tE)
    endif

    if (index_ti>0) then
       !----------------------------------------------------------------------------
       ! Ion temperature in to physics buffer.  If not found use neutral temperature
       !----------------------------------------------------------------------------
       allocate(tI(pcols,pver,begchunk:endchunk))
       call infld( 'TIon',ncid_ini,dim1name, 'lev', dim2name, 1, pcols, 1, pver, begchunk, endchunk, &
            tI, found, gridname='physgrid')

       if (.not.found) then
          if (masterproc) write(iulog,*) 'ion_electron_temp_inidat: Could not find ion temperature in ic file. ' &
                                      // 'Using neutral temperature'
          call infld( 'T',ncid_ini,dim1name, 'lev', dim2name, 1, pcols, 1, pver, begchunk, endchunk, &
               tI, found, gridname='physgrid')
       endif

       call pbuf_set_field(pbuf2d, index_ti, tI)

       deallocate(tI)
    endif

  end subroutine ion_electron_temp_inidat

!==============================================================================

  subroutine ion_electron_temp_tend(state, ptend, pbuf, ztodt)

    use physconst,        only : cpairv
    !-------------------------------------------------------------------------------------
    ! Calculate dry static energy and O+ tendency for extended ionosphere simulation 
    !-------------------------------------------------------------------------------------

!------------------------------Arguments--------------------------------

    use physics_types,       only : physics_ptend_sum
    
    type(physics_state), intent(in)    :: state               ! physics state structure
    type(physics_ptend), intent(inout) :: ptend               ! parameterization tendency structure
    type(physics_buffer_desc),pointer  :: pbuf(:)             ! physics buffer

    real(r8),            intent(in)    :: ztodt               ! Physics time step

!---------------------------Local storage-------------------------------

    type(physics_ptend)                :: ptend_loc           ! Local parameterization tendencies
    type(ionos_state)                  :: istate              ! ionosphere state structure

    integer :: lchnk      ! Chunk number 
    integer :: ncol       ! Number of columns in chunk 
                
    integer :: teTiBot                            ! bottom of ionosphere calculations

    real(r8), dimension(:,:), pointer   :: tE           ! Pointer to electron temperature in pbuf (K) 
    real(r8), dimension(:,:), pointer   :: tI           ! Pointer to ion temperature in pbuf (K) 
 
    logical :: ls         
    real(r8) :: dse_tend(pcols,pver) ! dry static energy tendency
    real(r8) :: qionelec(pcols,pver) ! diagnostic heating rate (neutrals)

    call t_startf ('ion_electron_temp_tend')

    !----------------------------------------------------------------
    !  Get number of this chunk
    !----------------------------------------------------------------
    lchnk = state%lchnk
    ncol = state%ncol

    ls = .TRUE.   
    call physics_ptend_init(ptend_loc, state%psetcols, 'ionosphere', ls=ls)

    !-------------------------------------------------------------------------------------------------------------------
    !  Get electron and ion temperatures from physics buffer. 
    !-------------------------------------------------------------------------------------------------------------------
    call pbuf_get_field(pbuf, index_te, tE)
    call pbuf_get_field(pbuf, index_ti, tI)

    !------------------------------------------------------------
    !  Initialize data needed in the ionosphere calculations
    !------------------------------------------------------------
    call update_istate(state, pbuf, istate, teTiBot)

    if (steady_state_ion_elec_temp) then
       !-----------------------------------------------------------------
       ! steady-state solution
       !-----------------------------------------------------------------
       dse_tend(:ncol,:) = ptend%s(:ncol,:)
       call steady_state_tei_tend(state, istate, dse_tend, pbuf)
       ptend_loc%s(:ncol,:) = dse_tend(:ncol,:)
    else
       !-----------------------------------------------------------------
       !  Get electron temperature and update dry static energy tendency
       !-----------------------------------------------------------------
       call update_teti(state, ptend%s, ptend_loc%s, ztodt, istate, tE, tI, teTiBot)
    end if

    !--------------------------------------------------------------
    !  Make Te and Ti fields available for output to history files
    !--------------------------------------------------------------
    call outfld ('TElec'   , tE, pcols, lchnk)
    call outfld ('TIon'    , tI, pcols, lchnk)
    if (write_inithist()) then
       call outfld ('TElec&IC', tE, pcols, lchnk)
       call outfld ('TIon&IC' , tI, pcols, lchnk)
    endif

    qionelec(:ncol,:) = ptend_loc%s(:ncol,:)/cpairv(:ncol,:,lchnk)
    call outfld ('QIonElec' , qionelec, pcols, lchnk)

    call physics_ptend_sum(ptend_loc, ptend, ncol)

    call t_stopf ('ion_electron_temp_tend')

  end subroutine ion_electron_temp_tend

!===============================================================================

  subroutine update_istate(state, pbuf, istate, teTiBot)
  
    !---------------------------------------------------------------------------------------
    ! Time independent initialization for extended ionosphere simulation called in phys_init
    ! of physpkg module which is called in cam_comp module
    !---------------------------------------------------------------------------------------
    use mo_apex,          only : bnorth, beast, bdown             ! Magnetic field components
    use time_manager,     only : get_curr_calday                  ! Routine to get current calendar day
    use physconst,        only : rairv, mbarv, rearth             ! Constituent dependent rair and mbar
    use ref_pres,         only : press_lim_idx                    
    use orbit,            only : zenith
    
    use short_lived_species, only : slvd_index,slvd_pbf_ndx => pbf_idx ! Routines to access short lived species     

    type(physics_buffer_desc), pointer  :: pbuf(:)             ! physics buffer
    type(physics_state), intent(in),    target :: state        ! physics state structure
    type(ionos_state),   intent(inout), target :: istate       ! ionosphere state structure

    integer, intent(out) :: teTiBot  ! bottom of ionosphere calculations 

!---------------------------Local storage-------------------------------    
    integer,parameter :: nCnst = 9     ! Number of species needed from state%q or pbuf

    integer :: lchnk      ! Chunk number 
    integer :: ncol       ! Number of columns in current chunk 

    integer :: indxIR     ! pbuf index for ionization rates
    integer :: indxAIPRS  ! pbuf index for aurora ion production rate sum

    integer :: indxCnst   ! Constituent index used in cslculating densities

    integer :: indxSLvd   ! index of pbuf to access short lived species
    integer :: sIndx      ! index of adv_mass for any short lived species to access constituent mass

    integer :: iVer       ! Counter for vertical loops
    integer :: iCol       ! Counter for column loops
    integer :: iIonR      ! Counter for ionization rates loops
    integer :: iCnst      ! Counter for constituent loop
   
    integer :: indxSP     ! pbuf index for Pedersen Conductivity
    integer :: indxSH     ! pbuf index for Hall Conductivity

    real(r8), parameter :: teTiBotPres   = 50._r8 ! Pressure above which electron/ion temperature are calculated
                                                  ! in WACCM-X. (Pa)

    character(len = 3), dimension(nCnst) :: cCnst

    real(r8), dimension(:,:), pointer :: sigma_ped    ! Pointer to Pedersen Conductivity in pbuf (siemens/m) from module iondrag
    real(r8), dimension(:,:), pointer :: sigma_hall   ! Pointer to Hall Conductivity in pbuf (siemens/m)

    real(r8), dimension(:,:), pointer :: mmrP         ! Pointer to access short lived species in pbuf

    real(r8), dimension(:),pointer    :: geoLatR  ! Latitude (radians)  Make ncol because zenith aurora are ncol
    real(r8), dimension(:),pointer    :: geoLonR  ! Longitude (radians)

    real(r8), dimension(:,:),pointer  :: pMid     ! Midpoint pressure (Pa)
    real(r8), dimension(:,:),pointer  :: tN       ! Neutral temperature (K)

    real(r8), dimension(:,:),pointer  :: tNInt   ! Interface Temperture (K)

    real(r8), dimension(:),pointer    :: cosZenAngR            ! cosine of zenith angle (radians)
    real(r8), dimension(:),pointer    :: zenAngD               ! zenith angle (degrees)

    real(r8), dimension(:,:),pointer  :: bNorth3d  ! northward component of magnetic field units?
    real(r8), dimension(:,:),pointer  :: bEast3d   ! eastward component of magnetic field
    real(r8), dimension(:,:),pointer  :: bDown3d          ! downward component of magnetic field

    real(r8), dimension(pcols,pver)   :: sourceR       ! R term of source g4 calculation
    real(r8), dimension(pcols,pver)   :: sourceEff     ! Efficiency term of source g4 calculation
    
    real(r8), dimension(:,:),pointer  :: rairvi        ! Constituent dependent gas constant
 
    real(r8), dimension(:,:),pointer  :: dipMag  ! dip angle for each column (radians)

    real(r8), parameter :: rMassN2 = 28._r8       ! N2 molecular weight kg/kmol
    real(r8) :: rMass   ! Constituent molecular weight kg/kmol

    real(r8), dimension(:,:),pointer  :: mmrN2    ! N2 mass mixing ratio kg/kg
    real(r8), dimension(pcols,pver)   :: mmrO2    ! O2 mass mixing ratio kg/kg
    real(r8), dimension(pcols,pver)   :: mmrO1    ! O mass mixing ratio kg/kg

    real(r8), dimension(pcols,pver)   :: mmr      ! Constituent mass mixing ratio kg/kg

    real(r8), dimension(:,:),pointer  :: ndensN2  ! N2 number density (cm-3)
    real(r8), dimension(pcols,pver)   :: ndensO2  ! O2 number density (cm-3)
    real(r8), dimension(pcols,pver)   :: ndensO1  ! O number density (cm-3)
    real(r8), dimension(pcols,pver)   :: ndensN1  ! N number density  (cm-3)
    real(r8), dimension(pcols,pver)   :: ndensE   ! E electron number density (cm-3)

    real(r8), dimension(:,:)  ,pointer :: ndens    ! Constituent number density  (cm-3)

    real(r8), dimension(:,:,:),pointer :: ionRates     ! Pointer to ionization rates for O+,O2+,N+,N2+,NO+ in pbuf (s-1) 
                                                       !                               (from modules mo_jeuv and mo_jshort)

    real(r8), dimension(:,:,:),pointer :: ionPRates    ! ionization rates temporary array (s-1 cm-3)
    real(r8), dimension(:,:)  ,pointer :: sumIonPRates ! Sum of ionization rates for O+,O2+,N+,N2+,NO+ (s-1 cm-3)

    real(r8), dimension(:,:)  ,pointer :: aurIPRateSum ! Auroral ion production sum for O2+,O+,N2+
                                                       ! (s-1 cm-3 from module mo_aurora)

    real(r8), dimension(:,:)  ,pointer :: sourceg4     ! g4 source term for electron/ion temperature update

    real(r8)                          :: calDay                  ! current calendar day

   !real(r8), dimension(pcols,pver)   :: tempout         ! temporary scratch output array

    real(r8), dimension(pcols,pver)    :: zGeom        ! Geometric altitude (cm)
    real(r8), dimension(pcols,pver)    :: zThickness   ! Geometric altitude thickness (cm)
    real(r8), dimension(pcols)         :: eColDens     ! Electron column density (TECU = 1E16 m-2)

!--------------------------------------------------------------------------------

    sourceR           = 0._r8  
    sourceEff         = 0._r8

    mmrN2 => istate%n2_mmr
    mmrN2             = 0._r8
    mmrO2             = 0._r8
    mmrO1             = 0._r8
    mmr               = 0._r8

    ndensO2(:,:)      = 0._r8
    ndensO1(:,:)      = 0._r8
    ndensN1(:,:)      = 0._r8
    ndensE(:,:)       = 0._r8

    sourceR(:,:)      = 0._r8
    sourceEff(:,:)    = 0._r8
    
    !tempout(:,:)      = 0._r8
    
    !--------------------------------------------------------------------------------------
    !  Get lchnk from state 
    !--------------------------------------------------------------------------------------
    lchnk = state%lchnk
    ncol  = state%ncol

    !------------------------------------------------------------------------------------------------------
    !  Set the bottom of the ionosphere calculations at around 50 Pascals or 0.5 hectopascals(millibars).  
    !  teTiBotPres is in Pascals.
    !------------------------------------------------------------------------------------------------------
    teTiBot  = press_lim_idx(teTiBotPres, top=.false.)       

    !----------------------------------------------------------------
    !  Get latitude and longitude of each column in this chunk
    !----------------------------------------------------------------    
    geoLatR => state%lat(1:ncol)
    geoLonR => state%lon(1:ncol)

    !-------------------------------------------------------------------------------------------------------
    !  Need to get midpoint and interface pressure and neutral temperature from state structure (pcols,pver)
    !-------------------------------------------------------------------------------------------------------
    pMid => state%pmid(1:ncol,1:pver)
    tN   => state%t(1:ncol,1:pver)

    tNInt => istate%tNInt(1:ncol,1:pverp)
    cosZenAngR    => istate%cosZenAngR(1:ncol)
    zenAngD       => istate%zenAngD(1:ncol)

    bNorth3d => istate%bNorth3d(1:ncol,1:pver)
    bEast3d  => istate%bEast3d(1:ncol,1:pver)
    bDown3d  => istate%bDown3d(1:ncol,1:pver)

    rairvi => istate%rairvi(1:ncol,1:pverp)

    dipMag  => istate%dipMag(1:ncol,1:pver)

    ndensN2 => istate%ndensN2(1:ncol,1:pver)

    ionPRates    => istate%ionPRates(1:ncol,1:pver,1:nIonRates)
    sumIonPRates => istate%sumIonPRates(1:ncol,1:pver)

    sourceg4 => istate%sourceg4(1:ncol,1:pver)

    !-------------------------------------------------------------------------------------
    !  Calculate neutral temperature on interface levels.  tN vertical dimension is pver
    !-------------------------------------------------------------------------------------   
    do iVer = 2, pver
 
      do iCol = 1, ncol

        tNInt(iCol,iVer) = 0.5_r8 * tN(iCol,iVer) + 0.5_r8 * tN(iCol,iVer-1)

      enddo
    enddo

    do iCol = 1, ncol
        tNInt(iCol,1) = 1.5_r8 * tNInt(iCol,2) - 0.5_r8 * tNInt(iCol,3) 
    enddo
    do iCol = 1, ncol
        tNInt(iCol,pverp) = 1.5_r8 * tNInt(iCol,pver) - 0.5_r8 * tNInt(iCol,pver-1) 
    enddo

    !--------------------------------------------------------------
    !  Get zenith angle
    !-------------------------------------------------------------- 
    calDay = get_curr_calday()    
    call zenith(calDay,geoLatR(1:ncol),geoLonR(1:ncol),cosZenAngR(1:ncol),ncol)

    do iCol = 1, ncol

      zenAngD(iCol) = ACOS(cosZenAngR(iCol)) * rads2Degs
    
    enddo

    !---------------------------------------------------------------------------------------
    !  Expand magnetic field components in vertical to make 3D, pcols,pver,begchunk:endchunk
    !  These are used in calculation of magnetic dip angle and magnetic declination angle so 
    !  store in local ionosphere module structure.
    !---------------------------------------------------------------------------------------
    do iVer = 1, pver
    
      do iCol = 1, ncol

        bNorth3d(iCol,iVer) = bnorth(iCol,lchnk)
        bEast3d(iCol,iVer) = beast(iCol,lchnk)
        bDown3d(iCol,iVer) = bdown(iCol,lchnk)

      enddo
    
    enddo

    !------------------------------------------------------------------------
    !  Get constituent dependent gas constant and derive on interface levels
    !------------------------------------------------------------------------        
    do iVer = 2, pver
      do iCol = 1, ncol
        rairvi(iCol,iVer) = 0.5_r8 * rairv(iCol,iVer-1,lchnk) + 0.5_r8 * rairv(iCol,iVer,lchnk)
      enddo
    enddo

    do iCol = 1, ncol
       rairvi(iCol,1) = 1.5_r8 * rairvi(iCol,2) - 0.5_r8 * rairvi(iCol,3)
    enddo
    do iCol = 1, ncol
       rairvi(iCol,pverp) = 1.5_r8 * rairvi(iCol,pver) - 0.5_r8 * rairvi(iCol,pver-1)
    enddo

    !-------------------------------------------------------------------------------
    !  Need to get dip angle from magnetic field components
    !-------------------------------------------------------------------------------     
    do iVer = 1, pver
      do iCol = 1, ncol
        dipMag(iCol,iVer) = ATAN(bDown3d(iCol,iVer) / SQRT(bNorth3d(iCol,iVer)**2 + bEast3d(iCol,iVer)**2))
        if (dipMag(iCol,iVer) < 0.17_r8 .and. dipMag(iCol,iVer) > 0._r8 ) dipMag(iCol,iVer) = 0.17_r8
        if (dipMag(iCol,iVer) > -0.17_r8 .and. dipMag(iCol,iVer) < 0._r8 ) dipMag(iCol,iVer) = 0.17_r8
      enddo
    enddo

    !-------------------------------------------------------------------------------------------
    !  Set up constituents to be accessed here from pbuf or state%q. 
    !-------------------------------------------------------------------------------------------
    cCnst = (/'O  ','O2 ','NO ','H  ','N  ','e  ','Op ','O2p','NOp'/)

    do iCnst = 1, nCnst

      !--------------------------------------
      !  Assign density to istate array
      !-------------------------------------- 
      if (cCnst(iCnst) == 'O  ') ndens => istate%ndensO1(1:ncol,1:pver)
      if (cCnst(iCnst) == 'O2 ') ndens => istate%ndensO2(1:ncol,1:pver) 
      if (cCnst(iCnst) == 'NO ') ndens => istate%ndensNO(1:ncol,1:pver) 
      if (cCnst(iCnst) == 'N  ') ndens => istate%ndensN1(1:ncol,1:pver) 
      if (cCnst(iCnst) == 'e  ') ndens => istate%ndensE(1:ncol,1:pver)
      if (cCnst(iCnst) == 'Op ') ndens => istate%ndensOp(1:ncol,1:pver)
      if (cCnst(iCnst) == 'O2p') ndens => istate%ndensO2p(1:ncol,1:pver)
      if (cCnst(iCnst) == 'NOp') ndens => istate%ndensNOp(1:ncol,1:pver)

      !-------------------------------------------------------------------------------------------
      !  Set flag and get field mmr whether each constituent is short-lived(pbuf) or not(state%q). 
      !-------------------------------------------------------------------------------------------
      call cnst_get_ind( TRIM(cCnst(iCnst)), indxCnst, abort=.false. )
      if (indxCnst < 0) then
         indxSlvd = slvd_index( TRIM(cCnst(iCnst)) )
         if (indxSLvd > 0) then
            call pbuf_get_field(pbuf, slvd_pbf_ndx, mmrP, start=(/1,1,indxSLvd/), kount=(/pcols,pver,1/) )
            mmr(1:ncol,1:pver) = mmrP(1:ncol,1:pver)
            sIndx = get_spc_ndx( TRIM(cCnst(iCnst)) )
            rMass  = adv_mass(sIndx)
         endif
      else
         mmr(1:ncol,1:pver) = state%q(1:ncol,1:pver,indxCnst)
         rMass  = cnst_mw(indxCnst)
      endif

      !--------------------------------------------------------------------------------------------------------------
      !  Need to get number density (cgs units) from mass mixing ratio.  mbarv is kg/mole, same as rMass units
      !  kg/kg * (kg/mole)/(kg/mole) * (Pa or N/m*m)/((Joules/K or N*m/K) * (K)) = m-3 * 1E-06 = cm-3
      !--------------------------------------------------------------------------------------------------------------- 
      ndens(1:ncol,1:pver)  = mmr(1:ncol,1:pver) * mbarv(1:ncol,1:pver,lchnk) / rMass * &
                                   pMid(1:ncol,1:pver) / (kboltz * tN(1:ncol,1:pver)) * 1.E-06_r8
 
      if (cCnst(iCnst) == 'O  ') then
        mmrO1(1:ncol,1:pver) = mmr(1:ncol,1:pver)
        ndensO1(1:ncol,1:pver) = ndens(1:ncol,1:pver)
      endif
      if (cCnst(iCnst) == 'O2 ') then
        mmrO2(1:ncol,1:pver) = mmr(1:ncol,1:pver)
        ndensO2(1:ncol,1:pver) = ndens(1:ncol,1:pver)
      endif
      if (cCnst(iCnst) == 'N  ') ndensN1(1:ncol,1:pver) = ndens(1:ncol,1:pver)
      if (cCnst(iCnst) == 'e  ') ndensE(1:ncol,1:pver) = ndens(1:ncol,1:pver)

      !----------------------------------------------------------------------------
      !  Calculate N2 density from O2 and O and assign to istate array
      !----------------------------------------------------------------------------
      if (iCnst == nCnst) then

        mmrN2(1:ncol,1:pver) = 1._r8 - (mmrO2(1:ncol,1:pver) + mmrO1(1:ncol,1:pver)) 
        mmrN2(1:ncol,1:pver) = MAX(1.e-20_r8,mmrN2(1:ncol,1:pver))
        ndensN2(1:ncol,1:pver) = mmrN2(1:ncol,1:pver) * mbarv(1:ncol,1:pver,lchnk) / rMassN2 * &
                                            pMid(1:ncol,1:pver) / (kboltz * tN(1:ncol,1:pver)) * 1.E-06_r8              
         
      endif
 
    enddo ! nCnst

    if (hist_fld_active('ElecColDens')) then
       !---------------------------------------
       ! Calculate electron column density 
       !---------------------------------------
       !------------------------------------------------------------------------------
       ! Convert geopotential altitude in meters to geometric altitude in centimeters
       !------------------------------------------------------------------------------
       zGeom(1:ncol,1:pver) = state%zm(1:ncol,1:pver) * (1._r8 + state%zm(1:ncol,1:pver) / rearth) * 100._r8

       !------------------------------------------------------------
       ! Calculate vertical thickness at each level in centimeters
       !------------------------------------------------------------
       do iVer = 2, pver-1

          zThickness(1:ncol,iVer) = (zGeom(1:ncol,iVer-1) - zGeom(1:ncol,iVer+1)) / 2._r8

       enddo

       zThickness(1:ncol,1) = (1.5_r8 * zThickness(1:ncol,2)) - (0.5_r8 * zThickness(1:ncol,3))
       zThickness(1:ncol,pver) = (1.5_r8 * zThickness(1:ncol,pver-1)) - (0.5_r8 * zThickness(1:ncol,pver-2))       

       !----------------------------------------------------------------------------------
       ! Calculate electron column density converting from cm-2 to TEC units (1E16 m-2) 
       ! and make available for history output
       !----------------------------------------------------------------------------------
       eColDens(1:ncol) = sum(ndensE(1:ncol,:) * zThickness(1:ncol,:), dim=2) / 1.E12_r8

       call outfld('ElecColDens', eColDens, pcols, lchnk)
    endif

    !------------------------------------------------------------------------------------
    ! Get ionization rates from physics buffer which were calculated in mo_jeuv and 
    ! mo_jshort modules.  Rates array dimensions are pcols, pver, nIonRates.  Units s-1
    !------------------------------------------------------------------------------------
    indxIR = pbuf_get_index( 'IonRates' )
    call pbuf_get_field(pbuf, indxIR, ionRates)

    !----------------------------------------------------------------------------------------------
    !  Need to convert these ionization rates to ion production rates by multiplying number density 
    !  of neutral species appropriate from reactions in mo_jeuv(jeuv) and mo_jshort(jshort)(for NO)  
    !----------------------------------------------------------------------------------------------         
    do iVer = 1, pver
       do iCol = 1, ncol

          do iIonR = 1, nIonRates
             IF (iIonR <= 3) then
                ionPRates(iCol,iVer,iIonR) = ionRates(iCol,iVer,iIonR) * ndensO1(iCol,iVer)
             else IF (iIonR == 4) then
                ionPRates(iCol,iVer,iIonR) = ionRates(iCol,iVer,iIonR) * ndensN1(iCol,iVer)
             else IF ((iIonR == 5) .OR. (iIonR >= 7 .AND. iIonR <= 9)) then

                ionPRates(iCol,iVer,iIonR) = ionRates(iCol,iVer,iIonR) * ndensO2(iCol,iVer)
             else IF (iIonR == 6 .OR. iIonR == 10 .OR. iIonR == 11) then
                ionPRates(iCol,iVer,iIonR) = ionRates(iCol,iVer,iIonR) * ndensN2(iCol,iVer)
             endif
          enddo

        !----------------------------------------------
        !  Sum ion production rates all reactions
        !----------------------------------------------   
        sumIonPRates(iCol,iVer) = SUM(ionPRates(iCol,iVer,1:11))

       enddo
    enddo
    if (.not.steady_state_ion_elec_temp) then
       !-------------------------------------------------------------------------------------------
       ! Get aurora ion production rate sum from physics buffer which were calculated in mo_aurora 
       ! module.  Rate array dimensions are pcols, pver.  Units s-1 cm-3
       !-------------------------------------------------------------------------------------------
       indxAIPRS = pbuf_get_index( 'AurIPRateSum' )
       call pbuf_get_field(pbuf, indxAIPRS, aurIPRateSum)

       !-------------------------------------------------------------------------------------------------
       !  Calculate electron heating rate which is a source in electron/ion temperature derivation
       !-------------------------------------------------------------------------------------------------
       do iVer = 1, teTiBot
          do iCol = 1, ncol
             sourceR(iCol,iVer) = LOG( ndensE(iCol,iVer) / (ndensO2(iCol,iVer) + ndensN2(iCol,iVer) + &
                  0.1_r8 * ndensO1(iCol,iVer)) )
             sourceEff(iCol,iVer) = EXP( -(12.75_r8 + 6.941_r8 * sourceR(iCol,iVer) + 1.166_r8 * sourceR(iCol,iVer)**2 + &
                  0.08043_r8 * sourceR(iCol,iVer)**3 + 0.001996_r8 * sourceR(iCol,iVer)**4) )

             !-------------------------------------------------------------------------------
             !  Calculate g4 source term for electron temperature update
             !-------------------------------------------------------------------------------         
             sourceg4(iCol,iVer) = (sumIonPRates(iCol,iVer) + aurIPRateSum(iCol,iVer)) * sourceEff(iCol,iVer)

          enddo

       enddo

       call outfld ('SOURCER'     , sourceR     , pcols, lchnk)
       call outfld ('SOURCEEff'   , sourceEff   , pcols, lchnk)
       call outfld ('AURIPRATESUM', aurIPRateSum, pcols, lchnk)

       !----------------------------------------------------------------------------------------------
       ! Get Pedersen and Hall Conductivities from physics buffer which were calculated in iondrag 
       ! module.  Conductivity array dimensions are pcols, pver
       !-------------------------------------------------------------------------------
       indxSP = pbuf_get_index( 'PedConduct' )
       indxSH = pbuf_get_index( 'HallConduct' )
       call pbuf_get_field(pbuf, indxSP, sigma_ped)
       call pbuf_get_field(pbuf, indxSH, sigma_hall)

    endif

    return

  end subroutine update_istate
!
!===============================================================================

  subroutine update_teti(state, dSETendIn, dSETendOut, ztodt, istate, tE, tI, teTiBot)

  !-----------------------------------------------------------------------
  ! Routine to compute the electron and ion temperature 
  !-----------------------------------------------------------------------

    use physconst, only : gravit ! Gravity (m/s2)
    use physconst, only : rairv, mbarv  ! Constituent dependent rair and mbar
    use mo_apex, only: alatm
    
!------------------------------Arguments--------------------------------

    type(physics_state),   intent(in), target    :: state    ! physics state structure
    type(ionos_state),     intent(in), target    :: istate   ! ionosphere state structure

    real(r8), dimension(pcols,pver),   intent(in)      :: dSETendIn    ! dry static energy tendency
    real(r8), dimension(pcols,pver),   intent(out)     :: dSETendOut   ! dry static energy tendency

    real(r8), intent(in)                               :: ztodt     ! physics time step

    real(r8), dimension(:,:), pointer, intent(inout)   :: tE        ! Pointer to electron temperature in pbuf (K) 
    real(r8), dimension(:,:), pointer, intent(inout)   :: tI        ! Pointer to ion temperature in pbuf (K) 

    integer, intent(in) :: teTiBot  ! bottom of ionosphere calculations 

!---------------------------Local storage-------------------------------
    integer, parameter  :: maxIter  = 6                 ! maximum number of iterations to solve for electron/ion temperature
 
    integer :: lchnk                                    ! Chunk number 
    integer :: ncol                                     ! Number of atmospheric columns 
    integer :: teTiBotP                                 ! bottom of ionosphere calculations plus one more level 

    integer :: iVer                                     ! Counter for vertical loops
    integer :: iCol                                     ! Counter for column loops
    integer :: iter                                     ! Counter for iteration loop

    real(r8), parameter :: Kec1   = 7.5E5_r8            ! c1 constant for calculation of electron conductivity(Ke)
    real(r8), parameter :: Kec2   = 3.22E4_r8           ! c2 constant for calculation of electron conductivity(Ke)
    real(r8), parameter :: stepweight  = 1.0_r8         ! weight of previous and current times step for diagonals
    real(r8), parameter :: sToQConv  = 6.24E15_r8       ! Conversion from J/kg/s to ev/g/s

    real(r8), parameter :: lossc5  = 1.21E-4_r8         ! c5 constant needed for loss term g3 for electron temperature update
    real(r8), parameter :: lossc7  = 3.6E-2_r8          ! c7 constant needed for loss term g3 for electron temperature update
    real(r8), parameter :: lossc9  = 5.7E-4_r8          ! c9 constant needed for loss term g3 for electron temperature update
    real(r8), parameter :: lossc13 = 7.E-5_r8           ! c13 constant needed for loss term g3 for electron temperature update

    real(r8), parameter :: lossc4pCoef  = 1.77E-19_r8 
    real(r8), parameter :: lossc6pCoef  = 1.21E-18_r8 
    real(r8), parameter :: lossc8pCoef  = 7.9E-19_r8
    real(r8), parameter :: lossc10pCoef = 1.3E-4_r8 
    real(r8), parameter :: lossc11pCoef = 3.125E-21_r8
    real(r8), parameter :: lossc12pCoef = 3.4E-12_r8
    real(r8), parameter :: lossc14pCoef = 1.57E-12_r8 
    real(r8), parameter :: lossc15pCoef = 2.9E-14_r8
    real(r8), parameter :: lossc16pCoef = 6.9E-14_r8
    real(r8), parameter :: lossc3pC1 = 3.2E-8_r8
    real(r8), parameter :: lossc3pC2 = 15._r8
    real(r8), parameter :: lossc3pC3 = 0.53_r8

    real(r8), parameter :: losscinCoef1 = 6.6e-14_r8
    real(r8), parameter :: losscinCoef2 = 5.8e-14_r8
    real(r8), parameter :: losscinCoef3 = 0.21e-14_r8
    real(r8), parameter :: losscinCoef4 = 5.9e-14_r8
    real(r8), parameter :: losscinCoef5 = 5.45e-14_r8
    real(r8), parameter :: losscinCoef6 = 4.5e-14_r8
    real(r8), parameter :: losscinCoef7 = 5.8e-14_r8
    real(r8), parameter :: losscinCoef8 = 0.14e-14_r8
    real(r8), parameter :: losscinCoef9 = 4.4e-14_r8
    
    real(r8), parameter :: FeDCoef1 = -5.0E+7_r8
    real(r8), parameter :: FeDCoef2 = 4.0E+7_r8

    real(r8), parameter :: losscACoef1 = 5.71E-8_r8
    real(r8), parameter :: losscACoef2 = -3352.6_r8
    real(r8), parameter :: losscACoef3 = 2.0E-7_r8
    real(r8), parameter :: losscACoef4 = -4605.2_r8
    real(r8), parameter :: losscACoef5 = 2.53E-6_r8 
    real(r8), parameter :: losscACoef6 = -17620._r8

    real(r8), parameter :: loss10pCoef = 3200._r8
    real(r8), parameter :: lossc12pC1  = 0.4_r8
    real(r8), parameter :: lossc12pC2  = 150._r8

    real(r8), parameter :: losscf2dC1 = 2.4E+4_r8
    real(r8), parameter :: losscf2dC2 = 0.3_r8
    real(r8), parameter :: losscf2dC3 = 1500._r8
    real(r8), parameter :: losscf2dC4 = 1.947E-5_r8
    real(r8), parameter :: losscf2dC5 = 4000._r8

    real(r8), parameter :: losscf2C1 = 3000._r8

    real(r8), parameter :: losscf3c1 = -22713._r8

    real(r8), parameter :: f1Ted1C1 = 2.82E-17_r8
    real(r8), parameter :: f1Ted1C2 = 3.41E-21_r8

    real(r8), parameter :: f1Ted2C1 = 2.2E-16_r8
    real(r8), parameter :: f1Ted2C2 = 7.92E-18_r8

    real(r8), parameter :: f1Ted3C1 = 1.1E-16_r8
    real(r8), parameter :: f1Ted3C2 = 5.7E-4_r8

    real(r8) :: wrk1                                    ! 2/3/kboltz_ev
    real(r8) :: FeDB                                    ! B term of electron heat flux of UB
    real(r8) :: FeD                                     ! Day time flux
    real(r8) :: FeN                                     ! Night time flux
    real(r8) :: f1Ted1                                  ! d1 of f1(Te) calculation used to get electron conductivity
    real(r8) :: f1Ted2                                  ! d2 of f1(Te) calculation used to get electron conductivity
    real(r8) :: f1Ted3                                  ! d3 of f1(Te) calculation used to get electron conductivity
    real(r8) :: f1Te

    real(r8), dimension(:,:), pointer   :: pMid         ! Midpoint pressure (Pa)
    real(r8), dimension(:,:), pointer   :: tN           ! Neutral temperature (K)
    real(r8), dimension(pcols,pver)     :: tEPrevI      ! Electron temperature from previous iteration (K)

    real(r8), dimension(:,:), pointer   :: pInt         ! Interface pressure (Pa)
    real(r8), dimension(:,:), pointer   :: tNInt        ! Interface Temperture (K)
    real(r8), dimension(:,:), pointer   :: rairvi       ! Constituent dependent gas constant on interface levels

    real(r8), dimension(:,:), pointer    :: ndensN2     ! N2 number density (cm-3)
    real(r8), dimension(:,:), pointer    :: ndensO2     ! O2 number density (cm-3)
    real(r8), dimension(:,:), pointer    :: ndensO1     ! O number density (cm-3)
    real(r8), dimension(:,:), pointer    :: ndensE      ! E electron number density (cm-3)
    real(r8), dimension(:,:), pointer    :: ndensOp     ! O plus number density (cm-3)
    real(r8), dimension(:,:), pointer    :: ndensO2p    ! O2 plus ion number density (cm-3)
    real(r8), dimension(:,:), pointer    :: ndensNOp    ! NO plus ion number density  (cm-3)

    real(r8), dimension(:,:), pointer   :: sourceg4     ! g4 source term for electron/ion temperature update
 
    real(r8), dimension(:,:), pointer   :: dipMag       ! dip angle for each column (radians)
    real(r8), dimension(pcols)          :: dlatm        ! magnetic latitude of each phys column (degrees)

    real(r8), dimension(:),   pointer   :: zenAngD      ! zenith angle (degrees)

    real(r8), dimension(pcols)          :: FeUB         ! electron heat flux at upper boundary
 
    real(r8), dimension(pver)           :: sqrtTE       ! Square root of electron temperature
 
    real(r8), dimension(pver)           :: Ke           ! electron conductivity

    real(r8), dimension(pverp)          :: Kei          ! electron conductivity interface levels

    real(r8), dimension(pcols,pver)     :: lossc4p      ! c4 prime of Lc(eN2) component of loss term
    real(r8), dimension(pcols,pver)     :: lossceN2     ! Lc(eN2) component of loss term equation
 
    real(r8), dimension(pcols,pver)     :: lossc6p      ! c6 prime of Lc(eO2) component of loss term equation
    real(r8), dimension(pcols,pver)     :: lossceO2     ! Lc(eO2) component of loss term equation
 
    real(r8), dimension(pcols,pver)     :: lossc8p      ! c8 prime of Lc(eO) component of loss term equation
    real(r8), dimension(pcols,pver)     :: lossceO1     ! Lc(eO) component of loss term equation

    real(r8), dimension(pcols,pver)     :: lossc10p     ! c10 prime of Lc(eN2) component of loss term equation
    real(r8), dimension(pcols,pver)     :: losscA       ! A of Lc(eN2)v component of loss term equation
    real(r8), dimension(pcols,pver)     :: tENDiff      ! Difference between electron and neutral temperatures
    real(r8), dimension(pcols,pver)     :: lossceN2v    ! Lc(eN2)v component of loss term equation

    real(r8), dimension(pcols,pver)     :: lossc11p     ! c11 prime of Lc(eO2)v component of loss term equation
    real(r8), dimension(pcols,pver)     :: lossceO2v    ! Lc(eO2)v component of loss term equation
 
    real(r8), dimension(pcols,pver)     :: lossc12p     ! c12 prime of Lc(eO)f component of loss term equation
    real(r8), dimension(pcols,pver)     :: lossceOf     ! Lc(eO)f component of loss term equation
 
    real(r8), dimension(pcols,pver)     :: lossc14p     ! c14 prime of Lc(eO)1D component of loss term equation
    real(r8), dimension(pcols,pver)     :: losscf2d     ! d of f2 of Lc(eO)1D component of loss term equation
    real(r8), dimension(pcols,pver)     :: losscf2      ! f2 of Lc(eO)1D component of loss term equation
    real(r8), dimension(pcols,pver)     :: losscf3      ! f3 of Lc(eO)1D component of loss term equation
    real(r8), dimension(pcols,pver)     :: lossceO1D    ! Lc(eO)1D component of loss term equation
 
    real(r8), dimension(pcols,pver)     :: lossc15p     ! c15 prime of Lc(eN2)Rot component of loss term equation
    real(r8), dimension(pcols,pver)     :: lossceN2Rot  ! Lc(eN2)Rot component of loss term equation
 
    real(r8), dimension(pcols,pver)     :: lossc16p     ! c16 prime of Lc(eO2)Rot component of loss term equation
    real(r8), dimension(pcols,pver)     :: lossceO2Rot  ! Lc(eO2)Rot component of loss term equation
 
    real(r8), dimension(pcols,pver)     :: lossc3p      ! c3 prime of Lc(ei) component of loss term equation
    real(r8), dimension(pcols,pver)     :: losscei      ! Lc(ei) component of loss term equation
    real(r8), dimension(pcols,pver)     :: losscin      ! ion-neutral heating coeff.

    real(r8), dimension(pcols,pver)     :: lossg3       ! g3 loss term for Te tendency

    real(r8), dimension(pcols,pverp)    :: delZi        ! Delta z: interfaces
    real(r8), dimension(pcols,pver)     :: delZ         ! Delta z: midpoints
 
    real(r8), dimension(pcols,pver)     :: qjoule       ! joule heating
    real(r8), dimension(pcols,pver)     :: qen          ! electron-neutral heating
    real(r8), dimension(pcols,pver)     :: qei          ! electron-ion Coulomb heating
    real(r8), dimension(pcols,pver)     :: qin          ! ion-neutral heating
    real(r8), dimension(pcols,pver)     :: rho          ! mass density

    real(r8), dimension(pcols,pver)     :: wrk2

    real(r8), dimension(teTiBot)        :: subdiag      ! subdiagonal values for Te tendency solving
    real(r8), dimension(teTiBot)        :: superdiag    ! superdiagonal values for Te tendency solving
    real(r8), dimension(teTiBot)        :: diag         ! diagonal values for Te tendency solving
    real(r8), dimension(teTiBot)        :: rHS          ! RHS of electron temperature update
    real(r8), dimension(teTiBot)        :: tETemp       ! temporary electron temperature array for input to tridag

    logical, dimension(pcols)           :: colConv      ! flag for column converging
    logical                             :: converged    ! Flag for convergence in electron temperature
                                                        ! calculation iteration loop

    !--------------------------------------------------------------------------------------------------------- 
    !  Initialize arrays to zero and column convergence logical to .false.
    !---------------------------------------------------------------------------------------------------------

    sqrtTE(:)           = 0._r8
    Ke(:)               = 0._r8
    Kei(:)              = 0._r8
    lossc4p(:,:)        = 0._r8
    lossceN2(:,:)       = 0._r8
    lossc6p(:,:)        = 0._r8
    lossceO2(:,:)       = 0._r8
    lossc8p(:,:)        = 0._r8
    lossceO1(:,:)       = 0._r8
    lossc10p(:,:)       = 0._r8
    losscA(:,:)         = 0._r8
    tENDiff(:,:)        = 0._r8
    lossceN2v(:,:)      = 0._r8
    lossc11p(:,:)       = 0._r8
    lossceO2v(:,:)      = 0._r8
    lossc12p(:,:)       = 0._r8
    lossceOf(:,:)       = 0._r8
    lossc14p(:,:)       = 0._r8
    losscf2d(:,:)       = 0._r8
    losscf2(:,:)        = 0._r8
    losscf3(:,:)        = 0._r8
    lossceO1D(:,:)      = 0._r8
    lossc15p(:,:)       = 0._r8
    lossceN2Rot(:,:)    = 0._r8
    lossc16p(:,:)       = 0._r8
    lossceO2Rot(:,:)    = 0._r8
    lossc3p(:,:)        = 0._r8
    losscei(:,:)        = 0._r8
    losscin(:,:)        = 0._r8
    lossg3(:,:)         = 0._r8
    delZi(:,:)          = 0._r8
    delZ(:,:)           = 0._r8 
    subDiag(:)          = 0._r8         
    superDiag(:)        = 0._r8        
    diag(:)             = 0._r8 
    rHS(:)              = 0._r8 
    teTemp(:)           = 0._r8 
    qjoule(:,:)         = 0._r8
    qei(:,:)            = 0._r8
    qen(:,:)            = 0._r8
    qin(:,:)            = 0._r8
    rho(:,:)            = 0._r8
    dSETendOut          = 0._r8
    colConv(:)          = .false.

    !--------------------------------------------------------------------------------------
    !  Get lchnk and ncol from state
    !--------------------------------------------------------------------------------------
    lchnk = state%lchnk
    ncol  = state%ncol
    
    !-------------------------------------------
    !  Calculate some commonly used variables
    !-------------------------------------------
    wrk1 = 2._r8 / 3._r8/ kboltz_ev
    teTiBotP = teTiBot + 1

    !-------------------------------------------------------------------------------------------------------
    !  Need to get midpoint and interface pressure and neutral temperature from state structure (ncol,teTiBot)
    !-------------------------------------------------------------------------------------------------------
    pMid  => state%pmid(1:ncol,1:pver)
    tN    => state%t(1:ncol,1:pver)
    rho(1:ncol,1:pver) = pMid(1:ncol,1:pver)/rairv(1:ncol,1:pver,lchnk)/tN(1:ncol,1:pver) * 1.E-3_r8     ! convert to g/cm3

    qjoule(1:ncol,1:teTiBot) = dSETendIn(1:ncol,1:teTiBot) * sToQConv     ! convert from J/kg/s to ev/g/s

    pInt    => state%pint(1:ncol,1:pverp)
    tNInt   => istate%tNInt(1:ncol,1:pverp)        
    rairvi  => istate%rairvi(1:ncol,1:pverp)

    !----------------------------------------------------------------
    !  Get variables needed from the ionosphere state structure
    !----------------------------------------------------------------
    ndensO2  => istate%ndensO2(1:ncol,1:pver) 
    ndensO1  => istate%ndensO1(1:ncol,1:pver)
    ndensE   => istate%ndensE(1:ncol,1:pver)  
    ndensOp  => istate%ndensOp(1:ncol,1:pver) 
    ndensO2p => istate%ndensO2p(1:ncol,1:pver)
    ndensNOp => istate%ndensNOp(1:ncol,1:pver)
    ndensN2  => istate%ndensN2(1:ncol,1:pver) 

    sourceg4 => istate%sourceg4(1:ncol,1:pver)

    dipMag   => istate%dipMag(1:ncol,1:pver)

    zenAngD  => istate%zenAngD(1:ncol) 
      
    !-------------------------------------------------------------------------------------------------------------------
    !  Set electron temperature limits 
    !-------------------------------------------------------------------------------------------------------------------
    tE(1:ncol,1:pver) = MAX(tN(1:ncol,1:pver),tE(1:ncol,1:pver))
    tE(1:ncol,1:pver) = MIN(temax,tE(1:ncol,1:pver))

    tI(1:ncol,1:pver) = MAX(tN(1:ncol,1:pver),ti(1:ncol,1:pver))
    tI(1:ncol,1:pver) = MIN(ti(1:ncol,1:pver),tE(1:ncol,1:pver))

    ! set Te and Ti to Tn below the levels where this module applies
    tE(1:ncol,teTiBotP:pver) = tN(1:ncol,teTiBotP:pver)
    tI(1:ncol,teTiBotP:pver) = tN(1:ncol,teTiBotP:pver)

    wrk2(1:ncol,1:teTiBot) =  ndensE(1:ncol,1:teTiBot)/wrk1/(SIN(dipMag(1:ncol,1:teTiBot)))**2._r8
    
    dlatm(:ncol) = rads2Degs* alatm(:ncol,lchnk)
    
    !-----------------------------------------------------------------------------
    !  Get terms needed for loss term g3 for electron temperature update which do 
    !  not need to be updated in iteration loop.  
    !-----------------------------------------------------------------------------
    do iCol = 1, ncol

      if (.not. colConv(iCol)) then
        do iVer = 1, teTiBot

          lossc4p(iCol,iVer)  = lossc4pCoef * ndensN2(iCol,iVer) * ndensE(iCol,iVer)     ! e-N2 elastic collision
          lossc6p(iCol,iVer)  = lossc6pCoef * ndensO2(iCol,iVer) * ndensE(iCol,iVer)     ! e-O2 elastic collision
          lossc8p(iCol,iVer)  = lossc8pCoef * ndensO1(iCol,iVer) * ndensE(iCol,iVer)     ! e-O elastic collision
          lossc10p(iCol,iVer) = lossc10pCoef * ndensN2(iCol,iVer) * ndensE(iCol,iVer)    ! e-N2(vib)
          lossc11p(iCol,iVer) = lossc11pCoef * ndensO2(iCol,iVer) * ndensE(iCol,iVer)    ! e-O2(vib)
          lossc12p(iCol,iVer) = lossc12pCoef * ndensO1(iCol,iVer) * ndensE(iCol,iVer)    ! e-O (fine)
          lossc14p(iCol,iVer) = lossc14pCoef * ndensO1(iCol,iVer) * ndensE(iCol,iVer)    ! e-O(1D)
          lossc15p(iCol,iVer) = lossc15pCoef * ndensN2(iCol,iVer) * ndensE(iCol,iVer)    ! e-N2(rot)
          lossc16p(iCol,iVer) = lossc16pCoef * ndensO2(iCol,iVer) * ndensE(iCol,iVer)    ! e-O2(rot)
          lossc3p(iCol,iVer)  = lossc3pC1 * lossc3pC2 * (ndensOP(iCol,iVer) + &
                                0.5_r8 * ndensO2P(iCol,iVer) + lossc3pC3 * ndensNOP(iCol,iVer)) * ndensE(iCol,iVer)  ! e-i

          losscin(iCol,iVer) = (losscinCoef1*ndensN2(iCol,iVer) + losscinCoef2*ndensO2(iCol,iVer)                    &
                             + losscinCoef3*ndensO1(iCol,iVer)*SQRT(2._r8*tN(iCol,iVer)))*ndensOP(iCol,iVer)    &
                             +(losscinCoef4*ndensN2(iCol,iVer) + losscinCoef5*ndensO2(iCol,iVer)                   &
                             + losscinCoef6*ndensO1(iCol,iVer))*ndensNOP(iCol,iVer)                            &
                             +(losscinCoef7*ndensN2(iCol,iVer) + losscinCoef8*ndensO2(iCol,iVer)*SQRT(tN(iCol,iVer)) &
                             + losscinCoef9*ndensO1(iCol,iVer)) * ndensO2P(iCol,iVer)


        enddo !iVer loop

        !----------------------------------------------------------------------------------
        !  Calculate upper boundary heat flux 
        !----------------------------------------------------------------------------------
        if (ABS(dlatm(iCol)) < 40.0_r8) FeDB = 0.5_r8 * &
                                        (1._r8 + SIN(pi * (ABS(dlatm(iCol)) - 20.0_r8) /40.0_r8))

        if (ABS(dlatm(iCol)) >= 40.0_r8) FeDB = 1._r8


        FeD = FeDCoef1 * f107 * FeDB - FeDCoef2 * f107
        FeN = .5_r8 * FeD
        !---------------------------------------------------
        !  Set upper boundary condition for right hand side
        !---------------------------------------------------
        if (zenAngD(iCol) <= 80.0_r8) FeUB(iCol) = FeD
        if (zenAngD(iCol) > 80.0_r8 .AND. zenAngD(iCol) < 100.0_r8) FeUB(iCol) = 0.5_r8 * (FeD + FeN) &
                                                                         + 0.5_r8 * (FeD - FeN) * &
                                                                 COS(pi * ((zenAngD(iCol) - 80.0_r8) / 20.0_r8))
        if (zenAngD(iCol) >= 100.0_r8) FeUB(iCol) = FeN

        !------------------------------------------------------------------------------------------
        !  Calculate thickness terms for vertical derivative
        !------------------------------------------------------------------------------------------
        do iVer = 1, teTiBot

          delZ(iCol,iVer) = (pInt(iCol,iVer+1) - pInt(iCol,iVer)) * rairv(iCol,iVer,lchnk) * &
                                                                       tN(iCol,iVer) / pMid(iCol,iVer) / gravit

        enddo

        do iVer = 2, teTiBotP           ! Assuming teTiBotP < pverp
          delZi(iCol,iVer) = (pMid(iCol,iVer) - pMid(iCol,iVer-1)) * rairvi(iCol,iVer) * &
                                                                  tNInt(iCol,iVer) / pInt(iCol,iVer) / gravit
        enddo
        delZi(iCol,1) = 1.5_r8*delZi(iCol,2) - .5_r8*delZi(iCol,3)

        !----------------------------------------------------------
        !  Convert delZ variables from meters to centimeters
        !----------------------------------------------------------
        delZi(iCol,1:teTiBotP) = delZi(iCol,1:teTiBotP)*100._r8
        delZ(iCol,1:teTiBot) = delZ(iCol,1:teTiBot)*100._r8
  
      endif ! Column not converged

    enddo !iCol loop

    !-------------------------------------------------------------------------------------------------------
    !  Iterate to calculate new electron temperature. 
    !  Time splitting is used: first solve the heating/cooling equation, then solve the diffusion equations.
    !  Also, set convergence flag to false and iterate until true or 6 iterations, whichever comes first 
    !-------------------------------------------------------------------------------------------------------
    converged = .false. 
    iter = 0
    do while (.not. converged .and. iter < maxIter)
    
      !--------------------------------------------------------------------------------------------------------
      !  Increment iteration loop counter and save electron temperature from previous iteration for convergence 
      !  test at end of interation loop.  Also, take square root of electron temperature to be used later
      !--------------------------------------------------------------------------------------------------------        
      iter = iter + 1
      
      tEPrevI(1:ncol,1:teTiBot) = tE(1:ncol,1:teTiBot)

      !--------------------------------------------------------------------------------------------------------
      !  Loop over columns then vertical levels and call tridiagonal solver for each column to get electron 
      !  temperature
      !--------------------------------------------------------------------------------------------------------
      do iCol = 1, ncol

         if (.not. colConv(iCol)) then

            sqrtTE(1:teTiBot) = SQRT(tE(iCol,1:teTiBot))

            do iVer = 1, teTiBot

               !-----------------------------------------------------------------------------
               !  Get loss term g3 for electron temperature update.  Need to calculate 
               !  constituent dependent loss terms which make up g3
               !-----------------------------------------------------------------------------
               lossceN2(iCol,iVer) = lossc4p(iCol,iVer) * (1._r8 - lossc5 * tE(iCol,iVer)) * tE(iCol,iVer)
               lossceO2(iCol,iVer) = lossc6p(iCol,iVer) * (1._r8 + lossc7 * sqrtTE(iVer)) * sqrtTE(iVer)
               lossceO1(iCol,iVer) = lossc8p(iCol,iVer) * (1._r8 + lossc9 * tE(iCol,iVer)) * sqrtTE(iVer)

               if (tE(iCol,iVer) < 1000.0_r8) then
                  losscA(iCol,iVer) = losscACoef1 * EXP(losscACoef2 / tE(iCol,iVer)) 
               endif
               if (tE(iCol,iVer) >= 1000.0_r8 .AND. tE(iCol,iVer) <= 2000.0_r8) then
                  losscA(iCol,iVer) = losscACoef3 * EXP(losscACoef4 / tE(iCol,iVer))
               endif
               if (tE(iCol,iVer) > 2000.0_r8) then
                  losscA(iCol,iVer) = losscACoef5 * sqrtTE(iVer) * EXP(losscACoef6 / tE(iCol,iVer))
               endif

               tENDiff(iCol,iVer) = tE(iCol,iVer) -  tN(iCol,iVer)
               if (ABS(tENDiff(iCol,iVer)) < 0.1_r8) tENDiff(iCol,iVer) = 0.1_r8

               lossceN2v(iCol,iVer) = lossc10p(iCol,iVer) * losscA(iCol,iVer) * &
                    (1._r8 - EXP(loss10pCoef * (1._r8 / tE(iCol,iVer) - 1._r8 / tN(iCol,iVer)))) / tENDiff(iCol,iVer)
               lossceO2v(iCol,iVer) = lossc11p(iCol,iVer) * tE(iCol,iVer) * tE(iCol,iVer)
               lossceOf(iCol,iVer) = lossc12p(iCol,iVer) * (1._r8 - lossc13 * tE(iCol,iVer)) * &
                    (lossc12pC1 + lossc12pC2 / tE(iCol,iVer)) / tN(iCol,iVer)
               losscf2d(iCol,iVer) = losscf2dC1 + losscf2dC2 * (tE(iCol,iVer) - losscf2dC3) - &
                    losscf2dC4 * (tE(iCol,iVer) - losscf2dC3) * (tE(iCol,iVer) - losscf2dC5)
               losscf2(iCol,iVer) = losscf2d(iCol,iVer) * (1._r8 / losscf2C1 - 1._r8 / tE(iCol,iVer))
               losscf3(iCol,iVer) = losscf3c1 * (1._r8 / tN(iCol,iVer) - 1._r8 / tE(iCol,iVer))
               lossceO1D(iCol,iVer) = lossc14p(iCol,iVer) * EXP(losscf2(iCol,iVer)) * & 
                    (1._r8 - EXP(losscf3(iCol,iVer))) / tENDiff(iCol,iVer)
               lossceN2Rot(iCol,iVer) = lossc15p(iCol,iVer) / sqrtTE(iVer)
               lossceO2Rot(iCol,iVer) = lossc16p(iCol,iVer) / sqrtTE(iVer)

               losscei(iCol,iVer) = lossc3p(iCol,iVer) / tE(iCol,iVer)**1.5_r8

               !------------------------------------------------
               ! Loss term: lossg3*tE/sin(I)^2
               !------------------------------------------------
               lossg3(iCol,iVer) =  lossceN2(iCol,iVer) + lossceO2(iCol,iVer) + lossceO1(iCol,iVer) + lossceN2v(iCol,iVer)   &
                    + lossceO2v(iCol,iVer) + lossceOf(iCol,iVer) + lossceO1D(iCol,iVer)                        &
                    + lossceN2Rot(iCol,iVer) + lossceO2Rot(iCol,iVer)

            enddo !iVer loop

         endif ! Column not converged

      enddo ! End of column loop

      !-----------------------------------------------------
      !  Calculate thermal conductivity of electron gas   
      !-----------------------------------------------------
      do iCol = 1, ncol

         if (.not. colConv(iCol)) then

            sqrtTE(1:teTiBot) = SQRT(tE(iCol,1:teTiBot))

            do iVer = 1, teTiBot

               f1Ted1 = f1Ted1C1 * sqrtTE(iVer) - f1Ted1C2 * tE(iCol,iVer)**1.5_r8
               f1Ted2 = f1Ted2C1 + f1Ted2C2 * sqrtTE(iVer)
               f1Ted3 = f1Ted3C1 * (1._r8 + f1Ted3C2 * tE(iCol,iVer))

               f1Te = ndensN2(iCol,iVer) / ndensE(iCol,iVer) * f1Ted1 + ndensO2(iCol,iVer) / &
                    ndensE(iCol,iVer) * f1Ted2 + ndensO1(iCol,iVer) / ndensE(iCol,iVer) * f1Ted3

               !-----------------------------------------------------------------------------
               !  Calculate electron conductivity using parameters set in module and f1(Te)
               !-----------------------------------------------------------------------------
               Ke(iVer) = Kec1 * tE(iCol,iVer)**2.5_r8 / (1._r8 + Kec2 * tE(iCol,iVer)**2._r8 * f1Te)

            enddo !iVer loop

            !----------------------------------------------------------------------
            !  Get electron conductivity at interface levels to be used later
            !----------------------------------------------------------------------
            do iVer = 2,teTiBot
               Kei(iVer) = SQRT(Ke(iVer-1)*Ke(iVer))
            enddo
            Kei(1) = 1.5_r8*Kei(2)-.5_r8*Kei(3)
            Kei(teTiBotP) = 1.5_r8*Kei(teTiBot)-.5_r8*Kei(teTiBot-1)

            !------------------------------------------------------------------------------------------------------
            !  Derive subdiagonal, superdiagonal, and diagonal as input to solver for electron temperature tendency
            !------------------------------------------------------------------------------------------------------
            do iVer = 2, teTiBot-1
               subDiag(iVer) = -Kei(iVer) / delZi(iCol,iVer) / delZ(iCol,iVer)
               superDiag(iVer) = -Kei(iVer+1) / delZi(iCol,iVer+1) / delZ(iCol,iVer)
               diag(iVer) = wrk2(iCol,iVer)/ztodt + (lossg3(iCol,iVer)+losscei(iCol,iVer))/SIN(dipMag(iCol,iVer))**2._r8  &
                    -subDiag(iVer)-superDiag(iVer)
               rHS(iVer) = tE(iCol,iVer) * wrk2(iCol,iVer)/ztodt + sourceg4(iCol,iVer)/SIN(dipMag(iCol,iVer))**2._r8 &
                    +(lossg3(iCol,iVer)*tN(iCol,iVer)+losscei(iCol,iVer)*ti(iCol,iVer))/SIN(dipMag(iCol,iVer))**2._r8
            enddo !iVer loop

            !-------------------------------------------------------------------------------------
            !  Calculate diagonal, superdiagonal, and right hand side upper boundary values
            !-------------------------------------------------------------------------------------
            superDiag(1)  = -Kei(2) / delZi(iCol,2) / delZ(iCol,1)
            diag(1) = wrk2(iCol,1)/ztodt - superDiag(1)
            rHS(1) = tE(iCol,1) * wrk2(iCol,1) / ztodt - FeUB(iCol) / delZ(iCol,1)

            !---------------------------------------------------------------------------------------------
            !  Calculate subdiagonal, diagonal, superdiagonal, and right hand side lower boundary values
            !---------------------------------------------------------------------------------------------
            subDiag(teTiBot) = -Kei(teTiBot) / delZi(iCol,teTiBot) / delZ(iCol,teTiBot)
            superDiag(teTiBot) = -Kei(teTiBotP) / delZi(iCol,teTiBotP) / delZ(iCol,teTiBot)
            diag(teTiBot) = wrk2(iCol,teTiBot)/ztodt &
                 + (lossg3(iCol,teTiBot)+losscei(iCol,teTiBot))/SIN(dipMag(iCol,teTiBot))**2._r8 &
                 - subDiag(teTiBot)-superDiag(teTiBot)
            rHS(teTiBot) = tE(iCol,teTiBot) * wrk2(iCol,teTiBot)/ztodt &
                 + sourceg4(iCol,teTiBot)/SIN(dipMag(iCol,teTiBot))**2._r8  &
                 + (lossg3(iCol,teTiBot)*tN(iCol,teTiBot)+losscei(iCol,teTiBot)*ti(iCol,teTiBot)) &
                 / SIN(dipMag(iCol,teTiBot))**2._r8 - superDiag(teTiBot) * tN(iCol,teTiBotP)

            !-------------------------------------------------
            ! Call solver to get electron temperature update
            !-------------------------------------------------
            call tridag(subDiag,diag,superDiag,rHS,tETemp,teTiBot)

            tE(iCol,1:teTiBot) = tETemp(1:teTiBot)
            do iVer = 1,teTiBot
               tE(iCol,iVer) = min(temax,tE(iCol,iVer))
               tE(iCol,iVer) = max(tN(iCol,iVer),tE(iCol,iVer))
            enddo

            !---------------------------------------------------------------------------------------------------------
            !  Calculate ion temperature from electron temperature, ion-neutral and electron-ion loss terms, neutral 
            !  temperature, mass density and joule heating.  Set minimum value to neutral temperature and maximum 
            !  value to electron temperature for each column and vertical level
            !---------------------------------------------------------------------------------------------------------
            do iVer = 1,teTiBot
               ti(iCol,iVer) = (losscei(iCol,iVer) * tE(iCol,iVer) + losscin(iCol,iVer) * tN(iCol,iVer) +   &
                    rho(iCol,iVer) * qjoule(iCol,iVer) * mbarv(iCol,iVer,lchnk) / (mbarv(iCol,iVer,lchnk)+rMassOp)) / &
                    (losscei(iCol,iVer) + losscin(iCol,iVer))
               ti(iCol,iVer) = max(tN(iCol,iVer),ti(iCol,iVer))
               ti(iCol,iVer) = min(tE(iCol,iVer),ti(iCol,iVer))
            enddo

            !--------------------------------------------------------------------------------------------------------
            ! Check for convergence which is a change of electron temperature ratio to previous loop for all levels
            ! and columns of less than 0.05K.  Had to modify this to do convergence check on each column since 
            ! checking all columns in a chunk gives different answers depending on number of tasks and tasks per node.
            !--------------------------------------------------------------------------------------------------------
            if (ALL(ABS(tE(iCol,1:teTiBot) / tEPrevI(iCol,1:teTiBot) - 1._r8) < 0.05_r8)) then 

               colConv(iCol) = .true.

            endif

         endif ! Column not converged

      enddo ! iCol loop
      !--------------------------------------------------------------
      !  Check to see if all columns have converged and set flag
      !--------------------------------------------------------------
      if (ALL(colConv(1:ncol))) converged = .true.

   enddo ! End of iteration loop

   !--------------------------------------------------------------------------------------------------------
   ! Calculate electron-neutral heating and electron-ion Coulomb heating.  Then update dry static energy.
   !--------------------------------------------------------------------------------------------------------
   do iVer = 1, teTiBot
      do iCol = 1, ncol
         sqrtTE(iVer) = SQRT(tE(iCol,iVer))
         lossceN2(iCol,iVer) = lossc4p(iCol,iVer) * (1._r8 - lossc5 * tE(iCol,iVer)) * tE(iCol,iVer)
         lossceO2(iCol,iVer) = lossc6p(iCol,iVer) * (1._r8 + lossc7 * sqrtTE(iVer)) * sqrtTE(iVer)
         lossceO1(iCol,iVer) = lossc8p(iCol,iVer) * (1._r8 + lossc9 * tE(iCol,iVer)) * sqrtTE(iVer)
      enddo
   enddo

   qen(1:ncol,1:teTiBot) = (lossceN2(1:ncol,1:teTiBot)+lossceO2(1:ncol,1:teTiBot)+lossceO1(1:ncol,1:teTiBot)+ &
        lossceN2v(1:ncol,1:teTiBot)+lossceO2v(1:ncol,1:teTiBot)+lossceOf(1:ncol,1:teTiBot)+ &
        lossceO1D(1:ncol,1:teTiBot)+lossceN2Rot(1:ncol,1:teTiBot)+lossceO2Rot(1:ncol,1:teTiBot))*  &
        (tE(1:ncol,1:teTiBot)-tN(1:ncol,1:teTiBot)) / rho(1:ncol,1:teTiBot)
   qei(1:ncol,1:teTiBot) = losscei(1:ncol,1:teTiBot) * (tE(1:ncol,1:teTiBot)-ti(1:ncol,1:teTiBot)) / rho(1:ncol,1:teTiBot)
   qin(1:ncol,1:teTiBot) = losscin(1:ncol,1:teTiBot) * (tI(1:ncol,1:teTiBot)-tN(1:ncol,1:teTiBot)) / rho(1:ncol,1:teTiBot)

   dSETendOut(1:ncol,1:teTiBot) = (qei(1:ncol,1:teTiBot)+qen(1:ncol,1:teTiBot)) / sToQConv    ! J/kg/s 

   call outfld ('QEN', qen, pcols, lchnk)
   call outfld ('QEI', qei, pcols, lchnk)
   call outfld ('QIN', qin, pcols, lchnk)

   return

 end subroutine update_teti

!-----------------------------------------------------------------------
! Simple tridiagonal solver routine
!-----------------------------------------------------------------------

 SUBROUTINE tridag(a,b,c,r,u,n)

   INTEGER,INTENT(IN)      :: n
   REAL(r8),INTENT(IN)     :: a(n),b(n),c(n),r(n)
   REAL(r8),INTENT(INOUT)  :: u(n)
   !------------------------------
   !  Local variables
   !------------------------------
   INTEGER j
   REAL(r8) :: bet,gam(n)

   if(b(1).eq.0._r8) call endrun('ion_electron_temp: bt(1)=0 in tridag')
   bet=b(1)
   u(1)=r(1)/bet
   do j=2,n
      gam(j)=c(j-1)/bet
      bet=b(j)-a(j)*gam(j)
      if(bet.eq.0._r8) call endrun('ion_electron_temp: bet=0 in tridag')
      u(j)=(r(j)-a(j)*u(j-1))/bet
   end do

   do j=n-1,1,-1
      u(j)=u(j)-gam(j+1)*u(j+1)
   end do

   return

 END SUBROUTINE tridag

end module ion_electron_temp
