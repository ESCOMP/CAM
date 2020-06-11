! ====================================================================================
! Air-sea exchange module for trace gases
! Ref: Carpenter et al Chem Soc Rev (2012); Johnson, Ocean sci (2010)
! ------------------------------------------------------------------------------------
! Required inputs for the air-sea flux module:
!   - Seawater concentration (nanomoles per liter) and Sea surface salinity 
!     (parts per thousand) read from namelist (netCDF)
!   - Concentration in the gas-phase (pptv), air temperature (K), 10m windspeed (m/s),
!     surface pressure (atm), sea surface temperature (K): all from other modules
! ------------------------------------------------------------------------------------
! Key subroutines:
! ocean_emis_readnl(..):   Read salinity from namelist (user_nl_cam). 
!                          Salinity not time-dependent. Flux depends very weakly on it
! ocean_emis_init(...):    Interpolate salinity, initialize the library for the flux 
!                          reading time-dependent seawater conc. from user_nl_cam
! ocean_emis_advance(...): process the seawater concentration
! ocean_emis_getflux(...): calculate the air-sea flux (upward or downward), 
!                          then add to total surface flux (sflx)
! ------------------------------------------------------------------------------------
! Last built: 9 March 2018.
! Written by: Siyuan Wang (ACOM/NCAR)           siyuan@ucar.edu
! Acknowledgement: Francis Vitt (NCAR). and of course Dr. Peppurr too 
! ====================================================================================

module ocean_emis

  use shr_kind_mod,   only : r8=>shr_kind_r8, cl=>shr_kind_cl
  use ppgrid,         only : pcols, begchunk,endchunk
  use spmd_utils,     only : masterproc
  use cam_abortutils, only : endrun
  use cam_history,    only : addfld, horiz_only, outfld
  use constituents,   only : cnst_get_ind
  use tracer_data,    only : trfld,trfile
  use chem_mods,      only : gas_pcnst
  use cam_logfile,    only : iulog
  use ioFileMod,      only : getfil  

  implicit none

  private
  public :: ocean_emis_readnl
  public :: ocean_emis_init
  public :: ocean_emis_getflux
  public :: ocean_emis_advance

  type :: Csw
     integer                   :: spc_ndx
     real(r8)                  :: scalefactor
     character(len=256)        :: filename
     character(len=16)         :: species
     character(len=8)          :: units
     integer                   :: nsectors
     character(len=32),pointer :: sectors(:)
     type(trfld), pointer      :: fields(:)
     type(trfile)              :: file
  end type Csw

  logical                :: switch_bubble
  type(Csw), allocatable :: Csw_nM(:)
  integer                :: n_Csw_files 

  real(r8), allocatable :: salinity(:,:)      

  ! ================
  ! Air-sea exchange
  ! ================
  Integer, Parameter    :: HowManyMolecules = 16     ! Change this number if you wanna add more species
  Integer, Parameter    :: HowManyProperties = 16    ! Don't touch this (unless you wanna add more fields)
  Integer, Parameter    :: HowManySalts = 5          ! Change this number if you wanna add more salts
  Integer, Parameter    :: HowManySaltProperties = 7 ! Don't touch this (unless you wanna add more fields)

  Type GasLib                                                   
     Character(16)                          :: CmpdName
     Real(r8), Dimension(HowManyProperties) :: CmpdProperties
  End Type GasLib

  Type SaltLib                                                  
     Character(16)                              :: SaltName
     Real(r8), Dimension(HowManySaltProperties) :: SaltProperties 
  End Type SaltLib

  Type(GasLib), Dimension(HowManyMolecules) :: GasList  ! Library for the trace gas properties
  Type(SaltLib), Dimension(HowManySalts)    :: SaltList ! Library for the salt properties

  ! ===========================  
  ! seawater concentration:
  ! ===========================
  character(len=cl) :: csw_specifier(gas_pcnst) = ''                    
  character(len=24) :: csw_time_type = 'CYCLICAL' ! 'CYCLICAL' | 'SERIAL' | 'INTERP_MISSING_MONTHS'
  integer           :: csw_cycle_yr  = 0
  logical           :: bubble_mediated_transfer = .false.                                                                                                      
  character(len=cl) :: ocean_salinity_file = 'NONE'

contains

!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
  subroutine ocean_emis_readnl(nlfile)

    use namelist_utils, only : find_group_name
    use spmd_utils,     only : mpicom, masterprocid, mpi_character, mpi_integer, mpi_logical

    character(len=*), intent(in)  :: nlfile  ! filepath for file containing namelist input

    integer                       :: unitn, ierr
    character(len=*), parameter   :: subname = 'ocean_emis_readnl'

    ! ===================       
    ! Namelist definition
    ! ===================
    namelist /ocean_emis_nl/ ocean_salinity_file
    namelist /ocean_emis_nl/ csw_specifier, csw_time_type, csw_cycle_yr, bubble_mediated_transfer

    ! =============
    ! Read namelist
    ! =============
    if (masterproc) then
       open( newunit=unitn, file=trim(nlfile), status='old' )
       call find_group_name(unitn, 'ocean_emis_nl', status=ierr)
       if (ierr == 0) then
          read(unitn, ocean_emis_nl, iostat=ierr)
          if (ierr /= 0) then
             call endrun(subname // ':: ERROR reading namelist')
          end if
       end if
       close(unitn)
    end if
    
    ! ============================
    ! Broadcast namelist variables
    ! ============================
    call mpi_bcast(ocean_salinity_file, len(ocean_salinity_file),   mpi_character, masterprocid, mpicom, ierr)
    call mpi_bcast(csw_specifier,       cl*gas_pcnst,               mpi_character, masterprocid, mpicom, ierr)
    call mpi_bcast(csw_time_type,       len(csw_time_type),         mpi_character, masterprocid, mpicom, ierr)
    call mpi_bcast(csw_cycle_yr,        1,                          mpi_integer,   masterprocid, mpicom, ierr)
    call mpi_bcast(bubble_mediated_transfer, 1,                     mpi_logical,   masterprocid, mpicom, ierr)

    if (masterproc) then
      write(iulog,*) 'ocean_emis_readnl: ocean_salinity_file = '//trim(ocean_salinity_file)
      write(iulog,*) 'ocean_emis_readnl: bubble_mediated_transfer = ',bubble_mediated_transfer
    end if

  end subroutine ocean_emis_readnl

!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
  subroutine ocean_emis_init()

    use ioFileMod,        only : getfil
    use cam_pio_utils,    only : cam_pio_openfile
    use pio,              only : file_desc_t, pio_inq_dimid, pio_inq_dimlen, pio_inq_varid, pio_get_var
    use pio,              only : PIO_NOWRITE, PIO_NOERR
    use pio,              only : pio_seterrorhandling, PIO_BCAST_ERROR, pio_closefile
    use phys_grid,        only : get_ncols_p, get_rlon_all_p, get_rlat_all_p 
    use interpolate_data, only : lininterp_init, lininterp, interp_type, lininterp_finish
    use mo_constants,     only : pi

    integer               :: file_nlat, file_nlon
    integer               :: dimid, vid, ierr, c, ncols
    type(file_desc_t)     :: fid
    character(len=cl)     :: filen
    real(r8), allocatable :: file_lats(:), file_lons(:)
    real(r8), allocatable :: wrk2d(:,:)
    real(r8)              :: to_lats(pcols), to_lons(pcols)
    type(interp_type)     :: lon_wgts, lat_wgts 

    real(r8), parameter   :: zero=0_r8, twopi=2_r8*pi, degs2rads = pi/180._r8

    character(len=*), parameter :: subname = 'ocean_emis_init'
    
    if (trim(ocean_salinity_file) == 'NONE') return

    call getfil( ocean_salinity_file, filen, 0 )
    call cam_pio_openfile( fid, filen, PIO_NOWRITE)
    
    call pio_seterrorhandling(fid, PIO_BCAST_ERROR)
    
    ierr = pio_inq_dimid( fid, 'lon', dimid )
    if (ierr /= PIO_NOERR) then
       call endrun(subname//': pio_inq_dimid lon FAILED')
    endif
    ierr = pio_inq_dimlen( fid, dimid, file_nlon )
    if (ierr /= PIO_NOERR) then
       call endrun(subname//': pio_inq_dimlen file_nlon FAILED')
    endif

    allocate( file_lons(file_nlon) )
    ierr =  pio_inq_varid( fid, 'lon', vid )
    if (ierr /= PIO_NOERR) then
       call endrun(subname//': pio_inq_varid lon FAILED')
    endif
    ierr =  pio_get_var( fid, vid, file_lons )
    if (ierr /= PIO_NOERR) then
       call endrun(subname//': pio_get_var file_lons FAILED')
    endif
    file_lons = file_lons * degs2rads

    ierr = pio_inq_dimid( fid, 'lat', dimid )
    if (ierr /= PIO_NOERR) then
       call endrun(subname//': pio_inq_dimid lat FAILED')
    endif
    ierr = pio_inq_dimlen( fid, dimid, file_nlat )
    if (ierr /= PIO_NOERR) then
       call endrun(subname//': pio_inq_dimlen file_nlat FAILED')
    endif
    allocate( file_lats(file_nlat) )
    ierr =  pio_inq_varid( fid, 'lat', vid )
    if (ierr /= PIO_NOERR) then
       call endrun(subname//': pio_inq_varid lat FAILED')
    endif
    ierr =  pio_get_var( fid, vid, file_lats )
    if (ierr /= PIO_NOERR) then
       call endrun(subname//': pio_inq_varid SSS FAILED')
    endif
    file_lats = file_lats * degs2rads

    allocate( wrk2d( file_nlon, file_nlat ) )
    ierr =  pio_inq_varid( fid, 'SSS', vid )
    if (ierr /= PIO_NOERR) then
       call endrun(subname//': pio_inq_varid SSS FAILED')
    endif
    ierr = pio_get_var( fid, vid, wrk2d )
    if (ierr /= PIO_NOERR) then
       call endrun(subname//': pio_get_var SSS FAILED')
    endif

    allocate(salinity(pcols,begchunk:endchunk))

    do c=begchunk,endchunk

       ncols = get_ncols_p(c)
       call get_rlat_all_p(c, pcols, to_lats)
       call get_rlon_all_p(c, pcols, to_lons)

       call lininterp_init(file_lons, file_nlon, to_lons, ncols, 2, lon_wgts, zero, twopi)
       call lininterp_init(file_lats, file_nlat, to_lats, ncols, 1, lat_wgts)

       call lininterp(wrk2d, file_nlon, file_nlat, salinity(1:ncols,c), ncols, lon_wgts, lat_wgts)    

       call lininterp_finish(lon_wgts)
       call lininterp_finish(lat_wgts)

    end do

    deallocate( file_lons, file_lats )
    deallocate( wrk2d )

    call addfld('OCN_SALINITY', horiz_only, 'A', 'parts per thousands', 'ocean salinity' ) 

    ! ======================================================
    ! initializing the libraries for the air-sea flux module
    ! ======================================================
    Call CmpLibInitialization()
    Call SaltLibInitialization()

    ! ---------------------------------------------	
    ! Read seawater concentration: WSY
    ! ---------------------------------------------
    call cseawater_ini()

    call pio_closefile (fid)
    
  end subroutine ocean_emis_init

!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
  subroutine ocean_emis_advance( pbuf2d, state )
    ! -------------------------------
    ! check serial case for time span
    ! -------------------------------

    use physics_types,  only : physics_state
    use ppgrid,         only : begchunk, endchunk
    use tracer_data,    only : advance_trcdata
    use physics_buffer, only : physics_buffer_desc

    type(physics_state), intent(in)    :: state(begchunk:endchunk)                 
    type(physics_buffer_desc), pointer :: pbuf2d(:,:)
    integer                            :: m

    if (trim(ocean_salinity_file) == 'NONE') return

    do m = 1,n_Csw_files
       call advance_trcdata( Csw_nM(m)%fields, Csw_nM(m)%file, state, pbuf2d  )
    end do

  end subroutine ocean_emis_advance

!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
  subroutine ocean_emis_getflux(lchnk, ncol, state, u10, sst, ocnfrac, icefrac, sflx)

    use physics_types, only : physics_state
    use ppgrid,        only : pver                      

    integer, intent(in)                     :: lchnk, ncol
    type(physics_state), target, intent(in) :: state        ! Physics state variables
    real(r8), intent(in)                    :: u10(:)       ! Wind speed at 10m
    real(r8), intent(in)                    :: sst(:)       ! Sea surface temperature
    real(r8), intent(in)                    :: ocnfrac(:)   ! Ocean fraction
    real(r8), intent(in)                    :: icefrac(:)   ! Ice fraction
    real(r8), intent(inout)                 :: sflx(:,:)    ! Surface emissions (kg/m^2/s)

    integer             :: m, isec, SpeciesID
    real(r8)            :: Csw_col(ncol)
    real(r8)            :: MW_species
    real(r8) :: oceanflux_kg_m2_s(ncol)        

    if (trim(ocean_salinity_file) == 'NONE') return
    
    ! ==================================================
    ! Get seawater concentrations and calculate the flux
    ! ==================================================
    Csw_loop : do m = 1, n_Csw_files

       Csw_col(:) = 0._r8
       isec = 1
       Csw_col(:ncol) = Csw_nM(m)%scalefactor*Csw_nM(m)%fields(isec)%data(:ncol,1,lchnk)

       MW_species = MolecularWeight(SpeciesIndex( Csw_nM(m)%species )) 

       call cnst_get_ind( trim(Csw_nM(m)%species), SpeciesID, abort=.true. )

       oceanflux_kg_m2_s = 0.0_r8

       where (ocnfrac(:ncol) >= 0.2_r8 .and. Csw_col(:ncol) >= 0._r8)   ! calculate flux only for ocean
          oceanflux_kg_m2_s(:ncol) = Flux_kg_m2_s( &
               Csw_nM(m)%species,                  & ! name of species
               state%q(:ncol,pver,SpeciesID) * (28.97_r8/MW_species) * 1.0e+12_r8, & ! air concentration (ppt)
               Csw_col(:ncol),                     & ! sea water concentration (nM)
               state%t(:ncol,pver),                & ! air temperature (K)
               u10(:ncol),                         & ! wind speed at 10m (m/s) <- should use this
               state%ps(:ncol) / 101325.0_r8,      & ! surface pressure (atm)
               sst(:ncol),                         & ! sea surface temperautre (K)
               salinity(:ncol,lchnk),              & ! ocean salinity (parts per thousands)
               switch_bubble,                      & ! bubble-mediated transfer: on or off
               ncol  ) 
       end where

       ! ===========================================================================
       ! Add the ocean flux to the other fluxes 
       ! Make sure this ocean module is called after other surface emissions are set
       ! ===========================================================================
       sflx(:ncol,SpeciesID) = sflx(:ncol,SpeciesID) + oceanflux_kg_m2_s(:ncol) * ocnfrac(:ncol)

       ! ====================================
       ! Write stuff into the historial files
       ! ====================================
       call outfld('OCN_FLUX_' // trim(Csw_nM(m)%species), oceanflux_kg_m2_s(:ncol), ncol, lchnk)
       call outfld('Csw_' // trim(Csw_nM(m)%species), Csw_col(:ncol), ncol, lchnk)

    end do Csw_loop

    call outfld('OCN_SALINITY', salinity(:ncol,lchnk), ncol, lchnk)

  end subroutine ocean_emis_getflux

  
!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------

  Subroutine CmpLibInitialization()
    ! =====================================================================================
    ! This is the lookup table for molecular weight, Vb, and Henry's law constant
    ! NOTE: IF YOU ADDED NEW SPECIES, REMEMBER TO UPDATE THIS -> HowManyMolecules
    ! -------------------------------------------------------------------------------------
    ! Col 1:     molecular weight (g/mol)
    ! Col 2-10:  number of C, H, N, O, S, Br, Cl, F, and I atoms in the formula, respectively
    ! Col 11-13: number of double bound, triple bond, and rings in the molecule, respectively
    ! Col 14:    Vb (cm^3/mol) values if available.
    ! Col 15-16: Henry's law constant at 298 K (M/atm) and temperature dependency (in K^-1)
    ! -------------------------------------------------------------------------------------
    ! Ref: Henry's law constants from Sander ACP 2014 (compilation)
    ! =====================================================================================
    GasList(1) = GasLib('CH3OH',    (/ 32.05_r8, 1.0_r8, 4.0_r8, 0.0_r8, 1.0_r8, 0.0_r8, 0.0_r8, &
         0.0_r8, 0.0_r8, 0.0_r8, 0.0_r8, 0.0_r8, 0.0_r8, 42.7_r8, 203.0_r8, 5645.0_r8 /))
    GasList(2) = GasLib('C2H5OH',   (/ 46.07_r8, 2.0_r8, 6.0_r8, 0.0_r8, 1.0_r8, 0.0_r8, 0.0_r8, &
         0.0_r8, 0.0_r8, 0.0_r8, 0.0_r8, 0.0_r8, 0.0_r8, 0.0_r8, 190.0_r8, 6500.0_r8 /))
    GasList(3) = GasLib('CH2O',     (/ 30.03_r8, 1.0_r8, 2.0_r8, 0.0_r8, 1.0_r8, 0.0_r8, 0.0_r8, &
         0.0_r8, 0.0_r8, 0.0_r8, 1.0_r8, 0.0_r8, 0.0_r8, 0.0_r8, 3230.0_r8, 7100.0_r8 /))  
    GasList(4) = GasLib('CH3CHO',   (/ 44.05_r8, 2.0_r8, 4.0_r8, 0.0_r8, 1.0_r8, 0.0_r8, 0.0_r8, &
         0.0_r8, 0.0_r8, 0.0_r8, 1.0_r8, 0.0_r8, 0.0_r8, 0.0_r8, 12.9_r8, 5890.0_r8/))
    GasList(5) = GasLib('PROPANAL',  (/ 58.08_r8, 3.0_r8, 6.0_r8, 0.0_r8, 1.0_r8, 0.0_r8, 0.0_r8, &
         0.0_r8, 0.0_r8, 0.0_r8, 1.0_r8, 0.0_r8, 0.0_r8, 0.0_r8, 10.0_r8, 4330.0_r8 /))
    GasList(6) = GasLib('CH3COCH3', (/ 58.08_r8, 3.0_r8, 6.0_r8, 0.0_r8, 1.0_r8, 0.0_r8, 0.0_r8, &
         0.0_r8, 0.0_r8, 0.0_r8, 1.0_r8, 0.0_r8, 0.0_r8, 77.6_r8, 27.8_r8, 5530.0_r8 /))
    GasList(7) = GasLib('NO',       (/ 30.01_r8, 0.0_r8, 0.0_r8, 1.0_r8, 1.0_r8, 0.0_r8, 0.0_r8, &
         0.0_r8, 0.0_r8, 1.0_r8, 0.0_r8, 0.0_r8, 0.0_r8, 0.0_r8, 1.92e-3_r8, 1762.0_r8/)) ! NO
    GasList(8) = GasLib('HNO2',     (/ 47.01_r8, 0.0_r8, 1.0_r8, 1.0_r8, 2.0_r8, 0.0_r8, 0.0_r8, &
         0.0_r8, 0.0_r8, 0.0_r8, 1.0_r8, 0.0_r8, 0.0_r8, 0.0_r8, 49.0_r8, 4900.0_r8/))
    GasList(9) = GasLib('HCN',      (/ 27.03_r8, 1.0_r8, 1.0_r8, 1.0_r8, 0.0_r8, 0.0_r8, 0.0_r8, &
         0.0_r8, 0.0_r8, 0.0_r8, 0.0_r8, 1.0_r8, 0.0_r8, 0.0_r8, 9.02_r8, 8258.0_r8/))
    GasList(10) = GasLib('CH3CN',   (/ 41.05_r8, 2.0_r8, 3.0_r8, 1.0_r8, 0.0_r8, 0.0_r8, 0.0_r8, &
         0.0_r8, 0.0_r8, 0.0_r8, 0.0_r8, 1.0_r8, 0.0_r8, 55.5_r8, 52.8_r8, 3970.0_r8/))
    GasList(11) = GasLib('PAN',     (/ 121.05_r8, 2.0_r8, 3.0_r8, 1.0_r8, 5.0_r8, 0.0_r8, 0.0_r8, &
         0.0_r8, 0.0_r8, 0.0_r8, 1.0_r8, 0.0_r8, 0.0_r8, 0.0_r8, 2.8_r8, 5730.0_r8/))
    GasList(12) = GasLib('CH3BR',   (/ 94.94_r8, 1.0_r8, 3.0_r8, 0.0_r8, 0.0_r8, 0.0_r8, 1.0_r8, &
         0.0_r8, 0.0_r8, 0.0_r8, 0.0_r8, 0.0_r8, 0.0_r8, 0.0_r8, 0.16_r8, 3100.0_r8/))
    GasList(13) = GasLib('CH2BR2',  (/ 173.88_r8, 1.0_r8, 2.0_r8, 0.0_r8, 0.0_r8, 0.0_r8, 2.0_r8, &
         0.0_r8, 0.0_r8, 0.0_r8, 0.0_r8, 0.0_r8, 0.0_r8, 0.0_r8, 1.1_r8, 4000.0_r8/))
    GasList(14) = GasLib('CHBR3',   (/ 252.73_r8, 1.0_r8, 1.0_r8, 0.0_r8, 0.0_r8, 0.0_r8, 3.0_r8, &
         0.0_r8, 0.0_r8, 0.0_r8, 0.0_r8, 0.0_r8, 0.0_r8, 0.0_r8, 1.4_r8, 5000.0_r8/))
    GasList(15) = GasLib('DMS',     (/ 62.13_r8, 2.0_r8, 6.0_r8, 0.0_r8, 0.0_r8, 1.0_r8, 0.0_r8, &
         0.0_r8, 0.0_r8, 0.0_r8, 0.0_r8, 0.0_r8, 0.0_r8, 0.0_r8, 0.54_r8, 3460.0_r8/))
    GasList(16) = GasLib('CH3NO3',  (/ 74.04_r8, 1.0_r8, 3.0_r8, 1.0_r8, 3.0_r8, 0.0_r8, 0.0_r8, &
         0.0_r8, 0.0_r8, 0.0_r8, 0.0_r8, 0.0_r8, 0.0_r8, 0.0_r8, 2.0_r8, 4700.0_r8/))
    ! --------------------------------------------------------------------------------
    ! ADD NEW SPECIES HERE (molecular weight, formula info, and Henry's law constants)
    ! --------------------------------------------------------------------------------
  End Subroutine CmpLibInitialization

  Subroutine SaltLibInitialization()
    ! ================================================================================
    ! This is the lookup table for common solutes in seawater and the parameters to 
    ! calculate the dynamic viscosity of seawater. 
    ! You may add other solutes or change the mass fractions.
    ! --------------------------------------------------------------------------------
    ! Col 1:   mass fraction of solute
    ! Col 2-7: some parameter
    ! --------------------------------------------------------------------------------
    ! Ref: typical mass fractions of solute: Laliberte (2007)
    !      parameterizations: Laliberte (2007) and Millero et al (2008)
    ! ================================================================================
    SaltList(1) = SaltLib('NaCl',  (/ 0.798_r8, 16.22_r8,  1.3229_r8,   1.4849_r8, 0.0074691_r8,  30.78_r8,    2.0583_r8 /))
    SaltList(2) = SaltLib('KCl',   (/ 0.022_r8, 6.4883_r8, 1.3175_r8,  -0.7785_r8, 0.09272_r8,   -1.3_r8,      2.0811_r8 /))
    SaltList(3) = SaltLib('CaCl2', (/ 0.033_r8, 32.028_r8, 0.78792_r8, -1.1495_r8, 0.0026995_r8,  780860.0_r8, 5.8442_r8 /))
    SaltList(4) = SaltLib('MgCl2', (/ 0.047_r8, 24.032_r8, 2.2694_r8,   3.7108_r8, 0.021853_r8,  -1.1236_r8,   0.14474_r8/))
    SaltList(5) = SaltLib('MgSO4', (/ 0.100_r8, 72.269_r8, 2.2238_r8,   6.6037_r8, 0.0079004_r8,  3340.1_r8,   6.1304_r8 /))
    ! ---------------------------------------------
    ! ADD NEW SALT HERE (mass fraction, and others)
    ! ---------------------------------------------
  End Subroutine SaltLibInitialization

  Function SpeciesIndex(SpeciesName)
    ! ==============================================
    ! This function is to look for the species index
    ! ==============================================
    Integer             :: SpeciesIndex, i
    Character(Len=16)   :: SpeciesName

    SpeciesIndex = -1 ! return -1 if species is not found
 
    Do i = 1, HowManyMolecules
       If (trim(SpeciesName) == trim(GasList(i)%CmpdName)) Then
          SpeciesIndex = i
          Exit
       Endif
    End Do
  End Function SpeciesIndex

  Function Flux_kg_m2_s(SpeciesName,Cgas_ppt,Cwater_nM,T_air_K,u10_m_s,P_atm,T_water_K,&
                        Salinity_PartsPerThousand,switch_bubble,ncol)
    ! ===========================================================================
    ! This is the main module function. Input variables:
    ! ---------------------------------------------------------------------------
    !    - SpeciesName: name of species
    !    - Cgas_ppt: mixing ratio (parts per trillion) of trace gas of interest 
    !      in the gas-phase (lowest modeling layer)
    !    - Cwater_nM: concentration of trace gas of interest in the surface ocean
    !    - T_air_K: temperature in the lowest modeling layer
    !    - u10_m_s: wind speed at 10 m above sea surface level
    !    - P_atm: air pressure in atm at sea surface level
    !    - T_water_K: sea surface temperature
    !    - Salinity_PartsPerThousand: surface ocean salinity
    !    - switch_bubble: bubble-mediated transfer switch
    ! All must be 1D arrays with same dimension(ncol, so CESM-compatible)
    ! ===========================================================================
    Integer                   :: ncol, SpeciesID
    Character(16)             :: SpeciesName
    Real(r8), Dimension(ncol) :: Flux_kg_m2_s
    Real(r8), Dimension(ncol) :: Cgas_ppt, Cwater_nM, T_air_K, u10_m_s, P_atm, T_water_K, Salinity_PartsPerThousand
    Real(r8), Dimension(ncol) :: H_gas_over_liquid_dimless, kt_m_s
    Logical                   :: switch_bubble

    where(Salinity_PartsPerThousand .lt. 0.0_r8) Salinity_PartsPerThousand = 33.0_r8

    SpeciesID = SpeciesIndex(SpeciesName)       
    H_gas_over_liquid_dimless = 1.0_r8/(Henry_M_atm(SpeciesID,T_water_K,Salinity_PartsPerThousand,ncol)*&
                                        0.082_r8*T_water_K)
    If (switch_bubble) then
       ! --------------------------------------------------------
       ! k_water parameterization with bubble-induced enhancement
       ! --------------------------------------------------------
       kt_m_s = (1.0_r8/k_water_m_s_bubble(SpeciesID, T_water_K, Salinity_PartsPerThousand, &
                                           u10_m_s, Cgas_ppt, P_atm, T_air_K, ncol) &
               + 1.0_r8/k_air_m_s(SpeciesID, u10_m_s, T_air_K, P_atm, ncol)&
                       /H_gas_over_liquid_dimless)**(-1.0_r8)
    else
       ! ------------------------------------------------
       ! Original k_water parameterization, scaled to CO2
       ! ------------------------------------------------
       kt_m_s = (1.0_r8/k_water_m_s(SpeciesID, T_water_K, Salinity_PartsPerThousand, u10_m_s, ncol) &
            + 1.0_r8/k_air_m_s(SpeciesID, u10_m_s, T_air_K, P_atm, ncol)/H_gas_over_liquid_dimless)**(-1.0_r8)
    endif
    Flux_kg_m2_s = kt_m_s * (Cwater_nM*1E-9_r8*1000.0_r8                                               &
         - Cgas_ppt*1E-12_r8*(101325.0_r8*P_atm)/8.314_r8/T_air_K/H_gas_over_liquid_dimless) & ! g/m2/s
         * MolecularWeight(SpeciesIndex(SpeciesName)) / 1000.0_r8                              ! convert to kg/m2/s
  End Function Flux_kg_m2_s


  Function k_air_m_s(SpeciesIndex, u10_m_s, T_air_K, P_atm, ncol)
    use shr_const_mod, only: vonKarman=>SHR_CONST_KARMAN
    ! =============================================================================
    ! Air-side transfer velocity. Slightly modified NOAA COARE (Fairall et al 2003; 
    ! Feffery et al 2010), as recommended by Johnson Ocean Sci. 2010. 
    ! Dynamic viscosity of air: Tsilingiris 2008
    ! =============================================================================
    Integer                   :: ncol, SpeciesIndex
    Real(r8), Dimension(ncol) :: k_air_m_s 
    Real(r8), Dimension(ncol) :: u10_m_s, T_air_K, P_atm, ustar_m_s, DragCoeff
    Real(r8), Dimension(ncol) :: DynamicViscosityAir_kg_m_s, DensityAir_kg_m3, DiffusivityInAir, SchmidtNumberInAir

    ! WSY: If local friction velocity is available from the model, might as well use that?
    ustar_m_s = u10_m_s * sqrt(6.1E-4_r8 + 6.3E-5_r8 * u10_m_s)
    DragCoeff = (ustar_m_s / u10_m_s)**2.0_r8
    DynamicViscosityAir_kg_m_s = 1.715747771E-5_r8 + 4.722402075E-8_r8 * (T_air_K-273.15_r8) &
         - 3.663027156E-10_r8 * ((T_air_K-273.15_r8)**2.0_r8) &
         + 1.873236686E-12_r8 * ((T_air_K-273.15_r8)**3.0_r8) &
         - 8.050218737E-14_r8 * ((T_air_K-273.15_r8)**4.0_r8)    
    DensityAir_kg_m3 = 1.293393662_r8 - 5.538444326e-3_r8 * (T_air_K-273.15_r8) &
         + 3.860201577e-5_r8 * (T_air_K-273.15_r8)**2.0_r8 &
         - 5.2536065e-7_r8 * (T_air_K-273.15_r8)**3.0_r8
    DiffusivityInAir = DiffusivityInAir_cm2_s(SpeciesIndex, T_air_K, P_atm, ncol)       
    SchmidtNumberInAir = DynamicViscosityAir_kg_m_s / DensityAir_kg_m3 / (DiffusivityInAir/10000.0_r8)  
    k_air_m_s = 1E-3_r8 + ustar_m_s / (13.3_r8*(SchmidtNumberInAir**0.5_r8)+(DragCoeff**(-0.5_r8))-&
                5.0_r8+log(SchmidtNumberInAir)/2.0_r8/vonKarman)
  End Function k_air_m_s




  Function k_water_m_s(SpeciesIndex, T_water_K, Salinity_PartsPerThousand, u10_m_s, ncol)
    ! ================================================================================
    ! Water-side transfer velocity. Ref: Nightingale et al (2000). Salinity considered
    ! ================================================================================
    Integer                    :: ncol, SpeciesIndex
    Real(r8), Dimension(ncol)  :: k_water_m_s
    Real(r8), Dimension(ncol)  :: T_water_K, Salinity_PartsPerThousand, u10_m_s
    Real(r8), Dimension(ncol)  :: DiffusivityInWater, SchmidtNumberInWater
    Real(r8)                   :: SchmidtNumberInWater_CO2ref
    SchmidtNumberInWater_CO2ref = 660.0_r8              ! this is the Schmidt number of CO2 at 20 degC in fresh water
    DiffusivityInWater = DiffusivityInWater_cm2_s(SpeciesIndex, T_water_K, Salinity_PartsPerThousand, ncol)
    SchmidtNumberInWater = DynamicViscosityWater_g_m_s(T_water_K, Salinity_PartsPerThousand, ncol) / 1000.0_r8 &
         / DensityWater_kg_m3(T_water_K,Salinity_PartsPerThousand,ncol)/(DiffusivityInWater/10000.0_r8)
    k_water_m_s = ((0.222_r8*(u10_m_s**2.0_r8)+0.333_r8*u10_m_s)*&
                  ((SchmidtNumberInWater/SchmidtNumberInWater_CO2ref)**(-0.5_r8)))/360000.0_r8
  End Function k_water_m_s




  Function k_water_m_s_bubble(SpeciesIndex, T_water_K, Salinity_PartsPerThousand, u10_m_s, Cgas_ppt, P_atm, T_air_K, ncol)
    ! ==============================================================
    ! Water-side transfer velocity. Ref: Asher and Wanninkhof (1998).
    ! ==============================================================
    Integer                   :: ncol, SpeciesIndex
    Real(r8), Dimension(ncol) :: k_water_m_s_bubble
    Real(r8), Dimension(ncol) :: T_water_K, Salinity_PartsPerThousand, u10_m_s, Cgas_ppt, P_atm, T_air_K
    Real(r8), Dimension(ncol) :: DiffusivityInWater, SchmidtNumberInWater
    Real(r8), Dimension(ncol) :: FracCoverage_WhiteCaps, OstwaldSolubilityCoefficient
    DiffusivityInWater = DiffusivityInWater_cm2_s(SpeciesIndex, T_water_K, Salinity_PartsPerThousand, ncol)
    SchmidtNumberInWater = DynamicViscosityWater_g_m_s(T_water_K, Salinity_PartsPerThousand, ncol) / 1000.0_r8 &
         / DensityWater_kg_m3(T_water_K,Salinity_PartsPerThousand,ncol)/(DiffusivityInWater/10000.0_r8)
    FracCoverage_WhiteCaps = 2.56e-6_r8 * (u10_m_s - 1.77_r8)**3.0_r8      
    OstwaldSolubilityCoefficient = Henry_M_atm(SpeciesIndex,T_water_K,Salinity_PartsPerThousand,ncol) ! just Henry's law (M/atm)
    OstwaldSolubilityCoefficient = OstwaldSolubilityCoefficient * (Cgas_ppt*1.0E-12_r8*P_atm)         ! mol / L
    OstwaldSolubilityCoefficient = OstwaldSolubilityCoefficient * 0.082_r8 * T_air_K / P_atm          ! L / L
    k_water_m_s_bubble = ((47.0_r8*u10_m_s + FracCoverage_WhiteCaps*(115200.0_r8 - 47.0_r8* u10_m_s)) &
                       *(SchmidtNumberInWater**(-0.5_r8))  &
                       + FracCoverage_WhiteCaps * (-37.0_r8/OstwaldSolubilityCoefficient &
                       + 6120.0_r8*(OstwaldSolubilityCoefficient**(-0.37_r8)) *(SchmidtNumberInWater**(-0.18_r8)))) &
                       * 2.8e-6_r8
  End Function k_water_m_s_bubble




  Function DiffusivityInAir_cm2_s(SpeciesIndex, T_air_K, P_atm, ncol)
    ! ============================
    ! Ref: Johnson Ocean Sci. 2010
    ! ============================
    Integer                   :: ncol, SpeciesIndex
    Real(r8), Dimension(ncol) :: DiffusivityInAir_cm2_s, T_air_K, P_atm
    Real(r8), parameter       :: MW_air = 28.97_r8   ! molecular weight for air
    Real(r8), parameter       :: Va = 20.1_r8        ! molar volume for air
    Real(r8)                  :: Vb, MW_species
    Vb = LiquidMolarVolume_cm3_mol(SpeciesIndex)
    MW_species = MolecularWeight(SpeciesIndex)
    DiffusivityInAir_cm2_s = 0.001_r8 * (T_air_K**1.75_r8) &    ! oh f* me
         * (((MW_air + MW_species)/(MW_air*MW_species))**0.5_r8) &
         / ((P_atm*(Va**(1.0_r8/3.0_r8)+Vb**(1.0_r8/3.0_r8)))**2.0_r8)
  End Function DiffusivityInAir_cm2_s



  Function DiffusivityInWater_cm2_s(SpeciesIndex, T_water_K, Salinity_PartsPerThousand, ncol)
    ! =================================================
    ! Ref: Johnson Ocean Sci. 2010. Salinity considered
    ! =================================================
    Integer                   :: ncol, SpeciesIndex
    Real(r8), Dimension(ncol) :: DiffusivityInWater_cm2_s, DynamicViscosityWater, T_water_K, Salinity_PartsPerThousand
    Real(r8), parameter       :: AssociationFactor = 2.6_r8     ! ... for water
    Real(r8)                  :: Vb, MW_species
    Vb = LiquidMolarVolume_cm3_mol(SpeciesIndex)
    MW_species = MolecularWeight(SpeciesIndex)
    DynamicViscosityWater = DynamicViscosityWater_g_m_s(T_water_K, Salinity_PartsPerThousand, ncol)
    ! -------------------------------------------------
    ! Wilke and Chang 1955: this seems to be a bit high
    ! -------------------------------------------------
    ! DiffusivityInWater_cm2_s = 7.4E-8_r8 * T_water_K * sqrt(AssociationFactor*MW_species) /
    !                            DynamicViscosityWater / (Vb**0.6_r8)
    ! ----------------------
    ! Hayduk and Minhas 1982
    ! ----------------------
    DiffusivityInWater_cm2_s = 1.25E-8_r8 * (T_water_K**1.52_r8) * (DynamicViscosityWater**(9.58_r8/Vb - 1.12_r8)) &
                             * (Vb**(-0.19_r8)-0.292_r8)

  End Function DiffusivityInWater_cm2_s


  Function DynamicViscosityWater_g_m_s(T_water_K, Salinity_PartsPerThousand, ncol)
    ! =================================================
    ! Ref: Johnson Ocean Sci. 2010. Salinity considered
    ! =================================================
    Integer                    :: ncol
    Real(r8), Dimension(ncol)  :: DynamicViscosityWater_g_m_s, T_water_K, Salinity_PartsPerThousand 
    Real(r8), Dimension(ncol)  :: MassFrac_water, DynamicViscosityPureWater_g_m_s, SaltViscosity, sum_w_ln_SaltViscosity
    Integer                    :: j, n 
    sum_w_ln_SaltViscosity = 0.0_r8
    MassFrac_water = 1.0_r8 - Salinity_PartsPerThousand / 1000.0_r8
    DynamicViscosityPureWater_g_m_s = ((T_water_K-273.15_r8)+246.0_r8) &
         / (0.05594_r8*(T_water_K-273.15_r8)**2.0_r8+5.2842_r8*(T_water_K-273.15_r8)+137.37_r8)   
    Do j = 1, ncol      
       If (Salinity_PartsPerThousand(j) == 0.0_r8) Then          ! pure water
          DynamicViscosityWater_g_m_s(j) = DynamicViscosityPureWater_g_m_s(j)
       Else                                                      ! salty water
          Do n = 1, HowManySalts
             SaltViscosity(j) = exp((SaltList(n)%SaltProperties(2) * &
                  (Salinity_PartsPerThousand(j)/1000.0_r8)**SaltList(n)%SaltProperties(3) &
                  + SaltList(n)%SaltProperties(4)) &
                  / (SaltList(n)%SaltProperties(5)*(T_water_K(j)-273.15_r8) + 1.0_r8)) &
                  / (SaltList(n)%SaltProperties(6) * (Salinity_PartsPerThousand(j) / &
                  1000.0_r8)**SaltList(n)%SaltProperties(7) + 1.0_r8)
             sum_w_ln_SaltViscosity(j) = sum_w_ln_SaltViscosity(j) + (Salinity_PartsPerThousand(j)/1000.0_r8) &
                  * SaltList(n)%SaltProperties(1) * log(SaltViscosity(j))
          End Do
          DynamicViscosityWater_g_m_s(j) = exp(MassFrac_water(j) &
               * log(DynamicViscosityPureWater_g_m_s(j)) + sum_w_ln_SaltViscosity(j))
       Endif
    End Do
  End Function DynamicViscosityWater_g_m_s


  Function DensityWater_kg_m3(T_water_K, Salinity_PartsPerThousand, ncol)
    ! ====================================================
    ! Ref: Millero and Poisson (1981). Salinity considered
    ! ====================================================
    Integer                   :: ncol
    Real(r8), Dimension(ncol) :: DensityWater_kg_m3, T_water_K, Salinity_PartsPerThousand
    Real(r8), Dimension(ncol) :: DensityPureWater_kg_m3, FactorA, FactorB, FactorC
    DensityPureWater_kg_m3 = 999.842594_r8 + 0.06793952_r8*(T_water_K-273.15_r8) &
         - 0.00909529_r8*((T_water_K-273.15_r8)**2.0_r8) &
         + 0.0001001685_r8*((T_water_K-273.15_r8)**3.0_r8) &
         - 0.000001120083_r8*((T_water_K-273.15_r8)**4.0_r8) &
         + 0.000000006536332_r8*((T_water_K-273.15_r8)**5.0_r8)
    FactorA = 0.824493_r8 - 0.0040899_r8*(T_water_K-273.15_r8) + 0.000076438_r8*((T_water_K-273.15_r8)**2.0_r8) &
         - 0.00000082467_r8*((T_water_K-273.15_r8)**3.0_r8) + 0.0000000053875_r8*((T_water_K-273.15_r8)**4.0_r8)
    FactorB = - 0.00572466_r8 + 0.00010277_r8*(T_water_K-273.15_r8) - 0.0000016546_r8*((T_water_K-273.15_r8)**2.0_r8)
    FactorC = 0.00048314_r8
    DensityWater_kg_m3 = DensityPureWater_kg_m3 + FactorA*Salinity_PartsPerThousand &
         + FactorB*(Salinity_PartsPerThousand**(2.0_r8/3.0_r8)) + FactorC*Salinity_PartsPerThousand
  End Function DensityWater_kg_m3


  Function Henry_M_atm(SpeciesIndex, T_water_K, Salinity_PartsPerThousand, ncol)
    ! =========================================================================================
    ! Ref: Sander compilation 2015. Salt-in or salt-out estimated based on Setschenow constants
    ! =========================================================================================
    Integer                    :: ncol, j
    Integer                    :: SpeciesIndex
    Real(r8), Dimension(ncol)  :: Henry_M_atm, T_water_K, Salinity_PartsPerThousand
    Real(r8), Dimension(ncol)  :: Heff_M_atm_PureWater, Setschenow, Heff_M_atm_SaltyWater
    Heff_M_atm_PureWater = GasList(SpeciesIndex)%CmpdProperties(15) * &
         exp(GasList(SpeciesIndex)%CmpdProperties(16) * (1.0_r8/T_water_K - 1.0_r8/298.0_r8))
    Do j = 1, ncol      
       If (Salinity_PartsPerThousand(j)==0.0_r8) Then
          Henry_M_atm(j) = Heff_M_atm_PureWater(j) 
       Else
          Setschenow(j) = log(LiquidMolarVolume_cm3_mol(SpeciesIndex)) * &
               (7.33532E-4_r8 + 3.39615E-5_r8 * log(Heff_M_atm_PureWater(j)) &
               - 2.40888E-6_r8 * ((log(Heff_M_atm_PureWater(j)))**2.0_r8) &
               + 1.57114E-7_r8 * ((log(Heff_M_atm_PureWater(j)))**3.0_r8))
          Heff_M_atm_SaltyWater(j) = Heff_M_atm_PureWater(j) * 10.0_r8**(Setschenow(j)*Salinity_PartsPerThousand(j))
          Henry_M_atm(j) = Heff_M_atm_SaltyWater(j)
       Endif
    End Do
  End Function Henry_M_atm


  Function MolecularWeight(SpeciesIndex)
    Real(r8)  :: MolecularWeight
    Integer   :: SpeciesIndex
    MolecularWeight = GasList(SpeciesIndex)%CmpdProperties(1)
  End Function MolecularWeight


  Function LiquidMolarVolume_cm3_mol(SpeciesIndex)
    ! ===========================================================================
    ! If no measurements available, i.e. GasList(SpeciesIndex)%CmpdProperties(14)
    ! is zero, Schroeder's group contribution method is used
    ! ===========================================================================
    Real(r8)    :: LiquidMolarVolume_cm3_mol
    Integer     :: SpeciesIndex

    If (GasList(SpeciesIndex)%CmpdProperties(14)/=0.0_r8) Then  
       LiquidMolarVolume_cm3_mol = GasList(SpeciesIndex)%CmpdProperties(14)
    Else
       LiquidMolarVolume_cm3_mol = 7.0_r8*GasList(SpeciesIndex)%CmpdProperties(2)                                ! C
       LiquidMolarVolume_cm3_mol = LiquidMolarVolume_cm3_mol + 7.0_r8*GasList(SpeciesIndex)%CmpdProperties(3)    ! H
       LiquidMolarVolume_cm3_mol = LiquidMolarVolume_cm3_mol + 7.0_r8*GasList(SpeciesIndex)%CmpdProperties(4)    ! N
       LiquidMolarVolume_cm3_mol = LiquidMolarVolume_cm3_mol + 7.0_r8*GasList(SpeciesIndex)%CmpdProperties(5)    ! O
       LiquidMolarVolume_cm3_mol = LiquidMolarVolume_cm3_mol + 21.0_r8*GasList(SpeciesIndex)%CmpdProperties(6)   ! S
       LiquidMolarVolume_cm3_mol = LiquidMolarVolume_cm3_mol + 31.5_r8*GasList(SpeciesIndex)%CmpdProperties(7)   ! Br
       LiquidMolarVolume_cm3_mol = LiquidMolarVolume_cm3_mol + 24.5_r8*GasList(SpeciesIndex)%CmpdProperties(8)   ! Cl
       LiquidMolarVolume_cm3_mol = LiquidMolarVolume_cm3_mol + 10.5_r8*GasList(SpeciesIndex)%CmpdProperties(9)   ! F
       LiquidMolarVolume_cm3_mol = LiquidMolarVolume_cm3_mol + 38.5_r8*GasList(SpeciesIndex)%CmpdProperties(10)  ! I
       LiquidMolarVolume_cm3_mol = LiquidMolarVolume_cm3_mol - 7.0_r8*GasList(SpeciesIndex)%CmpdProperties(13)   ! Ring
       LiquidMolarVolume_cm3_mol = LiquidMolarVolume_cm3_mol + 7.0_r8*GasList(SpeciesIndex)%CmpdProperties(11)   ! Double bond
       LiquidMolarVolume_cm3_mol = LiquidMolarVolume_cm3_mol + 14.0_r8*GasList(SpeciesIndex)%CmpdProperties(12)  ! Triple bond
    Endif

  End Function LiquidMolarVolume_cm3_mol

  subroutine cseawater_ini()

    use mo_chem_utls,     only : get_spc_ndx 
    use tracer_data,      only : trcdata_init   
    use cam_pio_utils,    only : cam_pio_openfile       
    use pio,              only : pio_inquire, pio_nowrite, pio_closefile, pio_inq_varndims
    use pio,              only : pio_inq_varname, file_desc_t, pio_get_att, PIO_NOERR, PIO_GLOBAL
    use pio,              only : pio_seterrorhandling, PIO_BCAST_ERROR    
    use string_utils,     only : GLC    

    integer             :: i, j, l, m, n, nn, astat, vid, ierr, nvars, isec
    integer             :: indx(gas_pcnst)      
    type(file_desc_t)   :: ncid
    character(len=16)   :: csw_species(gas_pcnst)
    character(len=256)  :: csw_filenam(gas_pcnst)
    integer             :: csw_indexes(gas_pcnst)
    real(r8)            :: csw_scalefactor(gas_pcnst)

    character(len=80)   :: file_interp_type = ' '

    character(len=16)   :: spc_name
    character(len=256)  :: filename
    character(len=256)  :: tmp_string = ' '
    character(len=32)   :: xchr = ' '
    real(r8)            :: xdbl

    character(len=32)              :: varname
    character(len=256)             :: locfn
    integer, allocatable           :: vndims(:)
    character(len=1),   parameter  :: filelist = ''
    character(len=1),   parameter  :: datapath = ''
    logical,            parameter  :: rmv_file = .false.

    character(len=*), parameter :: subname = 'cseawater_ini'

    ! ========================================================  
    ! Read sea water concentration specifier from the namelist
    ! ========================================================

    nn = 0
    indx(:) = 0

    switch_bubble = bubble_mediated_transfer

    count_Csw: do n=1,gas_pcnst
       if ( len_trim(csw_specifier(n) ) == 0 ) then
          exit count_Csw
       endif

       i = scan(csw_specifier(n),'->')
       spc_name = trim(adjustl(csw_specifier(n)(:i-1)))

       ! --------------------------
       ! get the scaling factor
       ! --------------------------
       tmp_string = adjustl(csw_specifier(n)(i+2:))
       j = scan( tmp_string, '*' )
       if (j>0) then
          xchr = tmp_string(1:j-1)               ! get the multipler (left of the '*')
          read( xchr, * ) xdbl                   ! convert the string to a real
          tmp_string = adjustl(tmp_string(j+1:)) ! get the filepath name (right of the '*')
       else
          xdbl = 1._r8
       endif
       filename = trim(tmp_string)

       m = get_spc_ndx(spc_name)

       if ( m<1 ) then
          if (masterproc) write(iulog,*) 'cseawater_inti: spc_name ',spc_name,' is not included in the simulation'
          call endrun('cseawater_inti: invalid seawater concentration specification')
       endif

       nn = nn+1
       csw_species(nn) = spc_name
       csw_filenam(nn) = filename
       csw_indexes(nn) = m
       csw_scalefactor(nn) = xdbl

       indx(n)=n

    enddo count_Csw

    n_Csw_files = nn

    if (masterproc) write(iulog,*) subname//': n_Csw_files = ',n_Csw_files

    allocate( Csw_nM(n_Csw_files), stat=astat )
    if( astat/= 0 ) then
       write(iulog,*) subname//': failed to allocate Csw_nM array; error = ',astat
       call endrun(subname//': failed to allocate Csw_nM array')
    end if

    ! -------------------------------------------
    ! Setup the seawater concentration type array
    ! -------------------------------------------
    do m=1,n_Csw_files 
       Csw_nM(m)%spc_ndx = csw_indexes(indx(m))
       Csw_nM(m)%units = 'nM'
       Csw_nM(m)%species = csw_species(indx(m))
       Csw_nM(m)%filename = csw_filenam(indx(m))
       Csw_nM(m)%scalefactor = csw_scalefactor(indx(m))
    enddo

    !-----------------------------------------------
    ! read emis files to determine number of sectors
    !-----------------------------------------------
    spc_loop: do m = 1, n_Csw_files

       Csw_nM(m)%nsectors = 0

       call getfil (Csw_nM(m)%filename, locfn, 0)
       call cam_pio_openfile ( ncid, trim(locfn), PIO_NOWRITE)

       call pio_seterrorhandling(ncid, PIO_BCAST_ERROR)

       ierr = pio_inquire(ncid, nvariables=nvars)
       if (ierr /= PIO_NOERR) then
          call endrun(subname//': pio_inquire nvars FAILED')
       endif

       allocate(vndims(nvars))

       do vid = 1,nvars

          ierr = pio_inq_varndims (ncid, vid, vndims(vid))
          if (ierr /= PIO_NOERR) then
             call endrun(subname//': pio_inq_varndims FAILED')
          endif

          if( vndims(vid) < 3 ) then
             cycle
          elseif( vndims(vid) > 3 ) then
             ierr = pio_inq_varname (ncid, vid, varname)
             if (ierr /= PIO_NOERR) then
                call endrun(subname//': pio_inq_varname varname FAILED')
             endif
             write(iulog,*) subname//': Skipping variable ', trim(varname),', ndims = ',vndims(vid), &
                            ', species=',trim(Csw_nM(m)%species)
             cycle
          end if

          Csw_nM(m)%nsectors = Csw_nM(m)%nsectors+1

       enddo

       allocate( Csw_nM(m)%sectors(Csw_nM(m)%nsectors), stat=astat )
       if( astat/= 0 ) then
          write(iulog,*) subname//': failed to allocate Csw_nM(m)%sectors array; error = ',astat
          call endrun(subname//': failed to allocate Csw_nM(m)%sectors array')
       end if

       isec = 1

       do vid = 1,nvars
          if( vndims(vid) == 3 ) then
             ierr = pio_inq_varname(ncid, vid, Csw_nM(m)%sectors(isec))
             if (ierr /= PIO_NOERR) then
                call endrun(subname//': pio_inq_varname Csw_nM(m)%sectors(isec) FAILED')
             endif
             isec = isec+1
          endif

       enddo
       deallocate(vndims)

       ! Global attribute 'input_method' overrides the srf_emis_type namelist setting on
       ! a file-by-file basis.  If the emis file does not contain the 'input_method' 
       ! attribute then the srf_emis_type namelist setting is used.
       ierr = pio_get_att(ncid, PIO_GLOBAL, 'input_method', file_interp_type)   
       if ( ierr == PIO_NOERR) then
          l = GLC(file_interp_type)
          csw_time_type(1:l) = file_interp_type(1:l)
          csw_time_type(l+1:) = ' '
       endif

       call pio_closefile (ncid)

       allocate(Csw_nM(m)%file%in_pbuf(size(Csw_nM(m)%sectors)))
       Csw_nM(m)%file%in_pbuf(:) = .false.

       call trcdata_init( Csw_nM(m)%sectors, &
            Csw_nM(m)%filename, filelist, datapath, &
            Csw_nM(m)%fields,  &
            Csw_nM(m)%file, &
            rmv_file, csw_cycle_yr, 0, 0, trim(csw_time_type) )

    enddo spc_loop

    ! ===================================
    ! Output stuff to the history files
    ! ===================================
    do m = 1, n_Csw_files
       call addfld('OCN_FLUX_' // trim(Csw_nM(m)%species), horiz_only, 'A', 'kg/m2/s', &
            'ocean flux ' // trim(Csw_nM(m)%species)  )
       call addfld('Csw_' // trim(Csw_nM(m)%species), horiz_only, 'A', 'nanomole per liter (nM)', &
            'seeawater concentration ' // trim(Csw_nM(m)%species)  )
    end do

  end subroutine cseawater_ini


end module ocean_emis
