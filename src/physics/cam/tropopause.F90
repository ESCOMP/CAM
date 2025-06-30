! This is the CAM interface to the CCPP-ized tropopause_find scheme.
! Full compatibility, bit-for-bit, to old CAM approach is achieved through
! this module, however this module will not be necessary in CAM-SIMA.
!
! For science description of the underlying algorithms, refer to
! atmospheric_physics/tropopause_find/tropopause_find.F90.
! (hplin, 8/20/24)

module tropopause
  !---------------------------------------------------------------
  ! ... variables for the tropopause module
  !---------------------------------------------------------------

  use shr_kind_mod,         only : r8 => shr_kind_r8
  use shr_const_mod,        only : pi => shr_const_pi
  use ppgrid,               only : pcols, pver, pverp, begchunk, endchunk
  use cam_abortutils,       only : endrun
  use cam_logfile,          only : iulog
  use cam_history_support,  only : fillvalue
  use physics_types,        only : physics_state
  use spmd_utils,           only : masterproc

  implicit none

  private

  public  :: tropopause_readnl, tropopause_init, tropopause_find_cam, tropopause_output
  public  :: tropopause_findChemTrop
  public  :: TROP_ALG_NONE, TROP_ALG_ANALYTIC, TROP_ALG_CLIMATE
  public  :: TROP_ALG_STOBIE, TROP_ALG_HYBSTOB, TROP_ALG_TWMO, TROP_ALG_WMO
  public  :: TROP_ALG_CPP
  public  :: NOTFOUND

  save

  ! These parameters define and enumeration to be used to define the primary
  ! and backup algorithms to be used with the tropopause_find() method. The
  ! backup algorithm is meant to provide a solution when the primary algorithm
  ! fail. The algorithms that can't fail are: TROP_ALG_ANALYTIC, TROP_ALG_CLIMATE
  ! and TROP_ALG_STOBIE.
  integer, parameter    :: TROP_ALG_NONE      = 1    ! Don't evaluate
  integer, parameter    :: TROP_ALG_ANALYTIC  = 2    ! Analytic Expression
  integer, parameter    :: TROP_ALG_CLIMATE   = 3    ! Climatology
  integer, parameter    :: TROP_ALG_STOBIE    = 4    ! Stobie Algorithm
  integer, parameter    :: TROP_ALG_TWMO      = 5    ! WMO Definition, Reichler et al. [2003]
  integer, parameter    :: TROP_ALG_WMO       = 6    ! WMO Definition
  integer, parameter    :: TROP_ALG_HYBSTOB   = 7    ! Hybrid Stobie Algorithm
  integer, parameter    :: TROP_ALG_CPP       = 8    ! Cold Point Parabolic
  integer, parameter    :: TROP_ALG_CHEMTROP  = 9    ! Chemical tropopause

  ! Note: exclude CHEMTROP here as it is a new flag added in CCPP-ized routines to unify the chemTrop routine. (hplin, 8/20/24)
  integer, parameter    :: TROP_NALG          = 8    ! Number of Algorithms
  character,parameter   :: TROP_LETTER(TROP_NALG) = (/ ' ', 'A', 'C', 'S', 'T', 'W', 'H', 'F' /)
                                                     ! unique identifier for output, don't use P

  ! These variables should probably be controlled by namelist entries.
  logical ,parameter    :: output_all         = .False.              ! output tropopause info from all algorithms
  integer ,parameter    :: default_primary    = TROP_ALG_TWMO        ! default primary algorithm
  integer ,parameter    :: default_backup     = TROP_ALG_CLIMATE     ! default backup algorithm

  ! Namelist variables
  character(len=256)    :: tropopause_climo_file = 'trop_climo'      ! absolute filepath of climatology file

  ! These variables are used to store the climatology data.
  real(r8)              :: days(12)                                  ! days in the climatology
  real(r8), pointer     :: tropp_p_loc(:,:,:)                        ! climatological tropopause pressures

  integer, parameter :: NOTFOUND = -1

  ! physical constants
  ! These constants are set in module variables rather than as parameters
  ! to support the aquaplanet mode in which the constants have values determined
  ! by the experiment protocol
  real(r8) :: cnst_kap     ! = cappa
  real(r8) :: cnst_faktor  ! = -gravit/rair
  real(r8) :: cnst_ka1     ! = cnst_kap - 1._r8

!================================================================================================
contains
!================================================================================================

   ! Read namelist variables.
   subroutine tropopause_readnl(nlfile)

      use namelist_utils,  only: find_group_name
      use units,           only: getunit, freeunit
      use mpishorthand

      character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

      ! Local variables
      integer :: unitn, ierr
      character(len=*), parameter :: subname = 'tropopause_readnl'

      namelist /tropopause_nl/ tropopause_climo_file
      !-----------------------------------------------------------------------------

      if (masterproc) then
         unitn = getunit()
         open( unitn, file=trim(nlfile), status='old' )
         call find_group_name(unitn, 'tropopause_nl', status=ierr)
         if (ierr == 0) then
            read(unitn, tropopause_nl, iostat=ierr)
            if (ierr /= 0) then
               call endrun(subname // ':: ERROR reading namelist')
            end if
         end if
         close(unitn)
         call freeunit(unitn)
      end if

#ifdef SPMD
      ! Broadcast namelist variables
      call mpibcast(tropopause_climo_file, len(tropopause_climo_file), mpichar, 0, mpicom)
#endif

   end subroutine tropopause_readnl


  ! This routine is called during intialization and must be called before the
  ! other methods in this module can be used. Its main tasks are to read in the
  ! climatology from a file and to define the output fields. Much of this code
  ! is taken from mo_tropopause.
  subroutine tropopause_init()

    use cam_history,     only: addfld, horiz_only
    use tropopause_find, only: tropopause_find_init
    use physconst,       only: cappa, rair, gravit, pi

    character(len=512) :: errmsg
    integer            :: errflg

    ! Call underlying CCPP-initialization routine.
    call tropopause_find_init(cappa, rair, gravit, pi, errmsg, errflg)

    ! Define the output fields.
    call addfld('TROP_P',          horiz_only,  'A', 'Pa',          'Tropopause Pressure',              flag_xyfill=.True.)
    call addfld('TROP_T',          horiz_only,  'A', 'K',           'Tropopause Temperature',           flag_xyfill=.True.)
    call addfld('TROP_Z',          horiz_only,  'A', 'm',           'Tropopause Height',                flag_xyfill=.True.)
    call addfld('TROP_DZ',         (/ 'lev' /), 'A', 'm',           'Relative Tropopause Height')
    call addfld('TROP_PD',         (/ 'lev' /), 'A', 'probability', 'Tropopause Probabilty')
    call addfld('TROP_FD',         horiz_only,  'A', 'probability', 'Tropopause Found')

    call addfld('TROPP_P',         horiz_only,  'A', 'Pa',          'Tropopause Pressure (primary)',    flag_xyfill=.True.)
    call addfld('TROPP_T',         horiz_only,  'A', 'K',           'Tropopause Temperature (primary)', flag_xyfill=.True.)
    call addfld('TROPP_Z',         horiz_only,  'A', 'm',           'Tropopause Height (primary)',      flag_xyfill=.True.)
    call addfld('TROPP_DZ',        (/ 'lev' /), 'A', 'm',           'Relative Tropopause Height (primary)')
    call addfld('TROPP_PD',        (/ 'lev' /), 'A', 'probability', 'Tropopause Distribution (primary)')
    call addfld('TROPP_FD',        horiz_only,  'A', 'probability', 'Tropopause Found (primary)')

    call addfld('TROPF_P',         horiz_only,  'A',  'Pa',         'Tropopause Pressure (cold point)',    flag_xyfill=.True.)
    call addfld('TROPF_T',         horiz_only,  'A',  'K',          'Tropopause Temperature (cold point)', flag_xyfill=.True.)
    call addfld('TROPF_Z',         horiz_only,  'A',  'm',          'Tropopause Height (cold point)',      flag_xyfill=.True.)
    call addfld('TROPF_DZ',        (/ 'lev' /),  'A', 'm',          'Relative Tropopause Height (cold point)', flag_xyfill=.True.)
    call addfld('TROPF_PD',        (/ 'lev' /), 'A', 'probability', 'Tropopause Distribution (cold point)')
    call addfld('TROPF_FD',        horiz_only,  'A', 'probability', 'Tropopause Found (cold point)')

    call addfld( 'hstobie_trop',   (/ 'lev' /), 'I',  'fraction of model time', 'Lowest level with stratospheric chemsitry')
    call addfld( 'hstobie_linoz',  (/ 'lev' /), 'I',  'fraction of model time', 'Lowest possible Linoz level')
    call addfld( 'hstobie_tropop', (/ 'lev' /), 'I', 'fraction of model time', &
         'Troposphere boundary calculated in chemistry' )

    ! If requested, be prepared to output results from all of the methods.
    if (output_all) then
      call addfld('TROPA_P',  horiz_only,  'A',  'Pa',          'Tropopause Pressure (analytic)',        flag_xyfill=.True.)
      call addfld('TROPA_T',  horiz_only,  'A',  'K',          'Tropopause Temperature (analytic)',      flag_xyfill=.True.)
      call addfld('TROPA_Z',  horiz_only,  'A',  'm',          'Tropopause Height (analytic)',           flag_xyfill=.True.)
      call addfld('TROPA_PD', (/ 'lev' /), 'A', 'probability', 'Tropopause Distribution (analytic)')
      call addfld('TROPA_FD', horiz_only,  'A', 'probability', 'Tropopause Found (analytic)')

      call addfld('TROPC_P',  horiz_only,  'A',  'Pa',         'Tropopause Pressure (climatology)',      flag_xyfill=.True.)
      call addfld('TROPC_T',  horiz_only,  'A',  'K',          'Tropopause Temperature (climatology)',   flag_xyfill=.True.)
      call addfld('TROPC_Z',  horiz_only,  'A',  'm',          'Tropopause Height (climatology)',        flag_xyfill=.True.)
      call addfld('TROPC_PD', (/ 'lev' /), 'A', 'probability', 'Tropopause Distribution (climatology)')
      call addfld('TROPC_FD', horiz_only,  'A', 'probability', 'Tropopause Found (climatology)')

      call addfld('TROPS_P',  horiz_only,  'A',  'Pa',         'Tropopause Pressure (stobie)',           flag_xyfill=.True.)
      call addfld('TROPS_T',  horiz_only,  'A',  'K',          'Tropopause Temperature (stobie)',        flag_xyfill=.True.)
      call addfld('TROPS_Z',  horiz_only,  'A',  'm',          'Tropopause Height (stobie)',             flag_xyfill=.True.)
      call addfld('TROPS_PD', (/ 'lev' /), 'A', 'probability', 'Tropopause Distribution (stobie)')
      call addfld('TROPS_FD', horiz_only,  'A', 'probability', 'Tropopause Found (stobie)')

      call addfld('TROPT_P',  horiz_only,  'A',  'Pa',         'Tropopause Pressure (twmo)',             flag_xyfill=.True.)
      call addfld('TROPT_T',  horiz_only,  'A',  'K',          'Tropopause Temperature (twmo)',          flag_xyfill=.True.)
      call addfld('TROPT_Z',  horiz_only,  'A',  'm',          'Tropopause Height (twmo)',               flag_xyfill=.True.)
      call addfld('TROPT_PD', (/ 'lev' /), 'A', 'probability', 'Tropopause Distribution (twmo)')
      call addfld('TROPT_FD', horiz_only,  'A', 'probability', 'Tropopause Found (twmo)')

      call addfld('TROPW_P',  horiz_only,  'A',  'Pa',         'Tropopause Pressure (WMO)',              flag_xyfill=.True.)
      call addfld('TROPW_T',  horiz_only,  'A',  'K',          'Tropopause Temperature (WMO)',           flag_xyfill=.True.)
      call addfld('TROPW_Z',  horiz_only,  'A',  'm',          'Tropopause Height (WMO)',                flag_xyfill=.True.)
      call addfld('TROPW_PD', (/ 'lev' /), 'A', 'probability', 'Tropopause Distribution (WMO)')
      call addfld('TROPW_FD', horiz_only,  'A', 'probability', 'Tropopause Found (WMO)')

      call addfld('TROPH_P',  horiz_only,  'A',  'Pa',         'Tropopause Pressure (Hybrid Stobie)',    flag_xyfill=.True.)
      call addfld('TROPH_T',  horiz_only,  'A',  'K',          'Tropopause Temperature (Hybrid Stobie)', flag_xyfill=.True.)
      call addfld('TROPH_Z',  horiz_only,  'A',  'm',          'Tropopause Height (Hybrid Stobie)',      flag_xyfill=.True.)
      call addfld('TROPH_PD', (/ 'lev' /), 'A', 'probability', 'Tropopause Distribution (Hybrid Stobie)')
      call addfld('TROPH_FD', horiz_only,  'A', 'probability', 'Tropopause Found (Hybrid Stobie)')
    end if


    call tropopause_read_file()


  end subroutine tropopause_init


  subroutine tropopause_read_file
    !------------------------------------------------------------------
    ! ... initialize upper boundary values
    !------------------------------------------------------------------
    use interpolate_data,  only : lininterp_init, lininterp, interp_type, lininterp_finish
    use dyn_grid,     only : get_dyn_grid_parm
    use phys_grid,    only : get_ncols_p, get_rlat_all_p, get_rlon_all_p
    use ioFileMod,    only : getfil
    use time_manager, only : get_calday
    use physconst,    only : pi
    use cam_pio_utils, only: cam_pio_openfile
    use pio,          only : file_desc_t, var_desc_t, pio_inq_dimid, pio_inq_dimlen, &
         pio_inq_varid, pio_get_var, pio_closefile, pio_nowrite

    !------------------------------------------------------------------
    ! ... local variables
    !------------------------------------------------------------------
    integer :: i, j, n
    integer :: ierr
    type(file_desc_t) :: pio_id
    integer :: dimid
    type(var_desc_t) :: vid
    integer :: nlon, nlat, ntimes
    integer :: start(3)
    integer :: count(3)
    integer, parameter :: dates(12) = (/ 116, 214, 316, 415,  516,  615, &
         716, 816, 915, 1016, 1115, 1216 /)
    integer :: plon, plat
    type(interp_type) :: lon_wgts, lat_wgts
    real(r8), allocatable :: tropp_p_in(:,:,:)
    real(r8), allocatable :: lat(:)
    real(r8), allocatable :: lon(:)
    real(r8) :: to_lats(pcols), to_lons(pcols)
    real(r8), parameter :: d2r=pi/180._r8, zero=0._r8, twopi=pi*2._r8
    character(len=256) :: locfn
    integer  :: c, ncols


    plon = get_dyn_grid_parm('plon')
    plat = get_dyn_grid_parm('plat')


    !-----------------------------------------------------------------------
    !       ... open netcdf file
    !-----------------------------------------------------------------------
    call getfil (tropopause_climo_file, locfn, 0)
    call cam_pio_openfile(pio_id, trim(locfn), PIO_NOWRITE)

    !-----------------------------------------------------------------------
    !       ... get time dimension
    !-----------------------------------------------------------------------
    ierr = pio_inq_dimid( pio_id, 'time', dimid )
    ierr = pio_inq_dimlen( pio_id, dimid, ntimes )
    if( ntimes /= 12 )then
       write(iulog,*) 'tropopause_init: number of months = ',ntimes,'; expecting 12'
       call endrun
    end if
    !-----------------------------------------------------------------------
    !       ... get latitudes
    !-----------------------------------------------------------------------
    ierr = pio_inq_dimid( pio_id, 'lat', dimid )
    ierr = pio_inq_dimlen( pio_id, dimid, nlat )
    allocate( lat(nlat), stat=ierr )
    if( ierr /= 0 ) then
       write(iulog,*) 'tropopause_init: lat allocation error = ',ierr
       call endrun
    end if
    ierr = pio_inq_varid( pio_id, 'lat', vid )
    ierr = pio_get_var( pio_id, vid, lat )
    lat(:nlat) = lat(:nlat) * d2r
    !-----------------------------------------------------------------------
    !       ... get longitudes
    !-----------------------------------------------------------------------
    ierr = pio_inq_dimid( pio_id, 'lon', dimid )
    ierr = pio_inq_dimlen( pio_id, dimid, nlon )
    allocate( lon(nlon), stat=ierr )
    if( ierr /= 0 ) then
       write(iulog,*) 'tropopause_init: lon allocation error = ',ierr
       call endrun
    end if
    ierr = pio_inq_varid( pio_id, 'lon', vid )
    ierr = pio_get_var( pio_id, vid, lon )
    lon(:nlon) = lon(:nlon) * d2r

    !------------------------------------------------------------------
    !  ... allocate arrays
    !------------------------------------------------------------------
    allocate( tropp_p_in(nlon,nlat,ntimes), stat=ierr )
    if( ierr /= 0 ) then
       write(iulog,*) 'tropopause_init: tropp_p_in allocation error = ',ierr
       call endrun
    end if
    !------------------------------------------------------------------
    !  ... read in the tropopause pressure
    !------------------------------------------------------------------
    ierr = pio_inq_varid( pio_id, 'trop_p', vid )
    start = (/ 1, 1, 1 /)
    count = (/ nlon, nlat, ntimes /)
    ierr = pio_get_var( pio_id, vid, start, count, tropp_p_in )

    !------------------------------------------------------------------
    !  ... close the netcdf file
    !------------------------------------------------------------------
    call pio_closefile( pio_id )

    !--------------------------------------------------------------------
    !  ... regrid
    !--------------------------------------------------------------------

    allocate( tropp_p_loc(pcols,begchunk:endchunk,ntimes), stat=ierr )

    if( ierr /= 0 ) then
      write(iulog,*) 'tropopause_init: tropp_p_loc allocation error = ',ierr
      call endrun
    end if

    do c=begchunk,endchunk
       ncols = get_ncols_p(c)
       call get_rlat_all_p(c, pcols, to_lats)
       call get_rlon_all_p(c, pcols, to_lons)
       call lininterp_init(lon, nlon, to_lons, ncols, 2, lon_wgts, zero, twopi)
       call lininterp_init(lat, nlat, to_lats, ncols, 1, lat_wgts)
       do n=1,ntimes
          call lininterp(tropp_p_in(:,:,n), nlon, nlat, tropp_p_loc(1:ncols,c,n), ncols, lon_wgts, lat_wgts)
       end do
       call lininterp_finish(lon_wgts)
       call lininterp_finish(lat_wgts)
    end do
    deallocate(lon)
    deallocate(lat)
    deallocate(tropp_p_in)

    !--------------------------------------------------------
    ! ... initialize the monthly day of year times
    !--------------------------------------------------------

    do n = 1,12
       days(n) = get_calday( dates(n), 0 )
    end do
    if (masterproc) then
       write(iulog,*) 'tropopause_init : days'
       write(iulog,'(1p,5g15.8)') days(:)
    endif

  end subroutine tropopause_read_file

  ! Searches all the columns in the chunk and attempts to identify the tropopause.
  ! Two routines can be specifed, a primary routine which is tried first and a
  ! backup routine which will be tried only if the first routine fails. If the
  ! tropopause can not be identified by either routine, then a NOTFOUND is returned
  ! for the tropopause level, temperature and pressure.
  subroutine tropopause_find_cam(pstate, tropLev, tropP, tropT, tropZ, primary, backup)

    use tropopause_find, only: tropopause_findWithBackup

    use cam_history,     only: outfld
    use time_manager,    only: get_curr_calday

    implicit none

    type(physics_state), intent(in)     :: pstate
    integer, optional, intent(in)       :: primary                   ! primary detection algorithm
    integer, optional, intent(in)       :: backup                    ! backup detection algorithm
    integer,            intent(out)     :: tropLev(:)                ! tropopause level index
    real(r8), optional, intent(out)     :: tropP(:)                  ! tropopause pressure (Pa)
    real(r8), optional, intent(out)     :: tropT(:)                  ! tropopause temperature (K)
    real(r8), optional, intent(out)     :: tropZ(:)                  ! tropopause height (m)

    ! Local Variable
    integer       :: primAlg            ! Primary algorithm
    integer       :: backAlg            ! Backup algorithm

    real(r8)      :: calday
    integer       :: ncol

    real(r8)      :: hstobie_trop  (pcols, pver)
    real(r8)      :: hstobie_linoz (pcols, pver)
    real(r8)      :: hstobie_tropop(pcols, pver)

    character(len=512) :: errmsg
    integer            :: errflg

    ! Get compatibility variables for CCPP-ized routine
    ncol   = pstate%ncol
    calday = get_curr_calday()

    ! Initialize the results to a missing value, so that the algorithms will
    ! attempt to find the tropopause for all of them. Only do this for the active columns.
    tropLev(:ncol) = NOTFOUND
    if (present(tropP)) tropP(:ncol) = fillvalue
    if (present(tropT)) tropT(:ncol) = fillvalue
    if (present(tropZ)) tropZ(:ncol) = fillvalue

    ! Set the algorithms to be used, either the ones provided or the defaults.
    if (present(primary)) then
      primAlg = primary
    else
      primAlg = default_primary
    end if

    if (present(backup)) then
      backAlg = backup
    else
      backAlg = default_backup
    end if

    ! This does not call the tropopause_find_run routine directly, because it
    ! computes multiple needed tropopauses simultaneously. Instead, here we
    ! specify the algorithm needed directly to the algorithm driver routine.
    call tropopause_findWithBackup( &
         ncol           = ncol, &
         pver           = pver, &
         fillvalue      = fillvalue, &
         lat            = pstate%lat(:ncol), &
         pint           = pstate%pint(:ncol, :pverp), &
         pmid           = pstate%pmid(:ncol, :pver), &
         t              = pstate%t(:ncol, :pver), &
         zi             = pstate%zi(:ncol, :pverp), &
         zm             = pstate%zm(:ncol, :pver), &
         phis           = pstate%phis(:ncol), &
         calday         = calday, &
         tropp_p_loc    = tropp_p_loc(:ncol,pstate%lchnk,:), &  ! Subset into chunk as the underlying routines are no longer chunkized.
         tropp_days     = days, &
         tropLev        = tropLev(:ncol), &
         tropP          = tropP, &
         tropT          = tropT, &
         tropZ          = tropZ, &
         primary        = primAlg, &
         backup         = backAlg, &
         hstobie_trop   = hstobie_trop(:ncol, :pver), &    ! Only used if TROP_ALG_HYBSTOB
         hstobie_linoz  = hstobie_linoz(:ncol, :pver), &   ! Only used if TROP_ALG_HYBSTOB
         hstobie_tropop = hstobie_tropop(:ncol, :pver), &  ! Only used if TROP_ALG_HYBSTOB
         errmsg         = errmsg, &
         errflg         = errflg &
    )

    ! Output hybridstobie specific fields
    if(primAlg == TROP_ALG_HYBSTOB) then
       call outfld('hstobie_trop',   hstobie_trop(:ncol,:),    ncol, pstate%lchnk )
       call outfld('hstobie_linoz',  hstobie_linoz(:ncol,:),   ncol, pstate%lchnk )
       call outfld('hstobie_tropop', hstobie_tropop(:ncol,:),  ncol, pstate%lchnk )
    endif
  end subroutine tropopause_find_cam

  ! Searches all the columns in the chunk and attempts to identify the "chemical"
  ! tropopause. This is the lapse rate tropopause, backed up by the climatology
  ! if the lapse rate fails to find the tropopause at pressures higher than a certain
  ! threshold. This pressure threshold depends on latitude. Between 50S and 50N,
  ! the climatology is used if the lapse rate tropopause is not found at P > 75 hPa.
  ! At high latitude (poleward of 50), the threshold is increased to 125 hPa to
  ! eliminate false events that are sometimes detected in the cold polar stratosphere.
  !
  ! NOTE: This routine was adapted from code in chemistry.F90 and mo_gasphase_chemdr.F90.
  subroutine tropopause_findChemTrop(pstate, tropLev)

    use tropopause_find, only: tropopause_findWithBackup

    use time_manager,    only: get_curr_calday

    implicit none

    type(physics_state), intent(in)     :: pstate
    integer,             intent(out)    :: tropLev(:)            ! tropopause level index

    ! Local Variable
    real(r8)            :: calday
    integer             :: i
    integer             :: ncol

    character(len=512) :: errmsg
    integer            :: errflg

    ! Get compatibility variables for CCPP-ized routine
    ncol   = pstate%ncol
    calday = get_curr_calday()

    ! Now call the unified routine with the CHEMTROP option, which has automatic
    ! backup fall to climatology.
    call tropopause_findWithBackup( &
         ncol           = ncol, &
         pver           = pver, &
         fillvalue      = fillvalue, &
         lat            = pstate%lat(:ncol), &
         pint           = pstate%pint(:ncol, :pverp), &
         pmid           = pstate%pmid(:ncol, :pver), &
         t              = pstate%t(:ncol, :pver), &
         zi             = pstate%zi(:ncol, :pverp), &
         zm             = pstate%zm(:ncol, :pver), &
         phis           = pstate%phis(:ncol), &
         calday         = calday, &
         tropp_p_loc    = tropp_p_loc(:ncol,pstate%lchnk,:), &  ! Subset into chunk as the underlying routines are no longer chunkized.
         tropp_days     = days, &
         tropLev        = tropLev(1:ncol), &
         primary        = TROP_ALG_CHEMTROP, &
         backup         = TROP_ALG_CLIMATE, &
         errmsg         = errmsg, &
         errflg         = errflg &
    )
  end subroutine tropopause_findChemTrop

  ! Output the tropopause pressure and temperature to the history files. Two sets
  ! of output will be generated, one for the default algorithm and another one
  ! using the default routine, but backed by a climatology when the default
  ! algorithm fails.
  subroutine tropopause_output(pstate)
    use cam_history,  only : outfld

    implicit none

    type(physics_state), intent(in)     :: pstate

    ! Local Variables
    integer       :: i
    integer       :: alg
    integer       :: ncol                     ! number of cloumns in the chunk
    integer       :: lchnk                    ! chunk identifier
    integer       :: tropLev(pcols)           ! tropopause level index
    real(r8)      :: tropP(pcols)             ! tropopause pressure (Pa)
    real(r8)      :: tropT(pcols)             ! tropopause temperature (K)
    real(r8)      :: tropZ(pcols)             ! tropopause height (m)
    real(r8)      :: tropFound(pcols)         ! tropopause found
    real(r8)      :: tropDZ(pcols, pver)      ! relative tropopause height (m)
    real(r8)      :: tropPdf(pcols, pver)     ! tropopause probability distribution

    ! Information about the chunk.
    lchnk = pstate%lchnk
    ncol  = pstate%ncol

    ! Find the tropopause using the default algorithm backed by the climatology.
    call tropopause_find_cam(pstate, tropLev, tropP=tropP, tropT=tropT, tropZ=tropZ)

    tropPdf(:,:) = 0._r8
    tropFound(:) = 0._r8
    tropDZ(:,:) = fillvalue
    do i = 1, ncol
      if (tropLev(i) /= NOTFOUND) then
        tropPdf(i, tropLev(i)) = 1._r8
        tropFound(i) = 1._r8
        tropDZ(i,:) = pstate%zm(i,:) - tropZ(i)
      end if
    end do

    call outfld('TROP_P',   tropP(:ncol),      ncol, lchnk)
    call outfld('TROP_T',   tropT(:ncol),      ncol, lchnk)
    call outfld('TROP_Z',   tropZ(:ncol),      ncol, lchnk)
    call outfld('TROP_DZ',  tropDZ(:ncol, :), ncol, lchnk)
    call outfld('TROP_PD',  tropPdf(:ncol, :), ncol, lchnk)
    call outfld('TROP_FD',  tropFound(:ncol),  ncol, lchnk)


    ! Find the tropopause using just the primary algorithm.
    call tropopause_find_cam(pstate, tropLev, tropP=tropP, tropT=tropT, tropZ=tropZ, backup=TROP_ALG_NONE)

    tropPdf(:,:) = 0._r8
    tropFound(:) = 0._r8
    tropDZ(:,:) = fillvalue

    do i = 1, ncol
      if (tropLev(i) /= NOTFOUND) then
        tropPdf(i, tropLev(i)) = 1._r8
        tropFound(i) = 1._r8
        tropDZ(i,:) = pstate%zm(i,:) - tropZ(i)
      end if
    end do

    call outfld('TROPP_P',   tropP(:ncol),      ncol, lchnk)
    call outfld('TROPP_T',   tropT(:ncol),      ncol, lchnk)
    call outfld('TROPP_Z',   tropZ(:ncol),      ncol, lchnk)
    call outfld('TROPP_DZ',  tropDZ(:ncol, :), ncol, lchnk)
    call outfld('TROPP_PD',  tropPdf(:ncol, :), ncol, lchnk)
    call outfld('TROPP_FD',  tropFound(:ncol),  ncol, lchnk)


    ! Find the tropopause using just the cold point algorithm.
    call tropopause_find_cam(pstate, tropLev, tropP=tropP, tropT=tropT, tropZ=tropZ, primary=TROP_ALG_CPP, backup=TROP_ALG_NONE)

    tropPdf(:,:) = 0._r8
    tropFound(:) = 0._r8
    tropDZ(:,:) = fillvalue

    do i = 1, ncol
      if (tropLev(i) /= NOTFOUND) then
        tropPdf(i, tropLev(i)) = 1._r8
        tropFound(i) = 1._r8
        tropDZ(i,:) = pstate%zm(i,:) - tropZ(i)
      end if
    end do

    call outfld('TROPF_P',   tropP(:ncol),      ncol, lchnk)
    call outfld('TROPF_T',   tropT(:ncol),      ncol, lchnk)
    call outfld('TROPF_Z',   tropZ(:ncol),      ncol, lchnk)
    call outfld('TROPF_DZ',  tropDZ(:ncol, :), ncol, lchnk)
    call outfld('TROPF_PD',  tropPdf(:ncol, :), ncol, lchnk)
    call outfld('TROPF_FD',  tropFound(:ncol),  ncol, lchnk)


    ! If requested, do all of the algorithms.
    if (output_all) then

      do alg = 2, TROP_NALG

        ! Find the tropopause using just the analytic algorithm.
        call tropopause_find_cam(pstate, tropLev, tropP=tropP, tropT=tropT, tropZ=tropZ, primary=alg, backup=TROP_ALG_NONE)

        tropPdf(:,:) = 0._r8
        tropFound(:) = 0._r8

        do i = 1, ncol
          if (tropLev(i) /= NOTFOUND) then
            tropPdf(i, tropLev(i)) = 1._r8
            tropFound(i) = 1._r8
          end if
        end do

        call outfld('TROP' // TROP_LETTER(alg) // '_P',   tropP(:ncol),      ncol, lchnk)
        call outfld('TROP' // TROP_LETTER(alg) // '_T',   tropT(:ncol),      ncol, lchnk)
        call outfld('TROP' // TROP_LETTER(alg) // '_Z',   tropZ(:ncol),      ncol, lchnk)
        call outfld('TROP' // TROP_LETTER(alg) // '_PD',  tropPdf(:ncol, :), ncol, lchnk)
        call outfld('TROP' // TROP_LETTER(alg) // '_FD',  tropFound(:ncol),  ncol, lchnk)
      end do
    end if
  end subroutine tropopause_output
end module tropopause
