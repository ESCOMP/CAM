module scamMod
  !----------------------------------------------------------------------
  !
  ! This module provides control variables and namelist functionality
  ! for SCAM.
  !
  ! As with global CAM, SCAM is initialized with state information
  ! from the initial and boundary files. For each succeeding timestep
  ! SCAM calculates the physics tendencies and combines these with
  ! the advective tendencies provided by the IOP file to produce
  ! a forcast. Generally, the control variables in this module
  ! determine what data and parameterizations will be used to make
  ! the forecast. For instance, many of the control variables in
  ! this module provide flexibility to affect the forecast by overriding
  ! parameterization prognosed tendencies with observed tendencies
  ! of a particular field program recorded on the IOP file.
  !
  ! Public functions/subroutines:
  !   scam_readnl
  !-----------------------------------------------------------------------

use shr_kind_mod,   only: r8 => shr_kind_r8, cl => shr_kind_cl
use spmd_utils,     only: masterproc,npes
use pmgrid,         only: plon, plat, plev, plevp
use constituents,   only: cnst_get_ind, pcnst, cnst_name
use netcdf,         only: NF90_NOERR,NF90_CLOSE,NF90_GET_VAR,NF90_INQUIRE_DIMENSION, &
                          NF90_INQ_DIMID, NF90_INQ_VARID, NF90_NOWRITE, NF90_OPEN, &
                          NF90_GET_ATT,NF90_GLOBAL,NF90_INQUIRE_ATTRIBUTE, &
                          NF90_INQUIRE_VARIABLE, NF90_MAX_VAR_DIMS, nf90_get_var
use shr_scam_mod,   only: shr_scam_getCloseLatLon
use cam_logfile,    only: iulog
use cam_abortutils, only: endrun
use time_manager,   only: get_curr_date, get_nstep,is_first_step,get_start_date,timemgr_time_inc
use error_messages, only: handle_ncerr


implicit none
private

! PUBLIC INTERFACES:

public :: scam_readnl         ! read SCAM namelist options
public :: readiopdata         ! read iop boundary data
public :: setiopupdate        ! find index in iopboundary data for current time
public :: plevs0              ! Define the pressures of the interfaces and midpoints
public :: scmiop_flbc_inti
public :: setiopupdate_init

! PUBLIC MODULE DATA:

real(r8), public ::  pressure_levels(plev)
real(r8), public ::  scmlat   ! input namelist latitude for scam
real(r8), public ::  scmlon   ! input namelist longitude for scam
real(r8), public ::  closeioplat   ! closest iop latitude for scam
real(r8), public ::  closeioplon   ! closest iop longitude for scam
integer,  public ::  closeioplatidx   ! file array index of closest iop latitude for scam
integer,  public ::  closeioplonidx   ! file array index closest iop longitude for scam


integer, parameter :: num_switches = 20
integer, parameter :: max_path_len = 128

logical, public ::  single_column         ! Using IOP file or not
logical, public ::  use_iop               ! Using IOP file or not
logical, public ::  use_pert_init         ! perturb initial values
logical, public ::  use_pert_frc          ! perturb forcing
logical, public ::  switch(num_switches)  ! Logical flag settings from GUI
logical, public ::  l_uvphys              ! If true, update u/v after TPHYS
logical, public ::  l_uvadvect            ! If true, T, U & V will be passed to SLT
logical, public ::  l_conv                ! use flux divergence terms for T and q?
logical, public ::  l_divtr               ! use flux divergence terms for constituents?
logical, public ::  l_diag                ! do we want available diagnostics?

integer, public ::  error_code            ! Error code from netCDF reads
integer, public ::  initTimeIdx
integer, public ::  seedval
integer :: bdate, last_date, last_sec

character(len=max_path_len), public ::  modelfile
character(len=max_path_len), public ::  analysisfile
character(len=max_path_len), public ::  sicfile
character(len=max_path_len), public ::  userfile
character(len=max_path_len), public ::  sstfile
character(len=max_path_len), public ::  lsmpftfile
character(len=max_path_len), public ::  pressfile
character(len=max_path_len), public ::  topofile
character(len=max_path_len), public ::  ozonefile
character(len=max_path_len), public ::  iopfile
character(len=max_path_len), public ::  absemsfile
character(len=max_path_len), public ::  aermassfile
character(len=max_path_len), public ::  aeropticsfile
character(len=max_path_len), public ::  timeinvfile
character(len=max_path_len), public ::  lsmsurffile
character(len=max_path_len), public ::  lsminifile

! note that scm_zadv_q is set to slt to be consistent with CAM BFB testing


character(len=16), public    :: scm_zadv_T  = 'eulc            '
character(len=16), public    :: scm_zadv_q  = 'slt             '
character(len=16), public    :: scm_zadv_uv = 'eulc            '

real(r8), public ::  fixmascam
real(r8), public ::  betacam
real(r8), public ::  alphacam(pcnst)
real(r8), public ::  dqfxcam(plon,plev,pcnst)

real(r8), public ::      divq3d(plev,pcnst)  ! 3D q advection
real(r8), public ::      divt3d(plev)        ! 3D T advection
real(r8), public ::      divu3d(plev)        ! 3D U advection
real(r8), public ::      divv3d(plev)        ! 3D V advection
real(r8), public ::      vertdivq(plev,pcnst)! vertical q advection
real(r8), public ::      vertdivt(plev)      ! vertical T advection
real(r8), public ::      vertdivu(plev)      ! vertical T advection
real(r8), public ::      vertdivv(plev)      ! vertical T advection
real(r8), public ::      ptend               ! surface pressure tendency
real(r8), public ::      qdiff(plev)         ! model minus observed humidity
real(r8), public ::      qobs(plev)          ! actual W.V. Mixing ratio
real(r8), public ::      qinitobs(plev,pcnst)! initial tracer field
real(r8), public ::      cldliqobs(plev)     ! actual W.V. Mixing ratio
real(r8), public ::      cldiceobs(plev)     ! actual W.V. Mixing ratio
real(r8), public ::      numliqobs(plev)     ! actual
real(r8), public ::      numiceobs(plev)     ! actual
real(r8), public ::      precobs(1)          ! observed precipitation
real(r8), public ::      lhflxobs(1)         ! observed surface latent heat flux
real(r8), public ::      heat_glob_scm(1)    ! observed heat total
real(r8), public ::      shflxobs(1)         ! observed surface sensible heat flux
real(r8), public ::      q1obs(plev)         ! observed apparent heat source
real(r8), public ::      q2obs(plev)         ! observed apparent heat sink
real(r8), public ::      tdiff(plev)         ! model minus observed temp
real(r8), public ::      tground(1)          ! ground temperature
real(r8), public ::      psobs               ! observed surface pressure
real(r8), public ::      tobs(plev)          ! observed temperature
real(r8), public ::      tsair(1)            ! air temperature at the surface
real(r8), public ::      udiff(plev)         ! model minus observed uwind
real(r8), public ::      uobs(plev)          ! actual u wind
real(r8), public ::      vdiff(plev)         ! model minus observed vwind
real(r8), public ::      vobs(plev)          ! actual v wind
real(r8), public ::      cldobs(plev)        ! observed cld
real(r8), public ::      clwpobs(plev)       ! observed clwp
real(r8), public ::      aldirobs(1)         ! observed aldir
real(r8), public ::      aldifobs(1)         ! observed aldif
real(r8), public ::      asdirobs(1)         ! observed asdir
real(r8), public ::      asdifobs(1)         ! observed asdif

real(r8), public ::      co2vmrobs(1)        ! observed co2vmr
real(r8), public ::      ch4vmrobs(1)        ! observed ch3vmr
real(r8), public ::      n2ovmrobs(1)        ! observed n2ovmr
real(r8), public ::      f11vmrobs(1)        ! observed f11vmr
real(r8), public ::      f12vmrobs(1)        ! observed f12vmr
real(r8), public ::      soltsiobs(1)        ! observed solar

real(r8), public ::      wfld(plev)          ! Vertical motion (slt)
real(r8), public ::      wfldh(plevp)        ! Vertical motion (slt)
real(r8), public ::      divq(plev,pcnst)    ! Divergence of moisture
real(r8), public ::      divt(plev)          ! Divergence of temperature
real(r8), public ::      divu(plev)          ! Horiz Divergence of E/W
real(r8), public ::      divv(plev)          ! Horiz Divergence of N/S
                                             ! mo_drydep algorithm
real(r8), public, pointer :: loniop(:)
real(r8), public, pointer :: latiop(:)

integer, public ::     iopTimeIdx            ! index into iop dataset
integer, public ::     steplength            ! Length of time-step
integer, public ::     base_date             ! Date in (yyyymmdd) of start time
integer, public ::     base_secs             ! Time of day of start time (sec)

! SCAM public data defaults

logical, public ::  doiopupdate            = .false. ! do we need to read next iop timepoint
logical, public ::  have_lhflx             = .false. ! dataset contains lhflx
logical, public ::  have_shflx             = .false. ! dataset contains shflx
logical, public ::  have_heat_glob         = .false. ! dataset contains heat total
logical, public ::  have_tg                = .false. ! dataset contains tg
logical, public ::  have_tsair             = .false. ! dataset contains tsair
logical, public ::  have_divq              = .false. ! dataset contains divq
logical, public ::  have_divt              = .false. ! dataset contains divt
logical, public ::  have_divq3d            = .false. ! dataset contains divq3d
logical, public ::  have_vertdivu          = .false. ! dataset contains vertdivu
logical, public ::  have_vertdivv          = .false. ! dataset contains vertdivv
logical, public ::  have_vertdivt          = .false. ! dataset contains vertdivt
logical, public ::  have_vertdivq          = .false. ! dataset contains vertdivq
logical, public ::  have_divt3d            = .false. ! dataset contains divt3d
logical, public ::  have_divu3d            = .false. ! dataset contains divu3d
logical, public ::  have_divv3d            = .false. ! dataset contains divv3d
logical, public ::  have_divu              = .false. ! dataset contains divu
logical, public ::  have_divv              = .false. ! dataset contains divv
logical, public ::  have_omega             = .false. ! dataset contains omega
logical, public ::  have_phis              = .false. ! dataset contains phis
logical, public ::  have_ptend             = .false. ! dataset contains ptend
logical, public ::  have_ps                = .false. ! dataset contains ps
logical, public ::  have_q                 = .false. ! dataset contains q
logical, public ::  have_q1                = .false. ! dataset contains Q1
logical, public ::  have_q2                = .false. ! dataset contains Q2
logical, public ::  have_prec              = .false. ! dataset contains prec
logical, public ::  have_t                 = .false. ! dataset contains t
logical, public ::  have_u                 = .false. ! dataset contains u
logical, public ::  have_v                 = .false. ! dataset contains v
logical, public ::  have_cld               = .false. ! dataset contains cld
logical, public ::  have_cldliq            = .false. ! dataset contains cldliq
logical, public ::  have_cldice            = .false. ! dataset contains cldice
logical, public ::  have_numliq            = .false. ! dataset contains numliq
logical, public ::  have_numice            = .false. ! dataset contains numice
logical, public ::  have_clwp              = .false. ! dataset contains clwp
logical, public ::  have_aldir             = .false. ! dataset contains aldir
logical, public ::  have_aldif             = .false. ! dataset contains aldif
logical, public ::  have_asdir             = .false. ! dataset contains asdir
logical, public ::  have_asdif             = .false. ! dataset contains asdif
logical, public ::  use_camiop             = .false. ! use cam generated forcing
logical, public ::  use_3dfrc              = .false. ! use 3d forcing
logical, public ::  isrestart              = .false. ! If this is a restart step or not

! SCAM namelist defaults

logical, public ::  scm_backfill_iop_w_init = .false. ! Backfill missing IOP data from initial file
logical, public ::  scm_relaxation         = .false. ! Use relaxation
logical, public ::  scm_crm_mode           = .false. ! Use column radiation mode
logical, public ::  scm_cambfb_mode        = .false. ! Use extra CAM IOP fields to assure bit for bit match with CAM run
logical, public ::  scm_use_obs_T          = .false. ! Use the SCAM-IOP observed T at each timestep instead of forecasting.
logical, public ::  scm_force_latlon       = .false. ! force scam to use the lat lon fields specified in the namelist not closest
real(r8), public              ::  scm_relaxation_low      ! lowest level to apply relaxation
real(r8), public              ::  scm_relaxation_high     ! highest level to apply relaxation
real(r8), public              ::  scm_relax_top_p         = 0._r8       ! upper bound for scm relaxation
real(r8), public              ::  scm_relax_bot_p         = huge(1._r8) ! lower bound for scm relaxation
real(r8), public              ::  scm_relax_tau_sec       = 10800._r8   ! relaxation time constant (sec)

! +++BPM:
! modification... allow a linear ramp in relaxation time scale:
logical, public :: scm_relax_linear = .false.
real(r8), public    :: scm_relax_tau_bot_sec = 10800._r8
real(r8), public    :: scm_relax_tau_top_sec = 10800._r8
character(len=26), public  :: scm_relax_fincl(pcnst)

!
! note that scm_use_obs_uv is set to true to be consistent with CAM BFB testing
!

logical, public ::  scm_use_obs_uv         = .true. ! Use the SCAM-IOP observed u,v at each time step instead of forecasting.

logical, public ::  scm_use_obs_qv         = .false. ! Use the SCAM-IOP observed qv  at each time step instead of forecasting.
logical, public ::  scm_use_3dfrc          = .false. ! Use CAMIOP 3d forcing if true, else use dycore vertical plus horizontal
logical, public ::  scm_iop_lhflxshflxTg   = .false. !turn off LW rad
logical, public ::  scm_iop_Tg             = .false. !turn off LW rad

character(len=200), public ::  scm_clubb_iop_name   ! IOP name for CLUBB

integer, allocatable, public :: tsec(:)
integer, public :: ntime

!=======================================================================
contains
!=======================================================================

subroutine scam_readnl(nlfile,single_column_in,scmlat_in,scmlon_in)

  use namelist_utils,  only: find_group_name
  use units,           only: getunit, freeunit
  use dycore,          only: dycore_is
  use wrap_nf,         only: wrap_open


!---------------------------Arguments-----------------------------------

  character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input  (nlfile=atm_in)
  logical,          intent(in) :: single_column_in
  real(r8),         intent(in) :: scmlat_in
  real(r8),         intent(in) :: scmlon_in

  ! Local variables
  character(len=*), parameter :: sub = 'scam_readnl'
  integer :: unitn, ierr, i
  integer  :: ncid
  integer  :: iatt
  logical  :: adv

! this list should include any variable that you might want to include in the namelist
  namelist /scam_nl/ iopfile, scm_iop_lhflxshflxTg, scm_iop_Tg, scm_relaxation, &
       scm_relax_top_p,scm_relax_bot_p,scm_relax_tau_sec, &
       scm_cambfb_mode,scm_crm_mode,scm_zadv_uv,scm_zadv_T,scm_zadv_q,&
       scm_use_obs_T, scm_use_obs_uv, scm_use_obs_qv, scm_use_3dfrc, &
       scm_relax_linear, scm_relax_tau_top_sec, &
       scm_relax_tau_bot_sec, scm_force_latlon, scm_relax_fincl, &
       scm_backfill_iop_w_init

  single_column=single_column_in

  iopfile            = ' '
  scm_clubb_iop_name = ' '
  scm_relax_fincl(:) = ' '
  if( single_column ) then
     if( npes>1) call endrun('SCAM_READNL: SCAM doesnt support using more than 1 pe.')

     if ( .not. (dycore_is('EUL') .or. dycore_is('SE')) .or. plon /= 1 .or. plat /=1 ) then
        call endrun('SCAM_SETOPTS: must compile model for SCAM mode when namelist parameter single_column is .true.')
     endif

     scmlat=scmlat_in
     scmlon=scmlon_in

     if( scmlat < -90._r8 .or. scmlat > 90._r8 ) then
        call endrun('SCAM_READNL: SCMLAT must be between -90. and 90. degrees.')
     elseif( scmlon < 0._r8 .or. scmlon > 360._r8 ) then
        call endrun('SCAM_READNL: SCMLON must be between 0. and 360. degrees.')
     end if

     ! Read namelist
     if (masterproc) then
        unitn = getunit()
        open( unitn, file=trim(nlfile), status='old' )
        call find_group_name(unitn, 'scam_nl', status=ierr)
        if (ierr == 0) then
           read(unitn, scam_nl, iostat=ierr)
           if (ierr /= 0) then
              call endrun(sub // ':: ERROR reading namelist')
           end if
        end if
        close(unitn)
        call freeunit(unitn)
     end if

     ! Error checking:

     iopfile = trim(iopfile)
     if( iopfile /= "" ) then
        use_iop = .true.
     else
        call endrun('SCAM_READNL: must specify IOP file for single column mode')
     endif

     call wrap_open( iopfile, NF90_NOWRITE, ncid )

     if( nf90_inquire_attribute( ncid, NF90_GLOBAL, 'CAM_GENERATED_FORCING', iatt ) ==  NF90_NOERR ) then
        use_camiop = .true.
     else
        use_camiop = .false.
     endif

     ! If we are not forcing the lat and lon from the namelist use the closest lat and lon that is found in the IOP file.
     if (.not.scm_force_latlon) then
        call shr_scam_GetCloseLatLon( ncid, scmlat, scmlon, closeioplat, closeioplon, closeioplatidx, closeioplonidx )
        write(iulog,*) 'SCAM_READNL: using closest IOP column to lat/lon specified in drv_in'
        write(iulog,*) '   requested lat,lon    =',scmlat,', ',scmlon
        write(iulog,*) '   closest IOP lat,lon  =',closeioplat,', ',closeioplon
        scmlat = closeioplat
        scmlon = closeioplon
     end if

     if (masterproc) then
        write (iulog,*) 'Single Column Model Options: '
        write (iulog,*) '============================='
        write (iulog,*) '  iopfile                     = ',trim(iopfile)
        write (iulog,*) '  scm_backfill_iop_w_init     = ',scm_backfill_iop_w_init
        write (iulog,*) '  scm_cambfb_mode             = ',scm_cambfb_mode
        write (iulog,*) '  scm_crm_mode                = ',scm_crm_mode
        write (iulog,*) '  scm_force_latlon            = ',scm_force_latlon
        write (iulog,*) '  scm_iop_Tg                  = ',scm_iop_Tg
        write (iulog,*) '  scm_iop_lhflxshflxTg        = ',scm_iop_lhflxshflxTg
        write (iulog,*) '  scm_relaxation              = ',scm_relaxation
        write (iulog,*) '  scm_relax_bot_p             = ',scm_relax_bot_p
        write (iulog,*) '  scm_relax_linear            = ',scm_relax_linear
        write (iulog,*) '  scm_relax_tau_bot_sec       = ',scm_relax_tau_bot_sec
        write (iulog,*) '  scm_relax_tau_sec           = ',scm_relax_tau_sec
        write (iulog,*) '  scm_relax_tau_top_sec       = ',scm_relax_tau_top_sec
        write (iulog,*) '  scm_relax_top_p             = ',scm_relax_top_p
        write (iulog,*) '  scm_use_obs_T               = ',scm_use_obs_T
        write (iulog,*) '  scm_use_3dfrc               = ',scm_use_3dfrc
        write (iulog,*) '  scm_use_obs_qv              = ',scm_use_obs_qv
        write (iulog,*) '  scm_use_obs_uv              = ',scm_use_obs_uv
        write (iulog,*) '  scm_zadv_T                  = ',trim(scm_zadv_T)
        write (iulog,*) '  scm_zadv_q                  = ',trim(scm_zadv_q)
        write (iulog,*) '  scm_zadv_uv                 = ',trim(scm_zadv_uv)
        write (iulog,*) '  scm_relax_finc: '
        ! output scm_relax_fincl character array
        do i=1,pcnst
           if (scm_relax_fincl(i) /= '') then
              adv = mod(i,4)==0
              if (adv) then
                 write (iulog, "(A18)") "'"//trim(scm_relax_fincl(i))//"',"
              else
                 write (iulog, "(A18)", ADVANCE="NO") "'"//trim(scm_relax_fincl(i))//"',"
              end if
           else
              exit
           end if
        end do
        print *
     end if
  end if

end subroutine scam_readnl
subroutine readiopdata(hyam, hybm, hyai, hybi, ps0)
!-----------------------------------------------------------------------
!
!     Open and read netCDF file containing initial IOP  conditions
!
!---------------------------Code history--------------------------------
!
!     Written by J.  Truesdale    August, 1996, revised January, 1998
!
!-----------------------------------------------------------------------
        use getinterpnetcdfdata, only: getinterpncdata
        use string_utils,        only: to_lower
        use wrap_nf,             only: wrap_inq_dimid,wrap_get_vara_realx
!-----------------------------------------------------------------------
   implicit none

   character(len=*), parameter ::  sub = "read_iop_data"
!
!------------------------------Input Arguments--------------------------
!
   real(r8),intent(in) :: hyam(plev),hybm(plev),hyai(plevp),hybi(plevp),ps0
!
!------------------------------Locals-----------------------------------
!
   integer :: NCID, status
   integer :: time_dimID, lev_dimID,  lev_varID, varid
   integer :: i,j
   integer :: nlev
   integer :: total_levs
   integer :: u_attlen

   integer :: k, m
   integer :: icldliq,icldice
   integer :: inumliq,inumice

   logical :: have_srf              ! value at surface is available
   logical :: fill_ends             !
   logical :: have_cnst(pcnst)
   real(r8) :: dummy
   real(r8) :: srf(1)                  ! value at surface
   real(r8) :: hyamiop(plev)  ! a hybrid coef midpoint
   real(r8) :: hybmiop(plev)  ! b hybrid coef midpoint
   real(r8) :: pmid(plev)  ! pressure at model levels (time n)
   real(r8) :: pint(plevp) ! pressure at model interfaces (n  )
   real(r8) :: pdel(plev)  ! pdel(k)   = pint  (k+1)-pint  (k)
   real(r8) :: weight
   real(r8) :: tmpdata(1)
   real(r8) :: coldata(plev)
   real(r8), allocatable :: dplevs( : )
   integer :: strt4(4),cnt4(4)
   integer :: nstep
   integer :: ios
   character(len=128) :: units ! Units

   nstep = get_nstep()
   fill_ends= .false.

!
!     Open IOP dataset
!
   call handle_ncerr( nf90_open (iopfile, 0, ncid),&
       'ERROR - scamMod.F90:readiopdata', __LINE__)

!
!     if the dataset is a CAM generated dataset set use_camiop to true
!       CAM IOP datasets have a global attribute called CAM_GENERATED_IOP
!
   if ( nf90_inquire_attribute( ncid, NF90_GLOBAL, 'CAM_GENERATED_FORCING', attnum=i )== NF90_NOERR ) then
      use_camiop = .true.
   else
      use_camiop = .false.
   endif

!=====================================================================
!
!     Read time variables


   status = nf90_inq_dimid (ncid, 'time', time_dimID )
   if (status /= NF90_NOERR) then
      status = nf90_inq_dimid (ncid, 'tsec', time_dimID )
      if (status /= NF90_NOERR) then
         if (masterproc) write(iulog,*) sub//':ERROR - Could not find dimension ID for time/tsec'
         status = NF90_CLOSE ( ncid )
         call endrun(sub // ':ERROR - time/tsec must be present on the IOP file.')
      end if
   end if

   call handle_ncerr( nf90_inquire_dimension( ncid, time_dimID, len=ntime ),&
         'Error - scamMod.F90:readiopdata unable to find time dimension', __LINE__)

!
!======================================================
!     read level data
!
   status = NF90_INQ_DIMID( ncid, 'lev', lev_dimID )
   if ( status /= nf90_noerr ) then
      if (masterproc) write(iulog,*) sub//':ERROR - Could not find variable dim ID  for lev'
      status = NF90_CLOSE ( ncid )
      call endrun(sub // ':ERROR - Could not find variable dim ID for lev')
   end if

   call handle_ncerr( nf90_inquire_dimension( ncid, lev_dimID, len=nlev ),&
         'Error - scamMod.f90:readiopdata unable to find level dimension', __LINE__)

   allocate(dplevs(nlev+1),stat=ios)
   if( ios /= 0 ) then
      write(iulog,*) sub//':ERROR: failed to allocate dplevs; error = ',ios
      call endrun(sub//':ERROR:readiopdata failed to allocate dplevs')
   end if

   status = NF90_INQ_VARID( ncid, 'lev', lev_varID )
   if ( status /= nf90_noerr ) then
      if (masterproc) write(iulog,*) sub//':ERROR - scamMod.F90:readiopdata:Could not find variable ID for lev'
      status = NF90_CLOSE ( ncid )
      call endrun(sub//':ERROR:ould not find variable ID for lev')
   end if

   call handle_ncerr( nf90_get_var (ncid, lev_varID, dplevs(:nlev)),&
                    'Error - scamMod.F90:readiopdata unable to read pressure levels', __LINE__)
!
!CAM generated forcing already has pressure on millibars convert standard IOP if needed.
!
   call handle_ncerr(nf90_inquire_attribute(ncid, lev_varID, 'units', len=u_attlen),&
                    'Error - scamMod.F90:readiopdata unable to find units attribute', __LINE__)
   call handle_ncerr(nf90_get_att(ncid, lev_varID, 'units', units),&
                    'Error - scamMod.F90:readiopdata unable to read units attribute', __LINE__)
   units=trim(to_lower(units(1:u_attlen)))

   if ( units=='pa' .or. units=='pascal' .or. units=='pascals' ) then
!
!     convert pressure from Pascals to Millibars ( lev is expressed in pascals in iop datasets )
!
      do i=1,nlev
         dplevs( i ) = dplevs( i )/100._r8
      end do
   endif

   status = nf90_inq_varid( ncid, 'Ps', varid   )
   if ( status /= nf90_noerr ) then
      have_ps= .false.
      if (masterproc) write(iulog,*) sub//':Could not find variable Ps'
      if ( .not. scm_backfill_iop_w_init ) then
         status = NF90_CLOSE( ncid )
         call endrun(sub//':ERROR :IOP file must contain Surface Pressure (Ps) variable')
      else
         if ( is_first_step() .and. masterproc) write(iulog,*) 'Using surface pressure value from IC file if present'
      endif
   else
      call get_start_count(ncid, varid, scmlat, scmlon, ioptimeidx, strt4, cnt4)
      status = nf90_get_var(ncid, varid, psobs, strt4)
      have_ps = .true.
   endif


!  If the IOP dataset has hyam,hybm,etc it is assumed to be a hybrid level
!  dataset

   status =  nf90_inq_varid( ncid, 'hyam', varid   )
   if ( status == nf90_noerr .and. have_ps) then
      call get_start_count(ncid, varid, scmlat, scmlon, ioptimeidx, strt4, cnt4)
      status = nf90_get_var(ncid, varid, hyamiop, strt4)
      status =  nf90_inq_varid( ncid, 'hybm', varid   )
      status = nf90_get_var(ncid, varid, hybmiop, strt4)
      do i = 1, nlev
         dplevs( i ) = 1000.0_r8 * hyamiop( i ) + psobs * hybmiop( i ) / 100.0_r8
      end do
   endif

!     add the surface pressure to the pressure level data, so that
!     surface boundary condition will be set properly,
!     making sure that it is the highest pressure in the array.
!

   total_levs = nlev+1
   dplevs(nlev+1) = psobs/100.0_r8 ! ps is expressed in pascals
   do i= nlev, 1, -1
      if ( dplevs(i) > psobs/100.0_r8) then
         total_levs = i
         dplevs(i) = psobs/100.0_r8
      end if
   end do
   if (.not. use_camiop ) then
      nlev = total_levs
   endif
   if ( nlev == 1 ) then
      if (masterproc) write(iulog,*) sub//':Error - scamMod.F90:readiopdata: Ps too low!'
      call endrun(sub//':ERROR:Ps value on datasets is incongurent with levs data - mismatch in units?')
   endif

!=====================================================================
!get global vmrs from camiop file
   status =  nf90_inq_varid( ncid, 'co2vmr', varid   )
   if ( status == nf90_noerr) then
      call get_start_count(ncid, varid, scmlat, scmlon, ioptimeidx, strt4, cnt4)
      call wrap_get_vara_realx (ncid,varid,strt4,cnt4,co2vmrobs)
   else
      if (is_first_step()) write(iulog,*)'using column value of co2vmr from boundary data as global volume mixing ratio'
   end if
   status =  nf90_inq_varid( ncid, 'ch4vmr', varid   )
   if ( status == nf90_noerr) then
      call get_start_count(ncid, varid, scmlat, scmlon, ioptimeidx, strt4, cnt4)
      call wrap_get_vara_realx (ncid,varid,strt4,cnt4,ch4vmrobs)
   else
      if (is_first_step()) write(iulog,*)'using column value of ch4vmr from boundary data as global volume mixing ratio'
   end if
   status =  nf90_inq_varid( ncid, 'n2ovmr', varid   )
   if ( status == nf90_noerr) then
      call get_start_count(ncid, varid, scmlat, scmlon, ioptimeidx, strt4, cnt4)
      call wrap_get_vara_realx (ncid,varid,strt4,cnt4,n2ovmrobs)
   else
      if (is_first_step()) write(iulog,*)'using column value of n2ovmr from boundary data as global volume mixing ratio'
   end if
   status =  nf90_inq_varid( ncid, 'f11vmr', varid   )
   if ( status == nf90_noerr) then
      call get_start_count(ncid, varid, scmlat, scmlon, ioptimeidx, strt4, cnt4)
      call wrap_get_vara_realx (ncid,varid,strt4,cnt4,f11vmrobs)
   else
      if (is_first_step()) write(iulog,*)'using column value of f11vmr from boundary data as global volume mixing ratio'
   end if
   status =  nf90_inq_varid( ncid, 'f12vmr', varid   )
   if ( status == nf90_noerr) then
      call get_start_count(ncid, varid, scmlat, scmlon, ioptimeidx, strt4, cnt4)
      call wrap_get_vara_realx (ncid,varid,strt4,cnt4,f12vmrobs)
   else
      if (is_first_step()) write(iulog,*)'using column value of f12vmr from boundary data as global volume mixing ratio'
   end if
   status =  nf90_inq_varid( ncid, 'soltsi', varid   )
   if ( status == nf90_noerr) then
      call get_start_count(ncid, varid, scmlat, scmlon, ioptimeidx, strt4, cnt4)
      call wrap_get_vara_realx (ncid,varid,strt4,cnt4,soltsiobs)
   else
      if (is_first_step()) write(iulog,*)'using column value of soltsi from boundary data as global solar tsi'
   end if
!=====================================================================
!get state variables from camiop file

   status =  nf90_inq_varid( ncid, 'Tsair', varid   )
   if ( status /= nf90_noerr ) then
      have_tsair = .false.
   else
      call get_start_count(ncid, varid, scmlat, scmlon, ioptimeidx, strt4, cnt4)
      call wrap_get_vara_realx (ncid,varid,strt4,cnt4,tsair)
      have_tsair = .true.
   endif
!
!      read in Tobs  For cam generated iop readin small t to avoid confusion
!      with capital T defined in cam
!
   tobs(:)= 0._r8

   if ( use_camiop ) then
     call getinterpncdata( ncid, scmlat, scmlon, ioptimeidx,'t', have_tsair, &
          tsair(1), fill_ends, scm_crm_mode, &
          dplevs, nlev,psobs, hyam, hybm,tobs, status )
   else
     call getinterpncdata( ncid, scmlat, scmlon, ioptimeidx,'T', have_tsair, &
          tsair(1), fill_ends, scm_crm_mode, &
          dplevs, nlev,psobs, hyam, hybm, tobs, status )
   endif
   if ( status /= nf90_noerr ) then
      have_t = .false.
      if (masterproc) write(iulog,*) sub//':Could not find variable T on IOP file'
      if ( scm_backfill_iop_w_init ) then
         if (masterproc) write(iulog,*) sub//':Using value of T(tobs) from IC file if it exists'
      else
         if (masterproc) write(iulog,*) sub//':set tobs to 0.'
      endif
!
!     set T3 to Tobs on first time step
!
   else
      have_t = .true.
   endif

   status = nf90_inq_varid( ncid, 'Tg', varid   )
   if (status /= nf90_noerr) then
      if (masterproc) write(iulog,*) sub//':Could not find variable Tg on IOP dataset'
      if ( have_tsair ) then
         if (masterproc) write(iulog,*) sub//':Using Tsair'
         tground = tsair     ! use surface value from T field
         have_Tg = .true.
      else
         have_Tg = .true.
         if (masterproc) write(iulog,*) sub//':Using T at lowest level from IOP dataset'
         tground = tobs(plev)
      endif
   else
      call get_start_count(ncid, varid, scmlat, scmlon, ioptimeidx, strt4, cnt4)
      call wrap_get_vara_realx (ncid,varid,strt4,cnt4,tground)
      have_Tg = .true.
   endif

   status = nf90_inq_varid( ncid, 'qsrf', varid   )

   if ( status /= nf90_noerr ) then
      have_srf = .false.
   else
      call get_start_count(ncid, varid, scmlat, scmlon, ioptimeidx, strt4, cnt4)
      status = nf90_get_var(ncid, varid, srf(1), strt4)
      have_srf = .true.
   endif

   qobs(:)= 0._r8
   call getinterpncdata( ncid, scmlat, scmlon, ioptimeidx,  'q', have_srf, &
      srf(1), fill_ends, scm_crm_mode, &
      dplevs, nlev,psobs, hyam, hybm, qobs, status )
   if ( status /= nf90_noerr ) then
      have_q = .false.
      if (masterproc) write(iulog,*) sub//':Could not find variable q on IOP file'
      if ( scm_backfill_iop_w_init ) then
         if (masterproc) write(iulog,*) sub//':Using values for q from IC file if available'
      else
         if (masterproc) write(iulog,*) sub//':Setting qobs to 0.'
      endif
   else
      have_q = .true.
   endif

   cldobs = 0._r8
   call getinterpncdata( ncid, scmlat, scmlon, ioptimeidx,  'cld', .false., &
      dummy, fill_ends, scm_crm_mode, dplevs, nlev,psobs, hyam, hybm, cldobs, status )
   if ( status /= nf90_noerr ) then
      have_cld = .false.
   else
      have_cld = .true.
   endif

   clwpobs = 0._r8
   call getinterpncdata( ncid, scmlat, scmlon, ioptimeidx,  'clwp', .false., &
      dummy, fill_ends, scm_crm_mode, dplevs, nlev,psobs, hyam, hybm, clwpobs, status )
   if ( status /= nf90_noerr ) then
      have_clwp = .false.
   else
      have_clwp = .true.
   endif

!
!	read divq (horizontal advection)
!
   status = nf90_inq_varid( ncid, 'divqsrf', varid   )
   if ( status /= nf90_noerr ) then
      have_srf = .false.
   else
      call get_start_count(ncid, varid, scmlat, scmlon, ioptimeidx, strt4, cnt4)
      status = nf90_get_var(ncid, varid, srf(1), strt4)
      have_srf = .true.
   endif

   divq(:,:)=0._r8

   call getinterpncdata( ncid, scmlat, scmlon, ioptimeidx, &
        'divq', have_srf, srf(1), fill_ends, scm_crm_mode, &
        dplevs, nlev,psobs, hyam, hybm, divq(:,1), status )
   if ( status /= nf90_noerr ) then
      have_divq = .false.
   else
      have_divq = .true.
   endif

!
!     read vertdivq if available
!
   status = nf90_inq_varid( ncid, 'vertdivqsrf', varid   )
   if ( status /= nf90_noerr ) then
      have_srf = .false.
   else
      call get_start_count(ncid, varid, scmlat, scmlon, ioptimeidx, strt4, cnt4)
      status = nf90_get_var(ncid, varid, srf(1), strt4)
      have_srf = .true.
   endif

   vertdivq=0._r8

   call getinterpncdata( ncid, scmlat, scmlon, ioptimeidx, 'vertdivq', &
        have_srf, srf(1), fill_ends, scm_crm_mode, &
        dplevs, nlev,psobs, hyam, hybm, vertdivq(:,1), status )
   if ( status /= nf90_noerr ) then
      have_vertdivq = .false.
   else
      have_vertdivq = .true.
   endif

   status = nf90_inq_varid( ncid, 'vertdivqsrf', varid   )
   if ( status /= nf90_noerr ) then
      have_srf = .false.
   else
      call get_start_count(ncid, varid, scmlat, scmlon, ioptimeidx, strt4, cnt4)
      status = nf90_get_var(ncid, varid, srf(1), strt4)
      have_srf = .true.
   endif
!
!   add calls to get dynamics tendencies for all prognostic consts
!
   divq3d=0._r8

   do m = 1, pcnst
      call getinterpncdata( ncid, scmlat, scmlon, ioptimeidx, trim(cnst_name(m))//'_dten', &
      have_srf, srf(1), fill_ends, scm_crm_mode, &
      dplevs, nlev,psobs, hyam, hybm, divq3d(:,m), status )
      write(iulog,*)'checking ',trim(cnst_name(m))//'_dten',status
      if ( status /= nf90_noerr ) then
         have_cnst(m) = .false.
         divq3d(1:,m)=0._r8
      else
         if (m==1) have_divq3d = .true.
         have_cnst(m) = .true.
      endif

      coldata = 0._r8
      call getinterpncdata( ncid, scmlat, scmlon, ioptimeidx, trim(cnst_name(m))//'_dqfx', &
      have_srf, srf(1), fill_ends, scm_crm_mode, &
      dplevs, nlev,psobs, hyam, hybm, coldata, status )
      if ( STATUS /= NF90_NOERR ) then
         dqfxcam(1,:,m)=0._r8
      else
         dqfxcam(1,:,m)=coldata(:)
      endif

      tmpdata = 0._r8
      call getinterpncdata( ncid, scmlat, scmlon, ioptimeidx, trim(cnst_name(m))//'_alph', &
      have_srf, srf(1), fill_ends, scm_crm_mode, &
      dplevs, nlev,psobs, hyam, hybm, tmpdata, status )
      if ( status /= nf90_noerr ) then
         alphacam(m)=0._r8
      else
          alphacam(m)=tmpdata(1)
      endif

   end do


   numliqobs = 0._r8
   call cnst_get_ind('NUMLIQ', inumliq, abort=.false.)
   if ( inumliq > 0 ) then
      have_srf = .false.
      call getinterpncdata( ncid, scmlat, scmlon, ioptimeidx, 'NUMLIQ', &
           have_srf, srf(1), fill_ends, scm_crm_mode, &
           dplevs, nlev,psobs, hyam, hybm, numliqobs, status )
      if ( status /= nf90_noerr ) then
         have_numliq = .false.
      else
         have_numliq = .true.
      endif
   else
         have_numliq = .false.
   end if

   have_srf = .false.

   cldliqobs = 0._r8
   call cnst_get_ind('CLDLIQ', icldliq, abort=.false.)
   if ( icldliq > 0 ) then
      call getinterpncdata( ncid, scmlat, scmlon, ioptimeidx, 'CLDLIQ', &
           have_srf, srf(1), fill_ends, scm_crm_mode, &
           dplevs, nlev,psobs, hyam, hybm, cldliqobs, status )
      if ( status /= nf90_noerr ) then
         have_cldliq = .false.
      else
         have_cldliq = .true.
      endif
   else
         have_cldliq = .false.
   endif

   cldiceobs = 0._r8
   call cnst_get_ind('CLDICE', icldice, abort=.false.)
   if ( icldice > 0 ) then
      call getinterpncdata( ncid, scmlat, scmlon, ioptimeidx, 'CLDICE', &
           have_srf, srf(1), fill_ends, scm_crm_mode, &
           dplevs, nlev,psobs, hyam, hybm, cldiceobs, status )
      if ( status /= nf90_noerr ) then
         have_cldice = .false.
      else
         have_cldice = .true.
      endif
   else
      have_cldice = .false.
   endif

   numiceobs = 0._r8
   call cnst_get_ind('NUMICE', inumice, abort=.false.)
   if ( inumice > 0 ) then
      have_srf = .false.
      call getinterpncdata( ncid, scmlat, scmlon, ioptimeidx, 'NUMICE', &
         have_srf, srf(1), fill_ends, scm_crm_mode, &
         dplevs, nlev,psobs, hyam, hybm, numiceobs, status )
      if ( status /= nf90_noerr ) then
         have_numice = .false.
      else
         have_numice = .true.
      endif
   else
      have_numice = .false.
   end if

!
!	read divu (optional field)
!
   status = nf90_inq_varid( ncid, 'divusrf', varid   )
   if ( status /= nf90_noerr ) then
      have_srf = .false.
   else
      call get_start_count(ncid, varid, scmlat, scmlon, ioptimeidx, strt4, cnt4)
      status = nf90_get_var(ncid, varid, srf(1), strt4)
      have_srf = .true.
   endif

   divu = 0._r8
   call getinterpncdata( ncid, scmlat, scmlon, ioptimeidx, 'divu', &
      have_srf, srf(1), fill_ends, scm_crm_mode, &
      dplevs, nlev,psobs, hyam, hybm, divu, status )
   if ( status /= nf90_noerr ) then
      have_divu = .false.
   else
      have_divu = .true.
   endif
!
!	read divv (optional field)
!
   status = nf90_inq_varid( ncid, 'divvsrf', varid   )
   if ( status /= nf90_noerr ) then
      have_srf = .false.
   else
      call get_start_count(ncid, varid, scmlat, scmlon, ioptimeidx, strt4, cnt4)
      status = nf90_get_var(ncid, varid, srf(1), strt4)
      have_srf = .true.
   endif

   divv = 0._r8
   call getinterpncdata( ncid, scmlat, scmlon, ioptimeidx, 'divv', &
      have_srf, srf(1), fill_ends, scm_crm_mode, &
      dplevs, nlev,psobs, hyam, hybm, divv, status )
   if ( status /= nf90_noerr ) then
      have_divv = .false.
   else
      have_divv = .true.
   endif
!
!	read divt (optional field)
!
   status = nf90_inq_varid( ncid, 'divtsrf', varid   )
   if ( status /= nf90_noerr ) then
      have_srf = .false.
   else
      call get_start_count(ncid, varid, scmlat, scmlon, ioptimeidx, strt4, cnt4)
      status = nf90_get_var(ncid, varid, srf(1), strt4)
      have_srf = .true.
   endif

   divt=0._r8
   call getinterpncdata( ncid, scmlat, scmlon, ioptimeidx, &
      'divT', have_srf, srf(1), fill_ends, scm_crm_mode, &
      dplevs, nlev,psobs, hyam, hybm, divt, status )
   if ( status /= nf90_noerr ) then
      have_divt = .false.
   else
      have_divt = .true.
   endif

!
!     read vertdivt if available
!
   status = nf90_inq_varid( ncid, 'vertdivTsrf', varid   )
   if ( status /= nf90_noerr ) then
      have_srf = .false.
   else
      call get_start_count(ncid, varid, scmlat, scmlon, ioptimeidx, strt4, cnt4)
      status = nf90_get_var(ncid, varid, srf(1), strt4)
      have_srf = .true.
   endif

   vertdivt=0._r8
   call getinterpncdata( ncid, scmlat, scmlon, ioptimeidx, 'vertdivTx', &
      have_srf, srf(1), fill_ends, scm_crm_mode, &
      dplevs, nlev,psobs, hyam, hybm, vertdivt, status )
   if ( status /= nf90_noerr ) then
      call getinterpncdata( ncid, scmlat, scmlon, ioptimeidx, 'vertdivT', &
           have_srf, srf(1), fill_ends, scm_crm_mode, &
           dplevs, nlev,psobs, hyam, hybm, vertdivt, status )
      if ( status /= nf90_noerr ) then
         have_vertdivt = .false.
      else
         have_vertdivt = .true.
      endif
   else
      have_vertdivt = .true.
   endif
!
!	read divt3d (combined vertical/horizontal advection)
!      (optional field)

   status = nf90_inq_varid( ncid, 'divT3dsrf', varid   )
   if ( status /= nf90_noerr ) then
      have_srf = .false.
   else
      call get_start_count(ncid, varid, scmlat, scmlon, ioptimeidx, strt4, cnt4)
      status = nf90_get_var(ncid, varid, srf(1), strt4)
      have_srf = .true.
   endif

   divT3d = 0._r8

   call getinterpncdata( ncid, scmlat, scmlon, ioptimeidx, 'divT3d', &
      have_srf, srf(1), fill_ends, scm_crm_mode, &
      dplevs, nlev,psobs, hyam, hybm, divt3d, status )
      write(iulog,*)'checking divT3d:',status,nf90_noerr
   if ( status /= nf90_noerr ) then
      have_divt3d = .false.
   else
      have_divt3d = .true.
   endif

   divU3d = 0._r8

   call getinterpncdata( ncid, scmlat, scmlon, ioptimeidx, 'divU3d', &
      have_srf, srf(1), fill_ends, scm_crm_mode, &
      dplevs, nlev,psobs, hyam, hybm, divu3d, status )
   if ( status /= nf90_noerr ) then
      have_divu3d = .false.
   else
      have_divu3d = .true.
   endif

   divV3d = 0._r8

   call getinterpncdata( ncid, scmlat, scmlon, ioptimeidx, 'divV3d', &
      have_srf, srf(1), fill_ends, scm_crm_mode, &
      dplevs, nlev,psobs, hyam, hybm, divv3d, status )
   if ( status /= nf90_noerr ) then
      have_divv3d = .false.
   else
      have_divv3d = .true.
   endif

   status = nf90_inq_varid( ncid, 'Ptend', varid   )
   if ( status /= nf90_noerr ) then
      have_ptend = .false.
      if (masterproc) write(iulog,*) sub//':Could not find variable Ptend. Setting to zero'
      ptend = 0.0_r8
   else
      call get_start_count(ncid, varid, scmlat, scmlon, ioptimeidx, strt4, cnt4)
      status = nf90_get_var(ncid, varid, srf(1), strt4)
      have_ptend = .true.
      ptend= srf(1)
   endif

   wfld=0._r8

   call getinterpncdata( ncid, scmlat, scmlon, ioptimeidx, &
      'omega', .true., ptend, fill_ends, scm_crm_mode, &
      dplevs, nlev,psobs, hyam, hybm, wfld, status )
   if ( status /= nf90_noerr ) then
      have_omega = .false.
      if (masterproc) write(iulog,*) sub//':Could not find variable omega on IOP'
      if ( scm_backfill_iop_w_init ) then
         if (masterproc) write(iulog,*) sub//'Using omega from IC file'
      else
         if (masterproc) write(iulog,*) sub//'setting Omega to 0. throughout the column'
      endif
   else
      have_omega = .true.
   endif
   call plevs0(plev, psobs, ps0, hyam, hybm, hyai, hybi, pint, pmid ,pdel)
!
! Build interface vector for the specified omega profile
! (weighted average in pressure of specified level values)
!
   wfldh(:) = 0.0_r8

   do k=2,plev
      weight = (pint(k) - pmid(k-1))/(pmid(k) - pmid(k-1))
      wfldh(k) = (1.0_r8 - weight)*wfld(k-1) + weight*wfld(k)
   end do

   status = nf90_inq_varid( ncid, 'usrf', varid   )
   if ( status /= nf90_noerr ) then
      have_srf = .false.
   else
      call get_start_count(ncid, varid, scmlat, scmlon, ioptimeidx, strt4, cnt4)
      call wrap_get_vara_realx (ncid,varid,strt4,cnt4,srf)
      have_srf = .true.
   endif

   uobs=0._r8

   call getinterpncdata( ncid, scmlat, scmlon, ioptimeidx, &
      'u', have_srf, srf(1), fill_ends, scm_crm_mode, &
      dplevs, nlev,psobs, hyam, hybm, uobs, status )
   if ( status /= nf90_noerr ) then
      have_u = .false.
   else
      have_u = .true.
   endif

   status = nf90_inq_varid( ncid, 'vsrf', varid   )
   if ( status /= nf90_noerr ) then
      have_srf = .false.
   else
      call get_start_count(ncid, varid, scmlat, scmlon, ioptimeidx, strt4, cnt4)
      call wrap_get_vara_realx (ncid,varid,strt4,cnt4,srf)
      have_srf = .true.
   endif

   vobs=0._r8

   call getinterpncdata( ncid, scmlat, scmlon, ioptimeidx, &
      'v', have_srf, srf(1), fill_ends, scm_crm_mode, &
      dplevs, nlev,psobs, hyam, hybm, vobs, status )
   if ( status /= nf90_noerr ) then
      have_v = .false.
   else
      have_v = .true.
   endif

   status = nf90_inq_varid( ncid, 'Prec', varid   )
   if ( status /= nf90_noerr ) then
      have_prec = .false.
   else
      call get_start_count(ncid, varid, scmlat, scmlon, ioptimeidx, strt4, cnt4)
      call wrap_get_vara_realx (ncid,varid,strt4,cnt4,precobs)
      have_prec = .true.
   endif

   q1obs = 0._r8

   call getinterpncdata( ncid, scmlat, scmlon, ioptimeidx, 'Q1', &
      .false., dummy, fill_ends, scm_crm_mode, & ! datasets don't contain Q1 at surface
      dplevs, nlev,psobs, hyam, hybm, q1obs, status )
   if ( status /= nf90_noerr ) then
      have_q1 = .false.
   else
      have_q1 = .true.
   endif

   q1obs = 0._r8

   call getinterpncdata( ncid, scmlat, scmlon, ioptimeidx, 'Q2', &
      .false., dummy, fill_ends, scm_crm_mode, & ! datasets don't contain Q2 at surface
      dplevs, nlev,psobs, hyam, hybm, q1obs, status )
   if ( status /= nf90_noerr ) then
      have_q2 = .false.
   else
      have_q2 = .true.
   endif

!  Test for BOTH 'lhflx' and 'lh' without overwriting 'have_lhflx'.
!  Analagous changes made for the surface heat flux

   status = nf90_inq_varid( ncid, 'lhflx', varid   )
   if ( status /= nf90_noerr ) then
      status = nf90_inq_varid( ncid, 'lh', varid   )
      if ( status /= nf90_noerr ) then
        have_lhflx = .false.
      else
         call get_start_count(ncid, varid, scmlat, scmlon, ioptimeidx, strt4, cnt4)
         call wrap_get_vara_realx (ncid,varid,strt4,cnt4,lhflxobs)
         have_lhflx = .true.
      endif
   else
      call get_start_count(ncid, varid, scmlat, scmlon, ioptimeidx, strt4, cnt4)
      call wrap_get_vara_realx (ncid,varid,strt4,cnt4,lhflxobs)
      have_lhflx = .true.
   endif

   status = nf90_inq_varid( ncid, 'shflx', varid   )
   if ( status /= nf90_noerr ) then
      status = nf90_inq_varid( ncid, 'sh', varid   )
      if ( status /= nf90_noerr ) then
        have_shflx = .false.
      else
         call get_start_count(ncid, varid, scmlat, scmlon, ioptimeidx, strt4, cnt4)
        call wrap_get_vara_realx (ncid,varid,strt4,cnt4,shflxobs)
        have_shflx = .true.
      endif
   else
      call get_start_count(ncid, varid, scmlat, scmlon, ioptimeidx, strt4, cnt4)
      call wrap_get_vara_realx (ncid,varid,strt4,cnt4,shflxobs)
      have_shflx = .true.
   endif

   ! If REPLAY is used, then need to read in the global
   !   energy fixer
   status = nf90_inq_varid( ncid, 'heat_glob', varid   )
   if (status /= nf90_noerr) then
      have_heat_glob = .false.
   else
      call wrap_get_vara_realx (ncid,varid,strt4,cnt4,heat_glob_scm)
      have_heat_glob = .true.
   endif

!
!     fill in 3d forcing variables if we have both horizontal
!     and vertical components, but not the 3d
!
   if ( .not. have_cnst(1) .and. have_divq .and. have_vertdivq ) then
      do k=1,plev
         do m=1,pcnst
            divq3d(k,m) = divq(k,m) + vertdivq(k,m)
         enddo
      enddo
      have_divq3d = .true.
   endif

   if ( .not. have_divt3d .and. have_divt .and. have_vertdivt ) then
      if (masterproc) write(iulog,*) sub//'Don''t have divt3d - using divt and vertdivt'
      do k=1,plev
         divt3d(k) = divt(k) + vertdivt(k)
      enddo
      have_divt3d = .true.
   endif
!
!     make sure that use_3dfrc flag is set to true if we only have
!     3d forcing available
!
   if (scm_use_3dfrc) then
      if (have_divt3d .and. have_divq3d) then
         use_3dfrc = .true.
      else
         call endrun(sub//':ERROR :IOP file must have both divt3d and divq3d forcing when scm_use_3dfrc is set to .true.')
      endif
   endif

   status =  nf90_inq_varid( ncid, 'beta', varid   )
   if ( status /= nf90_noerr ) then
      betacam = 0._r8
   else
      call get_start_count(ncid, varid, scmlat, scmlon, ioptimeidx, strt4, cnt4)
      status = nf90_get_var(ncid, varid, srf(1), strt4)
      betacam=srf(1)
   endif

   status =  nf90_inq_varid( ncid, 'fixmas', varid   )
   if ( status /= nf90_noerr ) then
      fixmascam=1.0_r8
   else
      call get_start_count(ncid, varid, scmlat, scmlon, ioptimeidx, strt4, cnt4)
      status = nf90_get_var(ncid, varid, srf(1), strt4)
      fixmascam=srf(1)
   endif

   status = nf90_close( ncid )

   deallocate(dplevs)

end subroutine readiopdata

subroutine setiopupdate

!-----------------------------------------------------------------------
!
! Open and read netCDF file to extract time information
!
!---------------------------Code history--------------------------------
!
! Written by John Truesdale    August, 1996
!
!-----------------------------------------------------------------------
  implicit none

  character(len=*), parameter ::  sub = "setiopupdate"

!------------------------------Locals-----------------------------------

   integer :: next_date, next_sec
   integer :: ncsec,ncdate                      ! current time of day,date
   integer :: yr, mon, day                      ! year, month, and day component
!------------------------------------------------------------------------------

   call get_curr_date(yr,mon,day,ncsec)
   ncdate=yr*10000 + mon*100 + day

!------------------------------------------------------------------------------
!     Check if iop data needs to be updated and set doiopupdate accordingly
!------------------------------------------------------------------------------

   if ( is_first_step() ) then
      doiopupdate = .true.

   else

      call timemgr_time_inc(bdate, 0, next_date, next_sec, inc_s=tsec(iopTimeIdx+1))
      if ( ncdate > next_date .or. (ncdate == next_date &
         .and. ncsec >= next_sec)) then
         doiopupdate = .true.
         ! check to see if we need to move iopindex ahead more than 1 step
         do while ( ncdate > next_date .or. (ncdate == next_date .and. ncsec >= next_sec))
            iopTimeIdx = iopTimeIdx + 1
            call timemgr_time_inc(bdate, 0, next_date, next_sec, inc_s=tsec(iopTimeIdx+1))
         end do
#if DEBUG > 2
         if (masterproc) write(iulog,*) sub//'nstep = ',get_nstep()
         if (masterproc) write(iulog,*) sub//'ncdate=',ncdate,' ncsec=',ncsec
         if (masterproc) write(iulog,*) sub//'next_date=',next_date,' next_sec=',next_sec
         if (masterproc) write(iulog,*) sub//':******* do iop update'
#endif
      else
         doiopupdate = .false.
      end if
   endif                     ! if (endstep = 1 )
!
!     make sure we're
!     not going past end of iop data
!
   if ( ncdate > last_date .or. (ncdate == last_date &
      .and. ncsec > last_sec))  then
      call endrun(sub//':ERROR: Reached the end of the time varient dataset')
   endif

#if DEBUG > 1
   if (masterproc) write(iulog,*) sub//':iop time index = ' , ioptimeidx
#endif

end subroutine setiopupdate

!===============================================================================

subroutine plevs0 (nver, ps, ps0, hyam, hybm, hyai, hybi, pint    ,pmid    ,pdel)

!-----------------------------------------------------------------------
!
! Purpose:
! Define the pressures of the interfaces and midpoints from the
! coordinate definitions and the surface pressure.
!
! Author: B. Boville
!
!-----------------------------------------------------------------------
  implicit none


!-----------------------------------------------------------------------
  integer , intent(in)  :: nver         ! vertical dimension
  real(r8), intent(in)  :: ps           ! Surface pressure (pascals)
  real(r8), intent(in)  :: ps0          ! reference pressure (pascals)
  real(r8), intent(in)  :: hyam(plev)   ! hybrid midpoint coef
  real(r8), intent(in)  :: hybm(plev)   ! hybrid midpoint coef
  real(r8), intent(in)  :: hyai(plevp)  ! hybrid interface coef
  real(r8), intent(in)  :: hybi(plevp)  ! hybrid interface coef
  real(r8), intent(out) :: pint(nver+1) ! Pressure at model interfaces
  real(r8), intent(out) :: pmid(nver)   ! Pressure at model levels
  real(r8), intent(out) :: pdel(nver)   ! Layer thickness (pint(k+1) - pint(k))
!-----------------------------------------------------------------------

!---------------------------Local workspace-----------------------------
  integer :: k             ! Longitude, level indices
!-----------------------------------------------------------------------
!
! Set interface pressures
!
!$OMP PARALLEL DO PRIVATE (K)
  do k=1,nver+1
     pint(k) = hyai(k)*ps0 + hybi(k)*ps
  end do
!
! Set midpoint pressures and layer thicknesses
!
!$OMP PARALLEL DO PRIVATE (K)
  do k=1,nver
     pmid(k) = hyam(k)*ps0 + hybm(k)*ps
     pdel(k) = pint(k+1) - pint(k)
  end do

end subroutine plevs0

subroutine scmiop_flbc_inti ( co2vmr, ch4vmr, n2ovmr, f11vmr, f12vmr )
  !-----------------------------------------------------------------------
  !
  ! Purpose:
  ! Get start count for variable
  !
  !-----------------------------------------------------------------------

  implicit none

  real(r8), intent(out)  :: co2vmr, ch4vmr, n2ovmr, f11vmr, f12vmr

  !-----------------------------------------------------------------------

  co2vmr=co2vmrobs(1)
  ch4vmr=ch4vmrobs(1)
  n2ovmr=n2ovmrobs(1)
  f11vmr=f11vmrobs(1)
  f12vmr=f12vmrobs(1)
end subroutine scmiop_flbc_inti

subroutine get_start_count (ncid    ,varid  ,scmlat, scmlon, timeidx, start    ,count)

  !-----------------------------------------------------------------------
  !
  ! Purpose:
  ! set global lower boundary conditions
  !
  !-----------------------------------------------------------------------

  implicit none

  character(len=*), parameter ::  sub = "get_start_count"

!-----------------------------------------------------------------------
  integer , intent(in)    :: ncid         ! file id
  integer , intent(in)    :: varid        ! variable id
  integer , intent(in)    :: TimeIdx      ! time index
  real(r8), intent(in)    :: scmlat,scmlon! scm lat/lon
  integer , intent(out) :: start(:),count(:)

!---------------------------Local workspace-----------------------------
  integer            :: dims_set,nlev,var_ndims
  logical            :: usable_var
  character(len=cl)  :: dim_name
  integer            :: var_dimIDs( NF90_MAX_VAR_DIMS )
  real(r8)           :: closelat,closelon
  integer            :: latidx,lonidx,status,i
!-----------------------------------------------------------------------

   call shr_scam_GetCloseLatLon(ncid,scmlat,scmlon,closelat,closelon,latidx,lonidx)

   STATUS = NF90_INQUIRE_VARIABLE( NCID, varID, ndims=var_ndims )
!
!     surface variables
!
   if ( var_ndims == 0 ) then
      call endrun(sub//':ERROR: var_ndims is 0 for varid:',varid)
   endif

   STATUS = NF90_INQUIRE_VARIABLE( NCID, varID, dimids=var_dimIDs)
   if ( STATUS /= NF90_NOERR ) then
      write(iulog,* ) sub//'ERROR - Cant get dimension IDs for varid', varid
      call endrun(sub//':ERROR: Cant get dimension IDs for varid',varid)
   endif
!
!     Initialize the start and count arrays
!
   dims_set = 0
   nlev = 1
   do i =  var_ndims, 1, -1

      usable_var = .false.
      STATUS = NF90_INQUIRE_DIMENSION( NCID, var_dimIDs( i ), dim_name )

      if ( trim(dim_name) == 'lat' ) then
         start( i ) =  latIdx
         count( i ) = 1           ! Extract a single value
         dims_set = dims_set + 1
         usable_var = .true.
      endif

      if ( trim(dim_name) == 'lon' .or. trim(dim_name) == 'ncol' .or. trim(dim_name) == 'ncol_d' ) then
         start( i ) = lonIdx
         count( i ) = 1           ! Extract a single value
         dims_set = dims_set + 1
         usable_var = .true.
      endif

      if ( trim(dim_name) == 'lev' ) then
         STATUS = NF90_INQUIRE_DIMENSION( NCID, var_dimIDs( i ), len=nlev )
         start( i ) = 1
         count( i ) = nlev       ! Extract all levels
         dims_set = dims_set + 1
         usable_var = .true.
      endif

      if ( trim(dim_name) == 'ilev' ) then
         STATUS = NF90_INQUIRE_DIMENSION( NCID, var_dimIDs( i ), len=nlev )
         start( i ) = 1
         count( i ) = nlev        ! Extract all levels
         dims_set = dims_set + 1
         usable_var = .true.
      endif

      if ( trim(dim_name) == 'time' .OR. trim(dim_name) == 'tsec' ) then
         start( i ) = TimeIdx
         count( i ) = 1           ! Extract a single value
         dims_set = dims_set + 1
         usable_var = .true.
      endif
   end do
 end subroutine get_start_count

!=========================================================================
subroutine setiopupdate_init

!-----------------------------------------------------------------------
!
! Open and read netCDF file to extract time information
!   This subroutine should be called at the first SCM time step
!
!---------------------------Code history--------------------------------
!
! Written by John Truesdale    August, 1996
! Modified for E3SM by  Peter Bogenschutz 2017 - onward
!
!-----------------------------------------------------------------------
  implicit none

!------------------------------Locals-----------------------------------

   integer :: NCID,i
   integer :: tsec_varID, time_dimID
   integer :: bdate_varID
   integer :: STATUS
   integer :: next_date, next_sec
   integer :: ncsec,ncdate                      ! current time of day,date
   integer :: yr, mon, day                      ! year, month, and day component
   integer :: start_ymd,start_tod

   character(len=*), parameter ::  sub = "setiopupdate_init"
!!------------------------------------------------------------------------------

   ! Open and read pertinent information from the IOP file

   call handle_ncerr( nf90_open (iopfile, 0, ncid),&
        'ERROR - scamMod.F90:setiopupdate_init Failed to open iop file', __LINE__)

   ! Read time (tsec) variable

    STATUS = NF90_INQ_VARID( NCID, 'tsec', tsec_varID )
    if ( STATUS /= NF90_NOERR ) then
       write(iulog,*)sub//':ERROR: Cant get variable ID for tsec'
       STATUS = NF90_CLOSE ( NCID )
       call endrun(sub//':ERROR: Cant get variable ID for tsec')
    end if

    STATUS = NF90_INQ_VARID( NCID, 'bdate', bdate_varID )
    if ( STATUS /= NF90_NOERR ) then
       STATUS = NF90_INQ_VARID( NCID, 'basedate', bdate_varID )
       if ( STATUS /= NF90_NOERR ) then
          write(iulog,*)'ERROR - setiopupdate:Cant get variable ID for base date'
          STATUS = NF90_CLOSE ( NCID )
          call endrun(sub//':ERROR: Cant get variable ID for base date')
       endif
    endif

    STATUS = NF90_INQ_DIMID( NCID, 'time', time_dimID )
    if ( STATUS /= NF90_NOERR )  then
       STATUS = NF90_INQ_DIMID( NCID, 'tsec', time_dimID )
       if ( STATUS /= NF90_NOERR )  then
          write(iulog,* )'ERROR - setiopupdate:Could not find variable dim ID for time'
          STATUS = NF90_CLOSE ( NCID )
          call endrun(sub//':ERROR:Could not find variable dim ID for time')
       end if
    end if

    if ( STATUS /= NF90_NOERR )  &
       write(iulog,*)'ERROR - setiopupdate:Cant get variable dim ID for time'

    STATUS = NF90_INQUIRE_DIMENSION( NCID, time_dimID, len=ntime )
    if ( STATUS /= NF90_NOERR )then
       write(iulog,*)'ERROR - setiopupdate:Cant get time dimlen'
    endif

    if (.not.allocated(tsec)) allocate(tsec(ntime))

    STATUS = NF90_GET_VAR( NCID, tsec_varID, tsec )
    if ( STATUS /= NF90_NOERR )then
       write(iulog,*)'ERROR - setiopupdate:Cant get variable tsec'
    endif
    STATUS = NF90_GET_VAR( NCID, bdate_varID, bdate )
    if ( STATUS /= NF90_NOERR )then
       write(iulog,*)'ERROR - setiopupdate:Cant get variable bdate'
    endif

    ! Close the netCDF file
    STATUS = NF90_CLOSE( NCID )

    ! determine the last date in the iop dataset

    call timemgr_time_inc(bdate, 0, last_date, last_sec, inc_s=tsec(ntime))

    ! set the iop dataset index
    iopTimeIdx=0
    do i=1,ntime           ! set the first ioptimeidx
       call timemgr_time_inc(bdate, 0, next_date, next_sec, inc_s=tsec(i))
       call get_start_date(yr,mon,day,start_tod)
       start_ymd = yr*10000 + mon*100 + day

       if ( start_ymd > next_date .or. (start_ymd == next_date &
          .and. start_tod >= next_sec)) then
          iopTimeIdx = i
       endif
    enddo

    call get_curr_date(yr,mon,day,ncsec)
    ncdate=yr*10000 + mon*100 + day

    if (iopTimeIdx == 0.or.iopTimeIdx >= ntime) then
       call timemgr_time_inc(bdate, 0, next_date, next_sec, inc_s=tsec(1))
       write(iulog,*) 'Error::setiopupdate: Current model time does not fall within IOP period'
       write(iulog,*) ' Current CAM Date is ',ncdate,' and ',ncsec,' seconds'
       write(iulog,*) ' IOP start is        ',next_date,' and ',next_sec,'seconds'
       write(iulog,*) ' IOP end is          ',last_date,' and ',last_sec,'seconds'
       call endrun(sub//':ERROR: Current model time does not fall within IOP period')
    endif

    doiopupdate = .true.

end subroutine setiopupdate_init

end module scamMod
