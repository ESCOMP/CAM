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

use shr_kind_mod,   only: r8 => shr_kind_r8
use pmgrid,         only: plon, plat, plev, plevp
use constituents,   only: pcnst
use shr_scam_mod,   only: shr_scam_getCloseLatLon
use dycore,         only: dycore_is
use cam_logfile,    only: iulog
use cam_abortutils, only: endrun

implicit none
private

! PUBLIC INTERFACES:

public scam_readnl   ! read SCAM namelist options 

! PUBLIC MODULE DATA:

real(r8), public ::  pressure_levels(plev)
real(r8), public ::  scmlat   ! input namelist latitude for scam
real(r8), public ::  scmlon   ! input namelist longitude for scam


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

character*(max_path_len), public ::  modelfile
character*(max_path_len), public ::  analysisfile
character*(max_path_len), public ::  sicfile
character*(max_path_len), public ::  userfile
character*(max_path_len), public ::  sstfile
character*(max_path_len), public ::  lsmpftfile
character*(max_path_len), public ::  pressfile
character*(max_path_len), public ::  topofile
character*(max_path_len), public ::  ozonefile
character*(max_path_len), public ::  iopfile
character*(max_path_len), public ::  absemsfile
character*(max_path_len), public ::  aermassfile
character*(max_path_len), public ::  aeropticsfile
character*(max_path_len), public ::  timeinvfile
character*(max_path_len), public ::  lsmsurffile
character*(max_path_len), public ::  lsminifile

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
real(r8), public ::      shflxobs(1)         ! observed surface sensible heat flux
real(r8), public ::      q1obs(plev)         ! observed apparent heat source
real(r8), public ::      q2obs(plev)         ! observed apparent heat sink
real(r8), public ::      tdiff(plev)         ! model minus observed temp 
real(r8), public ::      tground(1)          ! ground temperature
real(r8), public ::      tobs(plev)          ! actual temperature
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
logical, public ::  scm_use_obs_T          = .false. ! Use the SCAM-IOP specified observed T   at each time step instead of forecasting.
logical, public ::  scm_force_latlon       = .false. ! force scam to use the lat lon fields specified in the scam namelist not what is closest to iop avail lat lon
real*8, public              ::  scm_relax_top_p       = 1.e36_r8 ! upper bound for scm relaxation
real*8, public              ::  scm_relax_bot_p       = -1.e36_r8 !  lower bound for scm relaxation
real*8, public              ::  scm_relax_tau_sec       = 10800._r8  ! relaxation time constant (sec)

! +++BPM:
! modification... allow a linear ramp in relaxation time scale:
logical, public :: scm_relax_linear = .false.
real*8, public    :: scm_relax_tau_bot_sec = 10800._r8
real*8, public    :: scm_relax_tau_top_sec = 10800._r8
character(len=26), public  :: scm_relax_fincl(pcnst)

!
! note that scm_use_obs_uv is set to true to be consistent with CAM BFB testing
!

logical, public ::  scm_use_obs_uv         = .true. ! Use the SCAM-IOP specified observed u,v at each time step instead of forecasting.

logical, public ::  scm_use_obs_qv         = .false. ! Use the SCAM-IOP specified observed qv  at each time step instead of forecasting.
logical, public ::  scm_iop_lhflxshflxTg   = .false. !turn off LW rad
logical, public ::  scm_iop_Tg             = .false. !turn off LW rad

character(len=200), public ::  scm_clubb_iop_name   ! IOP name for CLUBB

!=======================================================================
contains
!=======================================================================

subroutine scam_readnl(nlfile,single_column_in,scmlat_in,scmlon_in)

  use namelist_utils,  only: find_group_name
  use units,           only: getunit, freeunit
  use dycore,          only: dycore_is
  use wrap_nf,         only: wrap_open
  use spmd_utils,      only : masterproc,npes
  use netcdf,          only : nf90_inquire_attribute,NF90_NOERR,NF90_GLOBAL,NF90_NOWRITE


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
  integer  :: latidx, lonidx
  logical  :: adv
  real(r8) :: ioplat,ioplon

! this list should include any variable that you might want to include in the namelist
  namelist /scam_nl/ iopfile, scm_iop_lhflxshflxTg, scm_iop_Tg, scm_relaxation, &
       scm_relax_top_p,scm_relax_bot_p,scm_relax_tau_sec, &
       scm_cambfb_mode,scm_crm_mode,scm_zadv_uv,scm_zadv_T,scm_zadv_q,&
       scm_use_obs_T, scm_use_obs_uv, scm_use_obs_qv, &
       scm_relax_linear, scm_relax_tau_top_sec, &
       scm_relax_tau_bot_sec, scm_force_latlon, scm_relax_fincl, scm_backfill_iop_w_init

  single_column=single_column_in

  iopfile            = ' '
  scm_clubb_iop_name = ' '
  scm_relax_fincl(:) = ' '
  
  if( single_column ) then
     if( npes.gt.1) call endrun('SCAM_READNL: SCAM doesnt support using more than 1 pe.')

     if (.not. dycore_is('EUL') .or. plon /= 1 .or. plat /=1 ) then 
        call endrun('SCAM_SETOPTS: must compile model for SCAM mode when namelist parameter single_column is .true.')
     endif

     scmlat=scmlat_in
     scmlon=scmlon_in
     
     if( scmlat .lt. -90._r8 .or. scmlat .gt. 90._r8 ) then
        call endrun('SCAM_READNL: SCMLAT must be between -90. and 90. degrees.')
     elseif( scmlon .lt. 0._r8 .or. scmlon .gt. 360._r8 ) then
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
     if( iopfile .ne. "" ) then 
        use_iop = .true.
     else
        call endrun('SCAM_READNL: must specify IOP file for single column mode')
     endif

     call wrap_open( iopfile, NF90_NOWRITE, ncid )

     if( nf90_inquire_attribute( ncid, NF90_GLOBAL, 'CAM_GENERATED_FORCING', iatt ) .EQ. NF90_NOERR ) then
        use_camiop = .true.
     else
        use_camiop = .false.
     endif
     
     ! If we are not forcing the lat and lon from the namelist use the closest lat and lon that is found in the IOP file.
     if (.not.scm_force_latlon) then
        call shr_scam_GetCloseLatLon( ncid, scmlat, scmlon, ioplat, ioplon, latidx, lonidx )
        write(iulog,*) 'SCAM_READNL: using closest IOP column to lat/lon specified in drv_in'
        write(iulog,*) '   requested lat,lon    =',scmlat,', ',scmlon
        write(iulog,*) '   closest IOP lat,lon  =',ioplat,', ',ioplon
     
        scmlat = ioplat
        scmlon = ioplon
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
        write (iulog,*) '  scm_use_obs_qv              = ',scm_use_obs_qv
        write (iulog,*) '  scm_use_obs_uv              = ',scm_use_obs_uv
        write (iulog,*) '  scm_zadv_T                  = ',trim(scm_zadv_T)
        write (iulog,*) '  scm_zadv_q                  = ',trim(scm_zadv_q)
        write (iulog,*) '  scm_zadv_uv                 = ',trim(scm_zadv_uv)
        write (iulog,*) '  scm_relax_finc: '
        ! output scm_relax_fincl character array
        do i=1,pcnst
           if (scm_relax_fincl(i) .ne. '') then
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

!===============================================================================

end module scamMod
