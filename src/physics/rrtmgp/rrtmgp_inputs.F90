module rrtmgp_inputs
 use machine, only: kind_phys
 use ccpp_gas_optics_rrtmgp, only: ty_gas_optics_rrtmgp_ccpp
 use ccpp_optical_props,     only: ty_optical_props_1scl_ccpp, ty_optical_props_2str_ccpp
 use ccpp_gas_concentrations, only: ty_gas_concs_ccpp
 use ccpp_source_functions,   only: ty_source_func_lw_ccpp
 use string_utils,          only: to_lower
 use radiation_utils,       only: radiation_utils_init, get_sw_spectral_boundaries_ccpp

 implicit none
 private

 public :: rrtmgp_inputs_init
 public :: rrtmgp_inputs_run

 contains
!> \section arg_table_rrtmgp_inputs_init Argument Table
!! \htmlinclude rrtmgp_inputs_init.html
!!
  subroutine rrtmgp_inputs_init(ktopcam, ktoprad, nlaycam, sw_low_bounds, sw_high_bounds, nswbands,         &
                   pref_edge, nlay, pver, pverp, kdist_sw, kdist_lw, qrl, is_first_step, use_rad_dt_cosz,   &
                   timestep_size, nstep, iradsw, dt_avg, irad_always, is_first_restart_step,             &
                   nlwbands, nradgas, gasnamelength, iulog, idx_sw_diag, idx_nir_diag, idx_uv_diag,      &
                   idx_sw_cloudsim, idx_lw_diag, idx_lw_cloudsim, gaslist, nswgpts, nlwgpts, nlayp,      &
                   nextsw_cday, current_cal_day, band2gpt_sw, errmsg, errflg)

     ! Inputs
     integer,                       intent(in) :: nswbands
     integer,                       intent(in) :: pverp
     integer,                       intent(in) :: pver
     integer,                       intent(in) :: iradsw
     integer,                       intent(in) :: timestep_size
     integer,                       intent(in) :: nstep
     integer,                       intent(in) :: nlwbands
     integer,                       intent(in) :: nradgas
     integer,                       intent(in) :: iulog
     integer,                       intent(in) :: gasnamelength
     real(kind_phys),               intent(in) :: current_cal_day
     real(kind_phys), dimension(:), intent(in) :: pref_edge
     type(ty_gas_optics_rrtmgp_ccpp),    intent(in) :: kdist_sw
     type(ty_gas_optics_rrtmgp_ccpp),    intent(in) :: kdist_lw
     logical,                       intent(in) :: is_first_step
     logical,                       intent(in) :: is_first_restart_step
     logical,                       intent(in) :: use_rad_dt_cosz
     character(len=*),  dimension(:), intent(in) :: gaslist

     ! Outputs
     integer,                       intent(out) :: ktopcam
     integer,                       intent(out) :: ktoprad
     integer,                       intent(out) :: nlaycam
     integer,                       intent(out) :: nlay
     integer,                       intent(out) :: nlayp
     integer,                       intent(out) :: idx_sw_diag
     integer,                       intent(out) :: idx_nir_diag
     integer,                       intent(out) :: idx_uv_diag
     integer,                       intent(out) :: idx_sw_cloudsim
     integer,                       intent(out) :: idx_lw_diag
     integer,                       intent(out) :: idx_lw_cloudsim
     integer,                       intent(out) :: nswgpts
     integer,                       intent(out) :: nlwgpts
     integer, dimension(:,:),       intent(out) :: band2gpt_sw
     real(kind_phys),               intent(out) :: nextsw_cday
     real(kind_phys), dimension(:), intent(out) :: sw_low_bounds
     real(kind_phys), dimension(:), intent(out) :: sw_high_bounds
     real(kind_phys), dimension(:,:), intent(out) :: qrl
     character(len=*),                intent(out) :: errmsg
     integer,                         intent(out) :: errflg
     integer,                       intent(inout) :: irad_always
     real(kind_phys),               intent(inout) :: dt_avg        ! averaging time interval for zenith angle

     ! Local variables
     real(kind_phys), target :: wavenumber_low_shortwave(nswbands)
     real(kind_phys), target :: wavenumber_high_shortwave(nswbands)
     real(kind_phys), target :: wavenumber_low_longwave(nlwbands)
     real(kind_phys), target :: wavenumber_high_longwave(nlwbands)
     character(len=gasnamelength) :: gaslist_lc(nradgas)

     ! Set error variables
     errflg = 0
     errmsg = ''

     ! Number of layers in radiation calculation is capped by the number of
     ! pressure interfaces below 1 Pa.  When the entire model atmosphere is
     ! below 1 Pa then an extra layer is added to the top of the model for
     ! the purpose of the radiation calculation.
     nlay = count( pref_edge(:) > 1._kind_phys ) ! pascals (0.01 mbar)
     nlayp = nlay + 1

     if (nlay == pverp) then
        ! Top model interface is below 1 Pa.  RRTMGP is active in all model layers plus
        ! 1 extra layer between model top and 1 Pa.
        ktopcam = 1
        ktoprad = 2
        nlaycam = pver
     else if (nlay == (pverp-1)) then
        ! Special case nlay == (pverp-1) -- topmost interface outside bounds (CAM MT config), treat as if it is ok.
        ktopcam = 1
        ktoprad = 2
        nlaycam = pver
        nlay = nlay+1 ! reassign the value so later code understands to treat this case like nlay==pverp
        write(iulog,*) 'RADIATION_INIT: Special case of 1 model interface at p < 1Pa. Top layer will be INCLUDED in radiation calculation.'
        write(iulog,*) 'RADIATION_INIT: nlay = ',nlay, ' same as pverp: ',nlay==pverp
     else
        ! nlay < pverp.  nlay layers are used in radiation calcs, and they are
        ! all CAM layers.
        ktopcam = pver - nlay + 1
        ktoprad = 1
        nlaycam = nlay
     end if

     ! Set the sw/lw band boundaries in radconstants.  Also sets
     ! indicies of specific bands for diagnostic output and COSP input.
     call set_wavenumber_bands(kdist_sw, kdist_lw, nswbands, nlwbands, idx_sw_diag, idx_nir_diag, &
                 idx_uv_diag, idx_sw_cloudsim, idx_lw_diag, idx_lw_cloudsim, nswgpts, nlwgpts,    &
                 wavenumber_low_shortwave, wavenumber_high_shortwave, wavenumber_low_longwave,    &
                 wavenumber_high_longwave, band2gpt_sw, errmsg, errflg)
     if (errflg /= 0) then
        return
     end if

     call radiation_utils_init(nswbands, nlwbands, wavenumber_low_shortwave, wavenumber_high_shortwave, &
             wavenumber_low_longwave, wavenumber_high_longwave, errmsg, errflg)
     if (errflg /= 0) then
        return
     end if

     ! Initialize the SW band boundaries
     call get_sw_spectral_boundaries_ccpp(sw_low_bounds, sw_high_bounds, 'cm^-1', errmsg, errflg)
     if (errflg /= 0) then
        return
     end if

     if (is_first_step) then
        qrl = 0._kind_phys
     end if

     ! Set the radiation timestep for cosz calculations if requested using
     ! the adjusted iradsw value from radiation
     if (use_rad_dt_cosz)  then
        dt_avg = iradsw*timestep_size
     end if

     ! "irad_always" is number of time steps to execute radiation continuously from
     ! start of initial OR restart run
     if (irad_always > 0) then
        irad_always = irad_always + nstep
     end if

     ! Surface components to get radiation computed today
     if (.not. is_first_restart_step) then
        nextsw_cday = current_cal_day
     end if

  end subroutine rrtmgp_inputs_init

!> \section arg_table_rrtmgp_inputs_run Argument Table
!! \htmlinclude rrtmgp_inputs_run.html
!!
  subroutine rrtmgp_inputs_run(dosw, dolw, snow_associated, graupel_associated, &
                  pmid, pint, t, nday, idxday, cldfprime, &
                  coszrs, kdist_sw, t_sfc, emis_sfc, t_rad, pmid_rad,     &
                  pint_rad, t_day, pmid_day, pint_day, coszrs_day,        &
                  alb_dir, alb_dif, lwup, stebol, ncol, ktopcam, ktoprad, &
                  nswbands, asdir, asdif, sw_low_bounds, sw_high_bounds,  &
                  aldir, aldif, nlay, pverp, pver, cld, cldfsnow,         &
                  cldfgrau, graupel_in_rad, gasnamelength, gaslist,       &
                  gas_concs_lw, aer_lw, atm_optics_lw, kdist_lw,          &
                  sources_lw, aer_sw, atm_optics_sw, gas_concs_sw,        &
                  errmsg, errflg)
     ! Inputs
     logical, intent(in) :: graupel_in_rad
     integer, intent(in) :: ncol
     integer, intent(in) :: pver
     integer, intent(in) :: pverp
     integer, intent(in) :: nlay
     integer, intent(in) :: nswbands
     integer, intent(in) :: ktopcam
     integer, intent(in) :: ktoprad
     integer, intent(in) :: gasnamelength
     integer, intent(in) :: nday
     logical, intent(in) :: dosw
     logical, intent(in) :: dolw
     logical, intent(in) :: snow_associated
     logical, intent(in) :: graupel_associated
     integer, dimension(:), intent(in)           :: idxday
     real(kind_phys), dimension(:,:), intent(in) :: pmid
     real(kind_phys), dimension(:,:), intent(in) :: pint
     real(kind_phys), dimension(:,:), intent(in) :: t
     real(kind_phys), dimension(:,:), intent(in) :: cldfsnow
     real(kind_phys), dimension(:,:), intent(in) :: cldfgrau
     real(kind_phys), dimension(:,:), intent(in) :: cld
     real(kind_phys), dimension(:),   intent(in) :: sw_low_bounds
     real(kind_phys), dimension(:),   intent(in) :: sw_high_bounds
     real(kind_phys), dimension(:),   intent(in) :: coszrs
     real(kind_phys), dimension(:),   intent(in) :: lwup
     real(kind_phys), dimension(:),   intent(in) :: asdir
     real(kind_phys), dimension(:),   intent(in) :: asdif
     real(kind_phys), dimension(:),   intent(in) :: aldir
     real(kind_phys), dimension(:),   intent(in) :: aldif
     real(kind_phys),                 intent(in) :: stebol    ! stefan-boltzmann constant
     type(ty_gas_optics_rrtmgp_ccpp),     intent(in) :: kdist_sw  ! spectral information
     type(ty_gas_optics_rrtmgp_ccpp),     intent(in) :: kdist_lw  ! spectral information
     character(len=*), dimension(:),  intent(in) :: gaslist
     ! Outputs
     real(kind_phys), dimension(:,:), intent(out) :: t_rad
     real(kind_phys), dimension(:,:), intent(out) :: pmid_rad
     real(kind_phys), dimension(:,:), intent(out) :: pint_rad
     real(kind_phys), dimension(:,:), intent(out) :: t_day
     real(kind_phys), dimension(:,:), intent(out) :: pint_day
     real(kind_phys), dimension(:,:), intent(out) :: pmid_day
     real(kind_phys), dimension(:,:), intent(out) :: emis_sfc
     real(kind_phys), dimension(:,:), intent(out) :: alb_dir
     real(kind_phys), dimension(:,:), intent(out) :: alb_dif
     real(kind_phys), dimension(:,:), intent(out) :: cldfprime ! modiifed cloud fraciton

     real(kind_phys), dimension(:),   intent(out) :: t_sfc
     real(kind_phys), dimension(:),   intent(out) :: coszrs_day
     type(ty_gas_concs_ccpp),              intent(out) :: gas_concs_lw
     type(ty_optical_props_1scl_ccpp),     intent(out) :: atm_optics_lw
     type(ty_optical_props_1scl_ccpp),     intent(out) :: aer_lw
     type(ty_source_func_lw_ccpp),         intent(out) :: sources_lw
     type(ty_gas_concs_ccpp),              intent(out) :: gas_concs_sw
     type(ty_optical_props_2str_ccpp),     intent(out) :: atm_optics_sw
     type(ty_optical_props_2str_ccpp),     intent(out) :: aer_sw
     character(len=*), intent(out) :: errmsg
     integer,          intent(out) :: errflg

     ! Local variables
     real(kind_phys) :: tref_min
     real(kind_phys) :: tref_max
     integer :: idx, kdx, iband
     character(len=gasnamelength) :: gaslist_lc(size(gaslist))

     ! Set error variables
     errmsg = ''
     errflg = 0

     if (.not. dosw .and. .not. dolw) then
        return
     end if

     ! RRTMGP set state
     t_sfc = sqrt(sqrt(lwup(:ncol)/stebol))  ! Surface temp set based on longwave up flux.

     ! Set surface emissivity to 1.0.
     ! The land model *does* have its own surface emissivity, but is not spectrally resolved.
     ! The LW upward flux is calculated with that land emissivity, and the "radiative temperature"
     ! t_sfc is derived from that flux. We assume, therefore, that the emissivity is unity
     ! to be consistent with t_sfc.
     emis_sfc(:,:) = 1._kind_phys

     ! Level ordering is the same for both CAM and RRTMGP (top to bottom)
     t_rad(:,ktoprad:) = t(:ncol,ktopcam:)
     pmid_rad(:,ktoprad:) = pmid(:ncol,ktopcam:)
     pint_rad(:,ktoprad:) = pint(:ncol,ktopcam:)

     ! Add extra layer values if needed.
     if (nlay == pverp) then
        t_rad(:,1) = t(:ncol,1)
        ! The top reference pressure from the RRTMGP coefficients datasets is 1.005183574463 Pa
        ! Set the top of the extra layer just below that.
        pint_rad(:,1) = 1.01_kind_phys

        ! next interface down in LT will always be > 1Pa
        ! but in MT we apply adjustment to have it be 1.02 Pa if it was too high
        where (pint_rad(:,2) <= pint_rad(:,1)) pint_rad(:,2) = pint_rad(:,1)+0.01_kind_phys

        ! set the highest pmid (in the "extra layer") to the midpoint (guarantees > 1Pa)
        pmid_rad(:,1)   = pint_rad(:,1) + 0.5_kind_phys * (pint_rad(:,2) - pint_rad(:,1))

        ! For case of CAM MT, also ensure pint_rad(:,2) > pint_rad(:,1) & pmid_rad(:,2) > max(pmid_rad(:,1), min_pressure)
        where (pmid_rad(:,2) <= kdist_sw%get_press_min()) pmid_rad(:,2) = pint_rad(:,2) + 0.01_kind_phys
     else
        ! nlay < pverp, thus the 1 Pa level is within a CAM layer.  Assuming the top interface of
        ! this layer is at a pressure < 1 Pa, we need to adjust the top of this layer so that it
        ! is within the valid pressure range of RRTMGP (otherwise RRTMGP issues an error).  Then
        ! set the midpoint pressure halfway between the interfaces.
        pint_rad(:,1) = 1.01_kind_phys
        pmid_rad(:,1) = 0.5_kind_phys * (pint_rad(:,1) + pint_rad(:,2))
     end if

     ! Limit temperatures to be within the limits of RRTMGP validity.
     tref_min = kdist_sw%get_temp_min()
     tref_max = kdist_sw%get_temp_max()
     t_rad = merge(t_rad, tref_min, t_rad > tref_min)
     t_rad = merge(t_rad, tref_max, t_rad < tref_max)

     ! Construct arrays containing only daylight columns
     do idx = 1, nday
        t_day(idx,:)    = t_rad(idxday(idx),:)
        pmid_day(idx,:) = pmid_rad(idxday(idx),:)
        pint_day(idx,:) = pint_rad(idxday(idx),:)
        coszrs_day(idx) = coszrs(idxday(idx))
     end do
     ! Assign albedos to the daylight columns (from E3SM implementation)
     ! Albedos are imported from the surface models as broadband (visible, and near-IR),
     ! and we need to map these to appropriate narrower bands used in RRTMGP. Bands
     ! are categorized broadly as "visible/UV" or "infrared" based on wavenumber.
     ! Loop over bands, and determine for each band whether it is broadly in the
     ! visible or infrared part of the spectrum based on a dividing line of
     ! 0.7 micron, or 14286 cm^-1
     do iband = 1,nswbands
        if (is_visible(sw_low_bounds(iband)) .and. &
           is_visible(sw_high_bounds(iband))) then

           ! Entire band is in the visible
           do idx = 1, nday
              alb_dir(iband,idx) = asdir(idxday(idx))
              alb_dif(iband,idx) = asdif(idxday(idx))
           end do

        else if (.not.is_visible(sw_low_bounds(iband)) .and. &
                 .not.is_visible(sw_high_bounds(iband))) then
           ! Entire band is in the longwave (near-infrared)
           do idx = 1, nday
              alb_dir(iband,idx) = aldir(idxday(idx))
              alb_dif(iband,idx) = aldif(idxday(idx))
           end do
        else
           ! Band straddles the visible to near-infrared transition, so we take
           ! the albedo to be the average of the visible and near-infrared
           ! broadband albedos
           do idx = 1, nday
              alb_dir(iband,idx) = 0.5_kind_phys * (aldir(idxday(idx)) + asdir(idxday(idx)))
              alb_dif(iband,idx) = 0.5_kind_phys * (aldif(idxday(idx)) + asdif(idxday(idx)))
           end do
        end if
     end do
     ! Strictly enforce albedo bounds
     where (alb_dir < 0)
        alb_dir = 0.0_kind_phys
     end where
     where (alb_dir > 1)
        alb_dir = 1.0_kind_phys
     end where
     where (alb_dif < 0)
        alb_dif = 0.0_kind_phys
     end where
     where (alb_dif > 1)
        alb_dif = 1.0_kind_phys
     end where

     ! modified cloud fraction
     ! Compute modified cloud fraction, cldfprime.
     ! 1. initialize as cld
     ! 2. modify for snow. use max(cld, cldfsnow)
     ! 3. modify for graupel if graupel_in_rad is true.
     !    use max(cldfprime, cldfgrau)
     if (snow_associated) then
        do kdx = 1, pver
           do idx = 1, ncol
              cldfprime(idx,kdx) = max(cld(idx,kdx), cldfsnow(idx,kdx))
           end do
        end do
     else
        cldfprime(:ncol,:) = cld(:ncol,:)
     end if

     if (graupel_associated .and. graupel_in_rad) then
        do kdx = 1, pver
           do idx = 1, ncol
              cldfprime(idx,kdx) = max(cldfprime(idx,kdx), cldfgrau(idx,kdx))
           end do
        end do
     end if

     ! Create lowercase version of the gaslist for RRTMGP.  The ty_gas_concs_ccpp objects
     ! work with CAM's uppercase names, but other objects that get input from the gas
     ! concs objects don't work.
     do idx = 1, size(gaslist)
        gaslist_lc(idx) = to_lower(gaslist(idx))
     end do

     ! If no daylight columns, can't create empty RRTMGP objects
     if (dosw .and. nday > 0) then
        ! Initialize object for gas concentrations.
         errmsg = gas_concs_sw%gas_concs%init(gaslist_lc)
        if (len_trim(errmsg) > 0) then
           errflg = 1
           return
        end if

        ! Initialize object for combined gas + aerosol + cloud optics.
        ! Allocates arrays for properties represented on g-points.
        errmsg = atm_optics_sw%alloc_2str(nday, nlay, kdist_sw)
        if (len_trim(errmsg) > 0) then
           errflg = 1
           return
        end if

        ! Initialize object for SW aerosol optics.  Allocates arrays
        ! for properties represented by band.
        errmsg = aer_sw%alloc_2str(nday, nlay, kdist_sw%get_band_lims_wavenumber())
        if (len_trim(errmsg) > 0) then
           errflg = 1
           return
        end if
     end if

     if (dolw) then
        ! Initialize object for gas concentrations
        errmsg = gas_concs_lw%gas_concs%init(gaslist_lc)
        if (len_trim(errmsg) > 0) then
           errflg = 1
           return
        end if

        ! Initialize object for combined gas + aerosol + cloud optics.
        errmsg = atm_optics_lw%alloc_1scl(ncol, nlay, kdist_lw)
        if (len_trim(errmsg) > 0) then
           errflg = 1
           return
        end if

        ! Initialize object for LW aerosol optics.
        errmsg = aer_lw%alloc_1scl(ncol, nlay, kdist_lw%get_band_lims_wavenumber())
        if (len_trim(errmsg) > 0) then
           errflg = 1
           return
        end if

        ! Initialize object for Planck sources.
        errmsg = sources_lw%sources%alloc(ncol, nlay, kdist_lw)
        if (len_trim(errmsg) > 0) then
           errflg = 1
           return
        end if
     end if

  end subroutine rrtmgp_inputs_run

!=========================================================================================
!                                     HELPER FUNCTIONS                                   !
!=========================================================================================
  subroutine set_wavenumber_bands(kdist_sw, kdist_lw, nswbands, nlwbands, idx_sw_diag, idx_nir_diag, &
                  idx_uv_diag, idx_sw_cloudsim, idx_lw_diag, idx_lw_cloudsim, nswgpts, nlwgpts,      &
                  wavenumber_low_shortwave, wavenumber_high_shortwave, wavenumber_low_longwave,      &
                  wavenumber_high_longwave, band2gpt_sw, errmsg, errflg)
   ! Set the low and high limits of the wavenumber grid for sw and lw.
   ! Values come from RRTMGP coefficients datasets, and are stored in the
   ! kdist objects.
   !
   ! Set band indices for bands containing specific wavelengths.

   ! Arguments
   type(ty_gas_optics_rrtmgp_ccpp), intent(in) :: kdist_sw
   type(ty_gas_optics_rrtmgp_ccpp), intent(in) :: kdist_lw
   integer,                    intent(in) :: nswbands
   integer,                    intent(in) :: nlwbands

   integer,                   intent(out) :: idx_sw_diag
   integer,                   intent(out) :: idx_nir_diag
   integer,                   intent(out) :: idx_uv_diag
   integer,                   intent(out) :: idx_sw_cloudsim
   integer,                   intent(out) :: idx_lw_diag
   integer,                   intent(out) :: idx_lw_cloudsim
   integer,                   intent(out) :: nswgpts
   integer,                   intent(out) :: nlwgpts
   integer, dimension(:,:),   intent(out) :: band2gpt_sw
   real(kind_phys), dimension(:), intent(out) :: wavenumber_low_shortwave
   real(kind_phys), dimension(:), intent(out) :: wavenumber_high_shortwave
   real(kind_phys), dimension(:), intent(out) :: wavenumber_low_longwave
   real(kind_phys), dimension(:), intent(out) :: wavenumber_high_longwave
   character(len=*),          intent(out) :: errmsg
   integer,                   intent(out) :: errflg

   ! Local variables
   integer :: istat
   real(kind_phys), allocatable :: values(:,:)

   character(len=*), parameter :: sub = 'set_wavenumber_bands'
   !----------------------------------------------------------------------------

   ! Initialize error variables
   errflg = 0
   errmsg = ''
   ! Check that number of sw/lw bands in gas optics files matches the parameters.
   if (kdist_sw%get_nband() /= nswbands) then
      write(errmsg,'(a, a,i4,a,i4)') sub, ': ERROR: number of sw bands in file, ', kdist_sw%get_nband(), &
         ", doesn't match parameter nswbands= ", nswbands
      errflg = 1
      return
   end if
   if (kdist_lw%get_nband() /= nlwbands) then
      write(errmsg,'(a, a,i4,a,i4)') sub, ': ERROR: number of lw bands in file, ', kdist_lw%get_nband(), &
         ", doesn't match parameter nlwbands= ", nlwbands
      errflg = 1
      return
   end if

   nswgpts = kdist_sw%get_ngpt()
   nlwgpts = kdist_lw%get_ngpt()

   ! SW band bounds in cm^-1
   allocate( values(2,nswbands), stat=istat )
   if (istat/=0) then
      write(errmsg, '(a,a)') sub, ': ERROR allocating array: values(2,nswbands)'
      errflg = 1
      return
   end if
   values = kdist_sw%get_band_lims_wavenumber()
   wavenumber_low_shortwave = values(1,:)
   wavenumber_high_shortwave = values(2,:)

   ! First and last g-point for each SW band:
   band2gpt_sw = kdist_sw%get_band_lims_gpoint()

   ! Indices into specific bands
   call get_band_index_by_value('sw', 500.0_kind_phys, 'nm', nswbands, &
        wavenumber_low_shortwave, wavenumber_high_shortwave, idx_sw_diag, errmsg, errflg)
   if (errflg /= 0) then
      return
   end if
   call get_band_index_by_value('sw', 1000.0_kind_phys, 'nm', nswbands, &
        wavenumber_low_shortwave, wavenumber_high_shortwave, idx_nir_diag, errmsg, errflg)
   if (errflg /= 0) then
      return
   end if
   call get_band_index_by_value('sw', 400._kind_phys, 'nm', nswbands, &
        wavenumber_low_shortwave, wavenumber_high_shortwave, idx_uv_diag, errmsg, errflg)
   if (errflg /= 0) then
      return
   end if
   call get_band_index_by_value('sw', 0.67_kind_phys, 'micron', nswbands, &
        wavenumber_low_shortwave, wavenumber_high_shortwave, idx_sw_cloudsim, errmsg, errflg)
   if (errflg /= 0) then
      return
   end if

   deallocate(values)

   ! LW band bounds in cm^-1
   allocate( values(2,nlwbands), stat=istat )
   if (istat/=0) then
      write(errmsg, '(a,a)') sub, ': ERROR allocating array: values(2,nlwbands)'
      errflg = 1
      return
   end if
   values = kdist_lw%get_band_lims_wavenumber()
   wavenumber_low_longwave = values(1,:)
   wavenumber_high_longwave = values(2,:)

   ! Indices into specific bands
   call get_band_index_by_value('lw', 1000.0_kind_phys, 'cm^-1', nlwbands, &
        wavenumber_low_longwave, wavenumber_high_longwave, idx_lw_diag, errmsg, errflg)
   if (errflg /= 0) then
      return
   end if
   call get_band_index_by_value('lw', 10.5_kind_phys, 'micron', nlwbands, &
        wavenumber_low_longwave, wavenumber_high_longwave, idx_lw_cloudsim, errmsg, errflg)
   if (errflg /= 0) then
      return
   end if

 end subroutine set_wavenumber_bands

!=========================================================================================

 subroutine get_band_index_by_value(swlw, targetvalue, units, nbnds, wavenumber_low, &
                wavenumber_high, ans, errmsg, errflg)

   ! Find band index for requested wavelength/wavenumber.

   character(len=*), intent(in) :: swlw        ! sw or lw bands
   real(kind_phys),         intent(in) :: targetvalue  
   character(len=*), intent(in) :: units       ! units of targetvalue
   integer,          intent(in) :: nbnds
   real(kind_phys), target, dimension(:), intent(in) :: wavenumber_low
   real(kind_phys), target, dimension(:), intent(in) :: wavenumber_high
   character(len=*), intent(out) :: errmsg
   integer,          intent(out) :: errflg
   integer,          intent(out) :: ans

   ! local
   real(kind_phys), pointer, dimension(:) :: lowboundaries, highboundaries
   real(kind_phys) :: tgt
   integer  :: idx

   character(len=*), parameter :: sub = 'get_band_index_by_value'
   !----------------------------------------------------------------------------

   ! Initialize error variables
   errflg = 0
   errmsg = ''
   lowboundaries => wavenumber_low
   highboundaries => wavenumber_high
   if (trim(swlw) /= 'sw' .and. trim(swlw) /= 'lw') then
      write(errmsg,'(a,a)') 'rrtmgp_inputs: get_band_index_by_value: type of bands not recognized: ', swlw
      errflg = 1
      return
   end if

   ! band info is in cm^-1 but target value may be other units,
   ! so convert targetvalue to cm^-1
   select case (units)
   case ('inv_cm','cm^-1','cm-1')
      tgt = targetvalue
   case('m','meter','meters')
      tgt = 1.0_kind_phys / (targetvalue * 1.e2_kind_phys)
   case('nm','nanometer','nanometers')
      tgt = 1.0_kind_phys / (targetvalue * 1.e-7_kind_phys)
   case('um','micrometer','micrometers','micron','microns')
      tgt = 1.0_kind_phys / (targetvalue * 1.e-4_kind_phys)
   case('cm','centimeter','centimeters')
      tgt = 1._kind_phys/targetvalue
   case default
      write(errmsg,'(a,a)') 'rrtmgp_inputs: get_band_index_by_value: units not recognized: ', units
      errflg = 1
   end select

   ! now just loop through the array
   ans = 0
   do idx = 1,nbnds
      if ((tgt > lowboundaries(idx)) .and. (tgt <= highboundaries(idx))) then
         ans = idx
         exit
      end if
   end do

   if (ans == 0) then
      write(errmsg,'(f10.3,a,a)') targetvalue, ' ', trim(units)
      errflg = 1
   end if
   
 end subroutine get_band_index_by_value

 !=========================================================================================

 pure logical function is_visible(wavenumber)

   ! Wavenumber is in the visible if it is above the visible threshold
   ! wavenumber, and in the infrared if it is below the threshold
   ! This function doesn't distinquish between visible and UV.

   ! wavenumber in inverse cm (cm^-1)
   real(kind_phys), intent(in) :: wavenumber

   ! Set threshold between visible and infrared to 0.7 micron, or 14286 cm^-1
   real(kind_phys), parameter :: visible_wavenumber_threshold = 14286._kind_phys  ! cm^-1

   if (wavenumber > visible_wavenumber_threshold) then
      is_visible = .true.
   else
      is_visible = .false.
   end if

 end function is_visible

end module rrtmgp_inputs
