!! This module handles reading the namelist and provides access to some other flags
!! that control CARMA's behavior.
!!
!! It needs to be in its own file to resolve some circular dependencies.
!!
!! @author  Chuck Bardeen
!! @version Aug-2010
module carma_flags_mod

  use shr_kind_mod,   only: r8 => shr_kind_r8
  use spmd_utils,     only: masterproc

  ! Flags for integration with CAM Microphysics

  implicit none
  public

  integer, parameter           :: carma_maxdiags = 100
  integer, protected           :: carma_ndiagpkgs  ! Number of diags_packages listed
  integer, protected           :: carma_ndebugpkgs  ! Number of diags_packages listed

  ! Namelist flags
  !
  ! NOTE: Setting the carma_flag to false prevents CARMA from doing any microphysics
  ! calculations, but it will still initialize itself. This allows the same build and
  ! namelist to be used, but the CARMA processing diabled. Use the configure option
  ! -carma none to totally disable CARMA and prevent even the register from happening.
  logical, protected           :: carma_flag        = .false.   ! If .true. then turn on CARMA microphysics in CAM
  logical, protected           :: carma_do_aerosol  = .true.    ! If .true. then CARMA is processed after surface coupling
  logical, protected           :: carma_do_coremasscheck = .false. ! If .true. then do coremasscheck and abort model after certain subroutines
  logical, protected           :: carma_do_cldice   = .false.   ! If .true. then do cloud ice
  logical, protected           :: carma_do_cldliq   = .false.   ! If .true. then do cloud liquid
  logical, protected           :: carma_do_clearsky = .false.   ! If .true. then do clear sky particle calculations
  logical, protected           :: carma_do_cloudborne = .false. ! If .true. then do then the carma groups can be cloudborne
  logical, protected           :: carma_do_coag     = .false.   ! If .true. then do coagulation
  logical, protected           :: carma_do_detrain  = .false.   ! If .true. then do detrain
  logical, protected           :: carma_do_drydep   = .false.   ! If .true. then do dry deposition
  logical, protected           :: carma_do_emission = .false.   ! If .true. then do emission
  logical, protected           :: carma_do_fixedinit= .false.   ! If .true. then do fixed initialization to a reference state
  logical, protected           :: carma_hetchem_feedback=.false.! If .true. then CARMA sulfate surface area density used in heterogeneous chemistry
  logical, protected           :: carma_rad_feedback= .false.   ! If .true. then CARMA sulfate mass mixing ratio & effective radius used in radiation
  logical, protected           :: carma_do_explised = .false.   ! If .true. then do sedimentation with substepping
  logical, protected           :: carma_do_incloud  = .false.   ! If .true. then do incloud particle calculations
  logical, protected           :: carma_do_budget_diags  = .false.   ! If .true. then do budget diagnostics
  logical, protected           :: carma_do_package_diags = .false.   ! If .true. then do package diagnostics
  logical, protected           :: carma_do_grow     = .false.   ! If .true. then do growth
  logical, protected           :: carma_do_optics   = .false.   ! If .true. then do optical properties file
  logical, protected           :: carma_do_partialinit= .false. ! If .true. then do initialization of coagulation to a reference state (requires fixedinit)
  logical, protected           :: carma_do_pheat    = .false.   ! If .true. then do particle heating
  logical, protected           :: carma_do_pheatatm = .false.   ! If .true. then do particle heating of atmosphere
  logical, protected           :: carma_do_substep  = .false.   ! If .true. then do substeping
  logical, protected           :: carma_do_thermo   = .false.   ! If .true. then do solve thermodynamics equation
  logical, protected           :: carma_do_wetdep   = .false.   ! If .true. then do wet deposition
  logical, protected           :: carma_do_vdiff    = .false.   ! If .true. then do vertical brownian diffusion
  logical, protected           :: carma_do_vtran    = .false.   ! If .true. then do vertical transport
  integer, protected           :: carma_diags_file  = 0         ! Default file for diagnostic output
  integer, protected           :: carma_maxsubsteps = 1         ! Maximum number of time substeps allowed
  integer, protected           :: carma_minsubsteps = 1         ! Minimum number of time substeps allowed
  integer, protected           :: carma_maxretries  = 8         ! Maximum number of time substeps allowed
  real(r8), protected          :: carma_conmax      = 0.1_r8    ! Minumum relative concentration to consider in substep
  real(r8), protected          :: carma_dgc_threshold  = 0.0_r8 ! When non-zero, the largest percentage change in gas concentration allowed per substep.
  real(r8), protected          :: carma_ds_threshold  = 0.0_r8  ! When non-zero, the largest percentage change in gas saturation allowed per substep.
  real(r8), protected          :: carma_dt_threshold  = 0.0_r8  ! When non-zero, the largest change in temperature (K) allowed per substep.
  real(r8), protected          :: carma_tstick      = 1.0_r8    ! Thermal accommodation coefficient
  real(r8), protected          :: carma_gsticki     = 0.93_r8   ! Growth accommodation coefficient for ice
  real(r8), protected          :: carma_gstickl     = 1.0_r8    ! Growth accommodation coefficient for liquid
  real(r8), protected          :: carma_cstick      = 1.0_r8    ! Coagulation accommodation coefficient
  real(r8), protected          :: carma_rhcrit      = 1.0_r8    ! Critical relative humidity for liquid clouds
  real(r8), protected          :: carma_vf_const    = 0.0_r8    ! If specified and non-zero, constant fall velocity for all particles [cm/s]
  character(len=32), protected :: carma_model       = "none"    ! String (no spaces) that identifies the model
  character(len=10), protected :: carma_sulfnuc_method = "none" ! Sulfate Nucleation method
  character(len=32), protected :: carma_diags_packages(carma_maxdiags) = " " ! Names of physics packages for which diagnostic output is desired
  character(len=12), protected :: carma_debug_packages(carma_maxdiags) = " " ! Names of physics packages for which debug output is desired


contains


  !! Read the CARMA runtime options from the namelist
  !!
  !! @author  Chuck Bardeen
  !! @version Aug-2010
  subroutine carma_readnl(nlfile)

    ! Read carma namelist group.

    use cam_abortutils,  only: endrun
    use namelist_utils,  only: find_group_name
    use spmd_utils,      only: mpicom, masterprocid, mpi_real8, mpi_integer, mpi_logical, mpi_character
    use carma_model_flags_mod, only: carma_model_readnl

    ! args

    character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

    ! local vars

    integer :: unitn, ierr, i

    ! read namelist for CARMA
    namelist /carma_nl/ &
      carma_flag, &
      carma_do_aerosol, &
      carma_do_coremasscheck, &
      carma_do_cldliq, &
      carma_do_cldice, &
      carma_do_clearsky, &
      carma_do_cloudborne, &
      carma_do_coag, &
      carma_do_detrain, &
      carma_do_drydep, &
      carma_do_emission, &
      carma_do_fixedinit, &
      carma_hetchem_feedback, &
      carma_rad_feedback, &
      carma_do_explised, &
      carma_do_incloud, &
      carma_do_budget_diags, &
      carma_do_package_diags, &
      carma_do_grow, &
      carma_do_optics, &
      carma_do_partialinit, &
      carma_do_pheat, &
      carma_do_pheatatm, &
      carma_do_substep, &
      carma_do_thermo, &
      carma_do_wetdep, &
      carma_do_vdiff, &
      carma_do_vtran, &
      carma_maxsubsteps, &
      carma_minsubsteps, &
      carma_maxretries, &
      carma_model, &
      carma_conmax, &
      carma_dgc_threshold, &
      carma_ds_threshold, &
      carma_dt_threshold, &
      carma_tstick, &
      carma_gsticki, &
      carma_gstickl, &
      carma_cstick, &
      carma_rhcrit, &
      carma_vf_const, &
      carma_sulfnuc_method, &
      carma_do_budget_diags, &
      carma_do_package_diags, &
      carma_diags_packages, &
      carma_debug_packages, &
      carma_diags_file

    if (masterproc) then
       open( newunit=unitn, file=trim(nlfile), status='old' )
       call find_group_name(unitn, 'carma_nl', status=ierr)
       if (ierr == 0) then
          read(unitn, carma_nl, iostat=ierr)
          if (ierr /= 0) then
             call endrun('carma_readnl: ERROR reading namelist')
          end if
       end if
       close(unitn)
    end if

    call mpi_bcast (carma_flag,            1 ,mpi_logical, masterprocid, mpicom, ierr)
    call mpi_bcast (carma_do_aerosol,      1 ,mpi_logical, masterprocid, mpicom, ierr)
    call mpi_bcast (carma_do_coremasscheck,1 ,mpi_logical, masterprocid, mpicom, ierr)
    call mpi_bcast (carma_do_cldliq,       1 ,mpi_logical, masterprocid, mpicom, ierr)
    call mpi_bcast (carma_do_cldice,       1 ,mpi_logical, masterprocid, mpicom, ierr)
    call mpi_bcast (carma_do_clearsky,     1 ,mpi_logical, masterprocid, mpicom, ierr)
    call mpi_bcast (carma_do_cloudborne,   1 ,mpi_logical, masterprocid, mpicom, ierr)
    call mpi_bcast (carma_do_coag,         1 ,mpi_logical, masterprocid, mpicom, ierr)
    call mpi_bcast (carma_do_detrain,      1 ,mpi_logical, masterprocid, mpicom, ierr)
    call mpi_bcast (carma_do_drydep,       1 ,mpi_logical, masterprocid, mpicom, ierr)
    call mpi_bcast (carma_do_emission,     1 ,mpi_logical, masterprocid, mpicom, ierr)
    call mpi_bcast (carma_do_fixedinit,    1 ,mpi_logical, masterprocid, mpicom, ierr)
    call mpi_bcast (carma_hetchem_feedback,1 ,mpi_logical, masterprocid, mpicom, ierr)
    call mpi_bcast (carma_rad_feedback,    1 ,mpi_logical, masterprocid, mpicom, ierr)
    call mpi_bcast (carma_do_explised,     1 ,mpi_logical, masterprocid, mpicom, ierr)
    call mpi_bcast (carma_do_budget_diags, 1 ,mpi_logical, masterprocid, mpicom, ierr)
    call mpi_bcast (carma_do_package_diags,1 ,mpi_logical, masterprocid, mpicom, ierr)
    call mpi_bcast (carma_do_incloud,      1 ,mpi_logical, masterprocid, mpicom, ierr)
    call mpi_bcast (carma_do_grow,         1 ,mpi_logical, masterprocid, mpicom, ierr)
    call mpi_bcast (carma_do_optics,       1 ,mpi_logical, masterprocid, mpicom, ierr)
    call mpi_bcast (carma_do_partialinit,  1 ,mpi_logical, masterprocid, mpicom, ierr)
    call mpi_bcast (carma_do_pheat,        1 ,mpi_logical, masterprocid, mpicom, ierr)
    call mpi_bcast (carma_do_pheatatm,     1 ,mpi_logical, masterprocid, mpicom, ierr)
    call mpi_bcast (carma_do_substep,      1 ,mpi_logical, masterprocid, mpicom, ierr)
    call mpi_bcast (carma_do_thermo,       1 ,mpi_logical, masterprocid, mpicom, ierr)
    call mpi_bcast (carma_do_wetdep,       1 ,mpi_logical, masterprocid, mpicom, ierr)
    call mpi_bcast (carma_do_vdiff,        1 ,mpi_logical, masterprocid, mpicom, ierr)
    call mpi_bcast (carma_do_vtran,        1 ,mpi_logical, masterprocid, mpicom, ierr)
    call mpi_bcast (carma_diags_file,      1 ,mpi_integer, masterprocid, mpicom, ierr)
    call mpi_bcast (carma_maxsubsteps,     1 ,mpi_integer, masterprocid, mpicom, ierr)
    call mpi_bcast (carma_minsubsteps,     1 ,mpi_integer, masterprocid, mpicom, ierr)
    call mpi_bcast (carma_maxretries,      1 ,mpi_integer, masterprocid, mpicom, ierr)
    call mpi_bcast (carma_conmax,          1 ,mpi_real8, masterprocid, mpicom, ierr)
    call mpi_bcast (carma_dgc_threshold,   1 ,mpi_real8, masterprocid, mpicom, ierr)
    call mpi_bcast (carma_ds_threshold,    1 ,mpi_real8, masterprocid, mpicom, ierr)
    call mpi_bcast (carma_dt_threshold,    1 ,mpi_real8, masterprocid, mpicom, ierr)
    call mpi_bcast (carma_tstick,          1 ,mpi_real8, masterprocid, mpicom, ierr)
    call mpi_bcast (carma_gsticki,         1 ,mpi_real8, masterprocid, mpicom, ierr)
    call mpi_bcast (carma_gstickl,         1 ,mpi_real8, masterprocid, mpicom, ierr)
    call mpi_bcast (carma_cstick,          1 ,mpi_real8, masterprocid, mpicom, ierr)
    call mpi_bcast (carma_rhcrit,          1 ,mpi_real8, masterprocid, mpicom, ierr)
    call mpi_bcast (carma_vf_const,        1 ,mpi_real8, masterprocid, mpicom, ierr)
    call mpi_bcast (carma_model, len(carma_model), mpi_character, masterprocid, mpicom, ierr)
    call mpi_bcast (carma_sulfnuc_method, len(carma_sulfnuc_method), mpi_character, masterprocid, mpicom, ierr)
    call mpibcast  (carma_diags_packages, len(carma_diags_packages(1))*carma_maxdiags, mpi_character, 0, mpicom)
    call mpibcast  (carma_debug_packages, len(carma_debug_packages(1))*carma_maxdiags, mpi_character, 0, mpicom)

    carma_ndiagpkgs = 0
    do i = 1, carma_maxdiags
       if (len_trim(carma_diags_packages(i)) > 0) then
          carma_ndiagpkgs = carma_ndiagpkgs + 1
       endif
    enddo

    carma_ndebugpkgs = 0
    do i = 1, carma_maxdiags
       if (len_trim(carma_debug_packages(i)) > 0) then
          carma_ndebugpkgs = carma_ndebugpkgs + 1
       endif
    enddo

    ! Also cause the CARMA model flags to be read in.
    call carma_model_readnl(nlfile)

  end subroutine carma_readnl

end module carma_flags_mod
