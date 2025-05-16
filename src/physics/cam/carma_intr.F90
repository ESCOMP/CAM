!! This module is stub for a coupler between the CAM model and the Community Aerosol
!! and Radiation Model for Atmospheres (CARMA) microphysics model. It is used when
!! CARMA is not being used, so that the CAM code that calls CARMA does not need to
!! be changed. The real version of this routine exists in the directory
!! physics/carma/cam. A CARMA model can be activated by using configure with the
!! option:
!!
!!  -carma <carma_pkg>
!!
!! where carma_pkg is the name for a particular microphysical model.
!!
!! @author  Chuck Bardeen
!! @version May 2009
module carma_intr

  use shr_kind_mod,   only: r8 => shr_kind_r8
  use pmgrid,         only: plat, plev, plevp, plon
  use ppgrid,         only: pcols, pver, pverp
  use constituents,   only: pcnst
  use physics_types,  only: physics_state, physics_ptend, physics_ptend_init
  use physics_buffer, only: physics_buffer_desc


  implicit none

  private
  save

  ! Public interfaces

  ! CAM Physics Interface
  public carma_register                 ! register consituents
  public carma_is_active                ! retrns true if this package is active (microphysics = .true.)
  public carma_implements_cnst          ! returns true if consituent is implemented by this package
  public carma_init_cnst                ! initialize constituent mixing ratios, if not read from initial file
  public carma_init                     ! initialize timestep independent variables
  public carma_final                    ! finalize the CARMA module
  public carma_timestep_init            ! initialize timestep dependent variables
  public carma_timestep_tend            ! interface to tendency computation
  public carma_accumulate_stats         ! collect stats from all MPI tasks

  ! Other Microphysics
  public carma_emission_tend            ! calculate tendency from emission source function
  public carma_calculate_cloudborne_diagnostics ! calculate model specific budget diagnostics for cloudborne aerosols
  public carma_output_cloudborne_diagnostics ! output model specific budget diagnostics for cloudborne aerosols
  public carma_output_budget_diagnostics ! calculate and output model specific aerosol budget terms
  public carma_wetdep_tend              ! calculate tendency from wet deposition

  public :: carma_restart_init
  public :: carma_restart_write
  public :: carma_restart_read

  public carma_get_bin
  public carma_get_bin_cld
  public carma_get_dry_radius
  public carma_get_elem_for_group
  public carma_get_group_by_name
  public carma_get_kappa
  public carma_get_number
  public carma_get_number_cld
  public carma_get_total_mmr
  public carma_get_total_mmr_cld
  public carma_get_wet_radius
  public carma_get_bin_rmass
  public carma_set_bin
  public carma_get_sad
  public :: carma_get_wght_pct
  public :: carma_effecitive_radius

  public :: carma_get_bin_radius

  integer, parameter, public  ::     MAXCLDAERDIAG = 16

contains


  subroutine carma_register
    implicit none

    return
  end subroutine carma_register


  function carma_is_active()
    implicit none

    logical :: carma_is_active

    carma_is_active = .false.

    return
  end function carma_is_active


  function carma_implements_cnst(name)
    implicit none

    character(len=*), intent(in) :: name   !! constituent name
    logical :: carma_implements_cnst       ! return value

    carma_implements_cnst = .false.

    return
  end function carma_implements_cnst


  subroutine carma_init(pbuf2d)
    implicit none
    type(physics_buffer_desc), pointer :: pbuf2d(:,:)

    return
  end subroutine carma_init


  subroutine carma_final
    implicit none

    return
  end subroutine carma_final


  subroutine carma_timestep_init
    implicit none

    return
  end subroutine carma_timestep_init


  subroutine carma_timestep_tend(state, cam_in, cam_out, ptend, dt, pbuf, dlf, rliq, prec_str, snow_str, &
    prec_sed, snow_sed, ustar, obklen)
    use hycoef,           only: hyai, hybi, hyam, hybm
    use time_manager,     only: get_nstep, get_step_size, is_first_step
    use camsrfexch,       only: cam_in_t, cam_out_t
    use scamMod,          only: single_column

    implicit none

    type(physics_state), intent(inout) :: state                 !! physics state variables
    type(cam_in_t), intent(in)         :: cam_in                !! surface inputs
    type(cam_out_t), intent(inout)     :: cam_out               !! cam output to surface models
    type(physics_ptend), intent(out)   :: ptend                 !! constituent tendencies
    real(r8), intent(in)               :: dt                    !! time step (s)
    type(physics_buffer_desc), pointer :: pbuf(:)               !! physics buffer
    real(r8), intent(in), optional     :: dlf(pcols,pver)       !! Detraining cld H20 from convection (kg/kg/s)
    real(r8), intent(inout), optional  :: rliq(pcols)           !! vertical integral of liquid not yet in q(ixcldliq)
    real(r8), intent(out), optional    :: prec_str(pcols)       !! [Total] sfc flux of precip from stratiform (m/s)
    real(r8), intent(out), optional    :: snow_str(pcols)       !! [Total] sfc flux of snow from stratiform   (m/s)
    real(r8), intent(out), optional    :: prec_sed(pcols)       !! total precip from cloud sedimentation (m/s)
    real(r8), intent(out), optional    :: snow_sed(pcols)       !! snow from cloud ice sedimentation (m/s)
    real(r8), intent(in), optional     :: ustar(pcols)          !! friction velocity (m/s)
    real(r8), intent(in), optional     :: obklen(pcols)         !! Obukhov length [ m ]

    call physics_ptend_init(ptend,state%psetcols,'none') !Initialize an empty ptend for use with physics_update

    if (present(prec_str))  prec_str(:)    = 0._r8
    if (present(snow_str))  snow_str(:)    = 0._r8
    if (present(prec_sed))  prec_sed(:)    = 0._r8
    if (present(snow_sed))  snow_sed(:)    = 0._r8

    return
  end subroutine carma_timestep_tend


  subroutine carma_init_cnst(name, latvals, lonvals, mask, q)
    implicit none

    character(len=*), intent(in)  :: name       !! constituent name
    real(r8),         intent(in)  :: latvals(:) !! lat in degrees (ncol)
    real(r8),         intent(in)  :: lonvals(:) !! lon in degrees (ncol)
    logical,          intent(in)  :: mask(:)    !! Only initialize where .true.
    real(r8),         intent(out) :: q(:,:)     !! mass mixing ratio

    if (name == "carma") then
      q = 0._r8
    end if

    return
  end subroutine carma_init_cnst

  subroutine carma_calculate_cloudborne_diagnostics(state, pbuf, aerclddiag)

    implicit none

    type(physics_state), intent(in)    :: state          !! Physics state variables - before CARMA
    type(physics_buffer_desc), pointer :: pbuf(:)        !! physics buffer
    real(r8), intent(out)              :: aerclddiag(pcols, MAXCLDAERDIAG)  !! previous cloudborne diagnostics

    return
  end subroutine carma_calculate_cloudborne_diagnostics


  subroutine carma_output_cloudborne_diagnostics(state, pbuf, pname, dt, oldaerclddiag)

    implicit none

    type(physics_state), intent(in)    :: state          !! Physics state variables - before CARMA
    type(physics_buffer_desc), pointer :: pbuf(:)        !! physics buffer
    character(*), intent(in)           :: pname          !! short name of the physics package
    real(r8), intent(in)               :: dt             !! timestep (s)
    real(r8), intent(in)               :: oldaerclddiag(pcols, MAXCLDAERDIAG)  !! previous cloudborne diagnostics

    return
  end subroutine carma_output_cloudborne_diagnostics


  subroutine carma_output_budget_diagnostics(state, ptend, old_cflux, cflux, dt, pname)

    implicit none

    type(physics_state), intent(in)    :: state          !! Physics state variables - before CARMA
    type(physics_ptend), intent(in)    :: ptend          !! indivdual parameterization tendencies
    real(r8)                           :: old_cflux(pcols,pcnst)  !! cam_in%clfux from before the timestep_tend
    real(r8)                           :: cflux(pcols,pcnst)  !! cam_in%clfux from after the timestep_tend
    real(r8), intent(in)               :: dt             !! timestep (s)
    character(*), intent(in)           :: pname          !! short name of the physics package


    return
  end subroutine carma_output_budget_diagnostics

  subroutine carma_emission_tend(state, ptend, cam_in, dt, pbuf)
    use camsrfexch,       only: cam_in_t

    implicit none

    type(physics_state), intent(in )    :: state                !! physics state
    type(physics_ptend), intent(inout)  :: ptend                !! physics state tendencies
    type(cam_in_t),      intent(inout)  :: cam_in               !! surface inputs
    real(r8),            intent(in)     :: dt                   !! time step (s)
    type(physics_buffer_desc), pointer  :: pbuf(:)              !! physics buffer

    return
  end subroutine carma_emission_tend


  subroutine carma_wetdep_tend(state, ptend, dt,  pbuf, dlf, cam_out)
    use camsrfexch,       only: cam_out_t

    implicit none

    real(r8),             intent(in)    :: dt             !! time step (s)
    type(physics_state),  intent(in )   :: state          !! physics state
    type(physics_ptend),  intent(inout) :: ptend          !! physics state tendencies
    type(physics_buffer_desc), pointer  :: pbuf(:)        !! physics buffer
    real(r8), intent(in)                :: dlf(pcols,pver)       !! Detraining cld H20 from convection (kg/kg/s)
    type(cam_out_t),      intent(inout) :: cam_out        !! cam output to surface models

    return
  end subroutine carma_wetdep_tend


  subroutine carma_accumulate_stats()
    implicit none

  end subroutine carma_accumulate_stats

  !---------------------------------------------------------------------------
  ! define fields for reference profiles in cam restart file
  !---------------------------------------------------------------------------
  subroutine CARMA_restart_init( File )
    use pio, only: file_desc_t

    ! arguments
    type(file_desc_t),intent(inout) :: File     ! pio File pointer

  end subroutine CARMA_restart_init

  !---------------------------------------------------------------------------
  ! write reference profiles to restart file
  !---------------------------------------------------------------------------
  subroutine CARMA_restart_write(File)
    use pio, only: file_desc_t

    ! arguments
    type(file_desc_t), intent(inout) :: File

  end subroutine CARMA_restart_write

  !---------------------------------------------------------------------------
  ! read reference profiles from restart file
  !---------------------------------------------------------------------------
  subroutine CARMA_restart_read(File)
    use pio, only: file_desc_t

    ! arguments
    type(file_desc_t),intent(inout) :: File     ! pio File pointer

  end subroutine CARMA_restart_read


  !! Get the mixing ratio for the specified element and bin.
  subroutine carma_get_bin(state, ielem, ibin, mmr, rc)
    type(physics_state), intent(in)   :: state                 !! physics state variables
    integer, intent(in)               :: ielem                 !! element index
    integer, intent(in)               :: ibin                  !! bin index
    real(r8), intent(out)             :: mmr(pcols,pver)       !! mass mixing ratio (kg/kg)
    integer, intent(out)              :: rc                    !! return code

  end subroutine carma_get_bin
  !! Get the mixing ratio for the specified element and bin.
  subroutine carma_get_bin_cld(pbuf, ielem, ibin, ncol, nlev, mmr, rc)
    type(physics_buffer_desc), pointer :: pbuf(:)               !! physics buffer
    integer, intent(in)               :: ielem                 !! element index
    integer, intent(in)               :: ibin                  !! bin index
    integer, intent(in)               :: ncol,nlev                  !! dimensions
    real(r8), intent(out)             :: mmr(:,:)       !! mass mixing ratio (kg/kg)
    integer, intent(out)              :: rc                    !! return code

  end subroutine carma_get_bin_cld
  !! Determine the dry radius and dry density for the particular bin.
  subroutine carma_get_dry_radius(state, igroup, ibin, rdry, rhopdry, rc)
    type(physics_state), intent(in)   :: state                 !! physics state variables
    integer, intent(in)               :: igroup                !! group index
    integer, intent(in)               :: ibin                  !! bin index
    real(r8), intent(out)             :: rdry(:,:)             !! dry radius (m)
    real(r8), intent(out)             :: rhopdry(:,:)          !! dry density (kg/m3)
    integer, intent(out)              :: rc                    !! return code

  end subroutine carma_get_dry_radius
  !! Get the number of elements and list of element ids for a group. This includes
  subroutine carma_get_elem_for_group(igroup, nelems, ielems, rc)
    integer, intent(in)               :: igroup                !! group index
    integer, intent(out)              :: nelems                !! number of elements in group
    integer, intent(out)              :: ielems(:)             !! indexes of elements in group
    integer, intent(out)              :: rc                    !! return code
  end subroutine carma_get_elem_for_group
  !! Get the CARMA group id a group name.
  subroutine carma_get_group_by_name(shortname, igroup, rc)
    character(len=*), intent(in)       :: shortname             !! the group short name
    integer, intent(out)               :: igroup                !! group index
    integer, intent(out)               :: rc                    !! return code

  end subroutine carma_get_group_by_name
  !! Get the CARMA group id and bin id from a compound name xxxxxxnn, where xxxxxx is the
  subroutine carma_get_group_and_bin_by_name(shortname, igroup, ibin, rc)
    character(len=*), intent(out)      :: shortname             !! the group short name
    integer, intent(out)               :: igroup                !! group index
    integer, intent(out)               :: ibin                  !! bin index
    integer, intent(out)               :: rc                    !! return code

  end subroutine carma_get_group_and_bin_by_name
  !! Determine a mass weighted kappa for the entire particle.
  subroutine carma_get_kappa(state, igroup, ibin, kappa, rc)
    type(physics_state), intent(in)   :: state                 !! physics state variables
    integer, intent(in)               :: igroup                !! group index
    integer, intent(in)               :: ibin                  !! bin index
    real(r8), intent(out)             :: kappa(:,:)     !! kappa value for the entire particle
    integer, intent(out)              :: rc                    !! return code
  end subroutine carma_get_kappa
  !! Get the number mixing ratio for the group. This is the number of particles per
  subroutine carma_get_number(state, igroup, ibin, nmr, rc)
    type(physics_state), intent(in)   :: state                 !! physics state variables
    integer, intent(in)               :: igroup                !! group index
    integer, intent(in)               :: ibin                  !! bin index
    real(r8), intent(out)             :: nmr(pcols,pver)       !! number mixing ratio (#/kg)
    integer, intent(out)              :: rc                    !! return code
  end subroutine carma_get_number

  subroutine carma_get_number_cld(pbuf, igroup, ibin, ncol, nlev, nmr, rc)
    type(physics_buffer_desc),pointer :: pbuf(:)               !! physics buffer
    integer, intent(in)               :: igroup                !! group index
    integer, intent(in)               :: ibin                  !! bin index
    integer, intent(in)               :: ncol,nlev             !! dimensions
    real(r8), intent(out)             :: nmr(pcols,pver)       !! number mixing ratio (#/kg)
    integer, intent(out)              :: rc                    !! return code
  end subroutine carma_get_number_cld
  !! Get the mixing ratio for the group. This is the total of all the elements that
  subroutine carma_get_total_mmr(state, igroup, ibin, totmmr, rc)
    type(physics_state), intent(in)   :: state                 !! physics state variables
    integer, intent(in)               :: igroup                !! group index
    integer, intent(in)               :: ibin                  !! bin index
    real(r8), intent(out)             :: totmmr(pcols,pver)    !! total mmr (kg/kg)
    integer, intent(out)              :: rc                    !! return code
  end subroutine carma_get_total_mmr

  subroutine carma_get_total_mmr_cld(pbuf, igroup, ibin, ncol, nlev, totmmr, rc)
    type(physics_buffer_desc),pointer :: pbuf(:)               !! physics buffer
    integer, intent(in)               :: igroup                !! group index
    integer, intent(in)               :: ibin                  !! bin index
    integer, intent(in)               :: ncol,nlev             !! dimensions
    real(r8), intent(out)             :: totmmr(pcols,pver)    !! total mmr (kg/kg)
    integer, intent(out)              :: rc                    !! return code

  end subroutine carma_get_total_mmr_cld

  subroutine carma_get_sad(state, igroup, ibin, sad, rc)
    type(physics_state), intent(in)   :: state                 !! physics state variables
    integer, intent(in)               :: igroup                !! group index
    integer, intent(in)               :: ibin                  !! bin index
    real(r8), intent(out)             :: sad(pcols,pver)       !! surface area dens (cm2/cm3)
    integer, intent(out)              :: rc                    !! return code
  end subroutine carma_get_sad

  !! Find the wet radius and wet density for the group and bin specified.
  subroutine carma_get_wet_radius(state, igroup, ibin, rwet, rhopwet, rc)
    type(physics_state), intent(in)   :: state                 !! physics state variables
    integer, intent(in)               :: igroup                !! group index
    integer, intent(in)               :: ibin                  !! bin index
    real(r8), intent(out)             :: rwet(pcols,pver)      !! wet radius (m)
    real(r8), intent(out)             :: rhopwet(pcols,pver)   !! wet density (kg/m3)
    integer, intent(inout)            :: rc                    !! return code

  end subroutine carma_get_wet_radius
  !! Provides the tendency (in kg/kg/s) required to change the element and bin from
  !! the current state to the desired mmr.
  subroutine carma_set_bin(state, ielem, ibin, mmr, dt, ptend, rc)
    type(physics_state), intent(in)    :: state                 !! physics state variables
    integer, intent(in)                :: ielem                 !! element index
    integer, intent(in)                :: ibin                  !! bin index
    real(r8), intent(in)               :: mmr(pcols,pver)       !! mass mixing ratio (kg/kg)
    integer                            :: dt                    !! timestep size (sec)
    type(physics_ptend), intent(inout) :: ptend                 !! constituent tendencies
    integer, intent(out)               :: rc                    !! return code
  end subroutine carma_set_bin

  subroutine carma_get_bin_rmass(igroup, ibin, mass, rc)

    integer, intent(in)               :: igroup                !! group index
    integer, intent(in)               :: ibin                  !! bin index
    real(r8),intent(out)              :: mass ! grams ???
    integer, intent(out)              :: rc                    !! return code

  end subroutine carma_get_bin_rmass

  function carma_get_wght_pct(icol,ilev,state) result(wtpct)

    integer, intent(in) ::  icol,ilev
    type(physics_state), intent(in)    :: state          !! Physics state variables - before CARMA

    real(r8) :: wtpct

  end function carma_get_wght_pct


  function carma_effecitive_radius(state) result(rad)

    type(physics_state), intent(in)   :: state                 !! physics state variables
    real(r8) :: rad(pcols,pver) ! effective radius (cm)
  end function carma_effecitive_radius

  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  subroutine carma_get_bin_radius(igroup, ibin, radius, rc)
    integer, intent(in)               :: igroup                !! group index
    integer, intent(in)               :: ibin                  !! bin index
    real(r8),intent(out)              :: radius ! cm ???
    integer, intent(out)              :: rc                    !! return code
  end subroutine carma_get_bin_radius

end module carma_intr
