!===============================================================================
! Seasalt for Modal Aerosol Model
!===============================================================================
module seasalt_model
  use shr_kind_mod,   only: r8 => shr_kind_r8
  use modal_aero_data,only: ntot_amode, nslt=>nSeaSalt

  implicit none
  private

  public :: seasalt_nbin
  public :: seasalt_nnum
  public :: seasalt_names
  public :: seasalt_indices
  public :: seasalt_init
  public :: seasalt_emis
  public :: seasalt_active

  integer, protected :: seasalt_nbin ! = nslt
  integer, protected :: seasalt_nnum ! = nnum

  character(len=6), protected, allocatable :: seasalt_names(:)
  integer, protected, allocatable :: seasalt_indices(:)

  logical :: seasalt_active = .false.

  real(r8):: emis_scale

contains
  
  !=============================================================================
  !=============================================================================
  subroutine seasalt_init(seasalt_emis_scale)
    use sslt_sections, only: sslt_sections_init
    use constituents,  only: cnst_get_ind
    use aerosol_instances_mod, only: aerosol_instances_get_props, aerosol_instances_get_num_models
    use aerosol_properties_mod, only: aerosol_properties

    real(r8), intent(in) :: seasalt_emis_scale
    integer :: m, l, nspec, ndx, iaermod
    character(len=32) :: spec_name
    class(aerosol_properties), pointer :: aero_props_modal
    
    seasalt_nbin = nslt
    seasalt_nnum = nslt
    allocate(seasalt_names(2*nslt))
    allocate(seasalt_indices(2*nslt))

    ! Find modal properties object from factory
    aero_props_modal => null()
    do iaermod = 1, aerosol_instances_get_num_models()
       aero_props_modal => aerosol_instances_get_props(iaermod, 0)
       if (aero_props_modal%model_is('MAM')) exit
       aero_props_modal => null()
    end do

    ndx=0
    do m = 1, ntot_amode
       nspec = aero_props_modal%nspecies(m)
       do l = 1, nspec
          call aero_props_modal%get(m, l, specname=spec_name)
          if (spec_name(:3) == 'ncl') then
             ndx=ndx+1
             seasalt_names(ndx) = spec_name
             seasalt_names(nslt+ndx) = 'num_'//spec_name(5:)
             call cnst_get_ind(seasalt_names(     ndx), seasalt_indices(     ndx))
             call cnst_get_ind(seasalt_names(nslt+ndx), seasalt_indices(nslt+ndx))
          endif
       enddo
    enddo

    seasalt_active = any(seasalt_indices(:) > 0)
    if (.not.seasalt_active) return

    call sslt_sections_init()

    emis_scale = seasalt_emis_scale

  end subroutine seasalt_init

  !=============================================================================
  !=============================================================================
  subroutine seasalt_emis( u_bottom, v_bottom, zmid_bottom, srf_temp, ocnfrc, ncol, cflx )

    use physconst, only: pi
    use modal_seasalt_emissions, only: modal_seasalt_emissions_run

    ! dummy arguments
    real(r8), intent(in) :: u_bottom(:)     ! bottom layer zonal wind (m/s)
    real(r8), intent(in) :: v_bottom(:)     ! bottom layer meridional wind (m/s)
    real(r8), intent(in) :: zmid_bottom(:)  ! bottom layer midpoint geopotential height above surface (m)
    real(r8), intent(in) :: srf_temp(:)
    real(r8), intent(in) :: ocnfrc(:)
    integer,  intent(in) :: ncol
    real(r8), intent(inout) :: cflx(:,:)

    call modal_seasalt_emissions_run( ncol=ncol, nslt=nslt,                      &
                                      seasalt_indices=seasalt_indices,           &
                                      emis_scale=emis_scale,                     &
                                      u_bottom=u_bottom, v_bottom=v_bottom,      &
                                      zmid_bottom=zmid_bottom,                   &
                                      srf_temp=srf_temp, ocnfrc=ocnfrc,          &
                                      pi=pi, cflx=cflx )

  end subroutine seasalt_emis

end module seasalt_model
