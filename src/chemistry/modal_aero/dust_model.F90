!===============================================================================
! Dust for Modal Aerosol Model
!===============================================================================
module dust_model
  use shr_kind_mod,     only: r8 => shr_kind_r8, cl => shr_kind_cl
  use spmd_utils,       only: masterproc
  use cam_abortutils,   only: endrun
  use modal_aero_data,  only: ntot_amode, ndst=>nDust
  use cam_logfile,      only: iulog
  use shr_dust_emis_mod,only: is_dust_emis_zender, is_zender_soil_erod_from_atm

  implicit none
  private

  public :: dust_names
  public :: dust_nbin
  public :: dust_nnum
  public :: dust_indices
  public :: dust_emis
  public :: dust_readnl
  public :: dust_init
  public :: dust_active

  integer, protected :: dust_nbin != 2
  integer, protected :: dust_nnum != 2
  character(len=6), protected, allocatable :: dust_names(:)

  real(r8), allocatable :: dust_emis_sclfctr(:)

  integer , protected, allocatable :: dust_indices(:)
  real(r8), allocatable :: dust_dmt_vwr(:)

  real(r8)          :: dust_emis_fact = 0._r8     ! tuning parameter for dust emissions
  character(len=cl) :: soil_erod_file = 'none'    ! full pathname for soil erodibility dataset

  logical :: dust_active = .false.

 contains

  !=============================================================================
  ! reads dust namelist options
  !=============================================================================
  subroutine dust_readnl(nlfile)

    use namelist_utils,  only: find_group_name
    use spmd_utils,      only: mpicom, masterprocid, mpi_character, mpi_real8, mpi_success
    use shr_dust_emis_mod, only: shr_dust_emis_readnl

    character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

    ! Local variables
    integer :: unitn, ierr
    character(len=*), parameter :: subname = 'dust_readnl'

    namelist /dust_nl/ dust_emis_fact, soil_erod_file

    !-----------------------------------------------------------------------------

    ! Read namelist
    if (masterproc) then
       open( newunit=unitn, file=trim(nlfile), status='old' )
       call find_group_name(unitn, 'dust_nl', status=ierr)
       if (ierr == 0) then
          read(unitn, dust_nl, iostat=ierr)
          if (ierr /= 0) then
             call endrun(subname // ':: ERROR reading namelist')
          end if
       end if
       close(unitn)
    end if

    ! Broadcast namelist variables
    call mpi_bcast(soil_erod_file, len(soil_erod_file), mpi_character, masterprocid, mpicom, ierr)
    if (ierr/=mpi_success) then
       call endrun(subname//' MPI_BCAST ERROR: soil_erod_file')
    end if
    call mpi_bcast(dust_emis_fact, 1, mpi_real8, masterprocid, mpicom, ierr)
    if (ierr/=mpi_success) then
       call endrun(subname//' MPI_BCAST ERROR: dust_emis_fact')
    end if

    call shr_dust_emis_readnl(mpicom, 'drv_flds_in')

    if ((soil_erod_file /= 'none') .and. (.not.is_zender_soil_erod_from_atm())) then
       call endrun(subname//': should not specify soil_erod_file if Zender soil erosion is not in CAM')
    end if

    if (masterproc) then
       if (is_dust_emis_zender()) then
          write(iulog,*) subname,': Zender_2003 dust emission method is being used.'
       end if
       if (is_zender_soil_erod_from_atm()) then
          write(iulog,*) subname,': Zender soil erod file is handled in atm'
          write(iulog,*) subname,': soil_erod_file = ',trim(soil_erod_file)
          write(iulog,*) subname,': dust_emis_fact = ',dust_emis_fact
       end if
    end if

  end subroutine dust_readnl

  !=============================================================================
  !=============================================================================
  subroutine dust_init()
    use soil_erod_mod, only: soil_erod_init
    use constituents,  only: cnst_get_ind
    use aerosol_instances_mod, only: aerosol_instances_get_props, aerosol_instances_get_num_models
    use aerosol_properties_mod, only: aerosol_properties
    use physconst,     only: pi, rair, gravit
    use modal_dust_emissions, only: modal_dust_emissions_init

    integer :: l, m, mm, ndx, nspec, iaermod
    character(len=32) :: spec_name
    integer, parameter :: mymodes(7) = (/ 2, 1, 3, 4, 5, 6, 7 /) ! tricky order ...
    class(aerosol_properties), pointer :: aero_props_modal
    character(len=256) :: errmsg
    integer :: errflg

    dust_nbin = ndst
    dust_nnum = ndst

    allocate( dust_names(2*ndst) )
    allocate( dust_indices(2*ndst) )
    allocate( dust_emis_sclfctr(ndst) )
    allocate( dust_dmt_vwr(ndst) )

    ! Find modal properties object from factory
    aero_props_modal => null()
    do iaermod = 1, aerosol_instances_get_num_models()
       aero_props_modal => aerosol_instances_get_props(iaermod, 0)
       if (aero_props_modal%model_is('MAM')) exit
       aero_props_modal => null()
    end do

    ndx = 0
    do mm = 1, ntot_amode
       m = mymodes(mm)
       nspec = aero_props_modal%nspecies(m)
       do l = 1, nspec
          call aero_props_modal%get(m, l, specname=spec_name)
          if (spec_name(:3) == 'dst') then
             ndx=ndx+1
             dust_names(ndx) = spec_name
             dust_names(ndst+ndx) = 'num_'//spec_name(5:)
             call cnst_get_ind(dust_names(     ndx), dust_indices(     ndx))
             call cnst_get_ind(dust_names(ndst+ndx), dust_indices(ndst+ndx))
          endif
       enddo
    enddo

    dust_active = any(dust_indices(:) > 0)
    if (.not.dust_active) return

    if (is_zender_soil_erod_from_atm()) then
       call  soil_erod_init( dust_emis_fact, soil_erod_file )
    end if

    call modal_dust_emissions_init( ntot_amode=ntot_amode, dust_nbin=dust_nbin,   &
                                    pi=pi, rair=rair, gravit=gravit,              &
                                    dust_emis_sclfctr=dust_emis_sclfctr,          &
                                    dust_dmt_vwr=dust_dmt_vwr,                    &
                                    errmsg=errmsg, errflg=errflg )
    if (errflg /= 0) then
       call endrun('dust_init: '//trim(errmsg))
    end if

  end subroutine dust_init

  !===============================================================================
  !===============================================================================
  subroutine dust_emis( ncol, lchnk, dust_flux_in, cflx, soil_erod )
    use soil_erod_mod, only : soil_erodibility
    use physconst,     only : pi
    use modal_dust_emissions, only: modal_dust_emissions_run

  ! args
    integer,  intent(in)    :: ncol, lchnk
    real(r8), intent(in)    :: dust_flux_in(:,:)
    real(r8), intent(inout) :: cflx(:,:)
    real(r8), intent(out)   :: soil_erod(:)

  ! local vars
    logical  :: zender_soil_erod_from_atm
    real(r8) :: soil_erod_in(ncol)

    zender_soil_erod_from_atm = is_zender_soil_erod_from_atm()

    ! soil_erod_mod storage exists only when Zender soil erosion is in atm
    if (zender_soil_erod_from_atm) then
       soil_erod_in(:ncol) = soil_erodibility(:ncol,lchnk)
    else
       soil_erod_in(:ncol) = 0._r8
    end if

    call modal_dust_emissions_run( ncol=ncol, dust_nbin=dust_nbin,                    &
                                   dust_indices=dust_indices,                         &
                                   dust_emis_sclfctr=dust_emis_sclfctr,               &
                                   dust_dmt_vwr=dust_dmt_vwr,                         &
                                   dust_emis_fact=dust_emis_fact,                     &
                                   zender_soil_erod_from_atm=zender_soil_erod_from_atm, &
                                   soil_erodibility=soil_erod_in,                     &
                                   dust_flux_in=dust_flux_in,                         &
                                   pi=pi, cflx=cflx, soil_erod=soil_erod )

  end subroutine dust_emis

end module dust_model
