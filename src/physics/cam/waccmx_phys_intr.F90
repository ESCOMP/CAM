module waccmx_phys_intr
  use shr_kind_mod,   only: r8 => shr_kind_r8
  use physics_types,  only: physics_state, physics_ptend
  use physics_buffer, only: physics_buffer_desc

#ifdef WACCMX_PHYS
  use majorsp_diffusion, only: mspd_init
  use majorsp_diffusion, only: mspd_intr
  use ion_electron_temp, only: ion_electron_temp_readnl
  use ion_electron_temp, only: ion_electron_temp_init
  use ion_electron_temp, only: ion_electron_temp_register
  use ion_electron_temp, only: ion_electron_temp_inidat
  use ion_electron_temp, only: ion_electron_temp_tend
  use iondrag, only: iondrag_inidat
#endif

  implicit none
  private

  public :: waccmx_phys_mspd_init
  public :: waccmx_phys_mspd_tend
  public :: waccmx_phys_ion_elec_temp_reg
  public :: waccmx_phys_ion_elec_temp_inidat
  public :: waccmx_phys_ion_elec_temp_init
  public :: waccmx_phys_ion_elec_temp_tend
  public :: waccmx_phys_ion_elec_temp_readnl

contains

  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  subroutine waccmx_phys_mspd_init

#ifdef WACCMX_PHYS
    call mspd_init()
#endif
  end subroutine waccmx_phys_mspd_init

  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  subroutine waccmx_phys_mspd_tend(ztodt, state, ptend)

    real(r8), intent(in) :: ztodt                  ! 2 delta-t
    type(physics_state), intent(in)     :: state   ! Physics state variables
    type(physics_ptend), intent(inout)  :: ptend   ! indivdual parameterization tendencies

#ifdef WACCMX_PHYS
    call mspd_intr(ztodt, state, ptend)
#endif
  end subroutine waccmx_phys_mspd_tend

  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  subroutine waccmx_phys_ion_elec_temp_reg

#ifdef WACCMX_PHYS
    call ion_electron_temp_register
#endif

  end subroutine waccmx_phys_ion_elec_temp_reg

  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  subroutine waccmx_phys_ion_elec_temp_readnl(nlfilename)
    character(len=*), intent(in) :: nlfilename
#ifdef WACCMX_PHYS
    call ion_electron_temp_readnl(nlfilename)
#endif
  end subroutine waccmx_phys_ion_elec_temp_readnl

  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  subroutine waccmx_phys_ion_elec_temp_inidat(ncid_ini,pbuf2d)

    use pio, only : file_desc_t

    type(file_desc_t), intent(inout)   :: ncid_ini      ! Initial condition file id
    type(physics_buffer_desc), pointer :: pbuf2d(:,:)   ! Physics buffer

#ifdef WACCMX_PHYS
    call ion_electron_temp_inidat(ncid_ini,pbuf2d)
    call iondrag_inidat(ncid_ini,pbuf2d)
#endif
  end subroutine waccmx_phys_ion_elec_temp_inidat

  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  subroutine waccmx_phys_ion_elec_temp_init(pbuf2d)

    type(physics_buffer_desc), pointer :: pbuf2d(:,:)

#ifdef WACCMX_PHYS
    call ion_electron_temp_init(pbuf2d)
#endif
  end subroutine waccmx_phys_ion_elec_temp_init

  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  subroutine waccmx_phys_ion_elec_temp_tend(state, ptend, pbuf, ztodt)

    type(physics_state), intent(in)    :: state               ! physics state structure
    type(physics_ptend), intent(inout) :: ptend               ! parameterization tendency structure
    type(physics_buffer_desc),pointer  :: pbuf(:)             ! physics buffer

    real(r8),            intent(in)    :: ztodt               ! Physics time step

#ifdef WACCMX_PHYS
    call ion_electron_temp_tend(state, ptend, pbuf, ztodt)
#endif
  end subroutine waccmx_phys_ion_elec_temp_tend

end module waccmx_phys_intr
