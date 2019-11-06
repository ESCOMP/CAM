!================================================================================================
! This is the "GEOS-Chem" chemistry module.
!================================================================================================

module chemistry
  use shr_kind_mod,        only: r8 => shr_kind_r8
  use physics_types,       only: physics_state, physics_ptend, physics_ptend_init
  use ppgrid,              only: begchunk, endchunk, pcols
  use ppgrid,              only: pver
  use constituents,        only: pcnst, cnst_add
  use mo_gas_phase_chemdr, only: map2chm
  use mo_constants,        only: pi
  use shr_const_mod,       only: molw_dryair=>SHR_CONST_MWDAIR
  use mo_chem_utls,        only : get_spc_ndx
  use chem_mods,           only : gas_pcnst, adv_mass
  use mo_sim_dat, only: set_sim_dat
  implicit none
  private
  save
  !
  ! Public interfaces
  !
  public :: chem_is                        ! identify which chemistry is being used
  public :: chem_register                  ! register consituents
  public :: chem_is_active                 ! returns true if this package is active (ghg_chem=.true.)
  public :: chem_implements_cnst           ! returns true if consituent is implemented by this package
  public :: chem_init_cnst                 ! initialize mixing ratios if not read from initial file
  public :: chem_init                      ! initialize (history) variables
  public :: chem_timestep_tend             ! interface to tendency computation
  public :: chem_final
  public :: chem_write_restart
  public :: chem_read_restart
  public :: chem_init_restart
  public :: chem_readnl                    ! read chem namelist 
  public :: chem_reset_fluxes

  interface chem_write_restart
     module procedure chem_write_restart_bin
     module procedure chem_write_restart_pio
  end interface
  interface chem_read_restart
     module procedure chem_read_restart_bin
     module procedure chem_read_restart_pio
  end interface

  ! Private data

!================================================================================================
contains
!================================================================================================

  logical function chem_is (name)

    character(len=*), intent(in) :: name

    chem_is = .false.
    if (name == 'geoschem' ) then
       chem_is = .true.
    end if

  end function chem_is

!================================================================================================

  subroutine chem_register
    !----------------------------------------------------------------------- 
    ! 
    ! Purpose: register advected constituents for parameterized greenhouse gas chemistry
    ! 
    !-----------------------------------------------------------------------

  end subroutine chem_register

  subroutine chem_readnl(nlfile)
    ! Gonna want to read in input.geos, etc. Somehow.

    ! args
    character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

  end subroutine chem_readnl

!================================================================================================

  function chem_is_active()
    !-----------------------------------------------------------------------
    logical :: chem_is_active
    !-----------------------------------------------------------------------
    chem_is_active = .true.
  end function chem_is_active

!================================================================================================

  function chem_implements_cnst(name)
    !----------------------------------------------------------------------- 
    ! 
    ! Purpose: return true if specified constituent is implemented by this package
    ! 
    ! Author: B. Eaton
    ! 
    !-----------------------------------------------------------------------
    implicit none
    !-----------------------------Arguments---------------------------------

    character(len=*), intent(in) :: name   ! constituent name
    logical :: chem_implements_cnst        ! return value

    integer :: i
    
    chem_implements_cnst = .false.

    do i = 1, nspecies
       if (trim(species(i)) .eq. trim(name)) then
          chem_implements_cnst = .true.
          exit
       end if
    end do

  end function chem_implements_cnst

!===============================================================================
  
  subroutine chem_init(phys_state, pbuf2d)
    !----------------------------------------------------------------------- 
    ! 
    ! Purpose: initialize GEOS-Chem parts (state objects, mainly)
    !          (and declare history variables)
    ! 
    !-----------------------------------------------------------------------

    ! Can add history output here too with the "addfld" & "add_default" routines

  end subroutine chem_init

!===============================================================================

  subroutine chem_timestep_tend( state, ptend, cam_in, cam_out, dt, pbuf,  fh2o )

    ! Here's where you'll call DO_CHEMISTRY

    return
  end subroutine chem_timestep_tend

!===============================================================================
  subroutine chem_init_cnst(name, latvals, lonvals, mask, q)

    ! Used to initialize tracer fields if desired.
    ! Will need a simple mapping structure as well as the CAM tracer registration
    ! routines.

  end subroutine chem_init_cnst

!===============================================================================
  subroutine chem_final
    
    ! Finalize GEOS-Chem

    return
  end subroutine chem_final
!===============================================================================
  subroutine chem_write_restart_bin( nrg )
    implicit none
    integer,intent(in) :: nrg     ! Unit number
    return
  end subroutine chem_write_restart_bin
!===============================================================================
  subroutine chem_read_restart_bin( nrg )
    implicit none
    integer,intent(in) :: nrg     ! Unit number
    return
  end subroutine chem_read_restart_bin
!===============================================================================
  subroutine chem_write_restart_pio( File )
    use pio, only : file_desc_t
    type(file_desc_t) :: File
    return
  end subroutine chem_write_restart_pio
!===============================================================================
  subroutine chem_read_restart_pio( File )
    use pio, only : file_desc_t
    type(file_desc_t) :: File
    return
  end subroutine chem_read_restart_pio
!===============================================================================
  subroutine chem_init_restart(File)
    use pio, only : file_desc_t
    type(file_desc_t) :: File
    return
  end subroutine chem_init_restart
!================================================================================
  subroutine chem_reset_fluxes( fptr, cam_in )
    use camsrfexch, only : cam_in_t     

    real(r8), pointer             :: fptr(:,:)        ! pointer into    array data
    type(cam_in_t), intent(inout) :: cam_in(begchunk:endchunk)

  end subroutine chem_reset_fluxes

end module chemistry
