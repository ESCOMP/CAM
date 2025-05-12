!--------------------------------------------------------------------------
! CAM interface layer for inline computation of atmosphere ionization rates
! due to medium energy electrons in the magnetosphere radiation belts impacting
! the atmosphere.  Fluxes of electrons incident on the upper atmosphere can
! be computed based on Ap or read from file.
!--------------------------------------------------------------------------
module mee_ionization
  use shr_kind_mod, only: r8 => shr_kind_r8
  use solar_parms_data, only: Ap=>solar_parms_ap ! geomag activity index
  use mo_apex, only: alatm ! mag latitude at each column (radians)
  use ppgrid, only: pcols, pver
  use cam_logfile,  only: iulog
  use spmd_utils,   only: masterproc
  use cam_abortutils, only: endrun
  use mee_ap_util_mod,only: mee_ap_init, mee_ap_error, mee_ap_iprs

  implicit none

  private
  public :: mee_ion_readnl
  public :: mee_ion_init
  public :: mee_ion_final
  public :: mee_ionpairs

  logical :: mee_ion_inline = .false.
  logical :: mee_ion_diagonly = .false.
  real(r8) :: mee_ion_blc = -huge(1._r8) ! bounce cone angle (degrees)

contains

  !-----------------------------------------------------------------------------
  ! reads namelist options
  !-----------------------------------------------------------------------------
  subroutine mee_ion_readnl(nlfile)

    use namelist_utils, only: find_group_name
    use spmd_utils,     only: mpicom, mpi_logical, mpi_real8, masterprocid
    use mee_fluxes,     only: mee_fluxes_readnl

    character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

    ! Local variables
    integer :: unitn, ierr
    character(len=*), parameter :: subname = 'mee_ion_readnl'

    namelist /mee_ion_nl/ mee_ion_inline, mee_ion_blc, mee_ion_diagonly


    ! Read namelist
    if (masterproc) then
       open( newunit=unitn, file=trim(nlfile), status='old' )
       call find_group_name(unitn, 'mee_ion_nl', status=ierr)
       if (ierr == 0) then
          read(unitn, mee_ion_nl, iostat=ierr)
          if (ierr /= 0) then
             call endrun(subname // ':: ERROR reading namelist')
          end if
       end if
       close(unitn)
    end if

    ! Broadcast namelist variables
    call mpi_bcast(mee_ion_inline, 1, mpi_logical, masterprocid, mpicom, ierr)
    call mpi_bcast(mee_ion_blc, 1, mpi_real8, masterprocid, mpicom, ierr)
    call mpi_bcast(mee_ion_diagonly, 1, mpi_logical, masterprocid, mpicom, ierr)
    if ( masterproc ) then
       write(iulog,*) subname//':: mee_ion_inline = ', mee_ion_inline
       if (mee_ion_inline) then
          write(iulog,*) subname//':: mee_ion_blc = ', mee_ion_blc
          write(iulog,*) subname//':: mee_ion_diagonly = ', mee_ion_diagonly
       endif
    endif

    if (mee_ion_inline) then
       call mee_fluxes_readnl(nlfile)
    end if

  end subroutine mee_ion_readnl

  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  subroutine mee_ion_init()
    use cam_history, only: addfld
    use mee_fluxes,  only: mee_fluxes_init

    integer :: err

    if (.not.mee_ion_inline) return

    call mee_fluxes_init()
    call mee_ap_init(mee_ion_blc,err)
    if (err==mee_ap_error) then
       call endrun('mee_ion_init: not able to initialize Ap based MEE ionization')
    endif

    call addfld( 'APMEEionprs', (/ 'lev' /), 'A', 'pairs/cm3/sec', 'Ap generated MEE ionization rate' )
  end subroutine mee_ion_init

  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  subroutine mee_ion_final()
    use mee_fluxes,      only: mee_fluxes_final
    use mee_ap_util_mod, only: mee_ap_final

    if (.not.mee_ion_inline) return

    call mee_fluxes_final()
    call mee_ap_final()

  end subroutine mee_ion_final

  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  subroutine mee_ionpairs(ncol, lchnk, pmid, alt, temp, ionpairs)

    use air_composition, only: mbarv  ! kg/kmole
    use physconst,       only: gravit
    use air_composition, only: rairv  ! composition dependent gas constant (J/K/kg)
    use physconst,       only: boltz  ! Boltzman's constant (J/K/molecule)
    use physconst,       only: avogad ! Avogadro's number (molecules/kmole)
    use physconst,       only: rearth ! radius of earth (m)
    use cam_history,     only: outfld

    integer,  intent(in) :: ncol,lchnk
    real(r8), intent(in) :: pmid(:,:)
    real(r8), intent(in) :: alt(:,:) ! meters
    real(r8), intent(in) :: temp(:,:)
    real(r8), intent(out) :: ionpairs(:,:) ! ion pairs /cm3/sec

    real(r8) :: rho(pcols,pver)
    real(r8) :: scaleh(pcols,pver)
    real(r8) :: grvty(pcols,pver)
    integer :: err

    if (.not.mee_ion_inline) then
       ionpairs(:,:) = 0._r8
       return
    end if

    rho(:ncol,:) = pmid(:ncol,:)/(rairv(:ncol,:,lchnk)*temp(:ncol,:)) ! kg/m3
    rho(:ncol,:) = rho(:ncol,:)*1.0e-3_r8 ! kg/m3 --> g/cm3

    grvty(:ncol,:) = gravit * ( (rearth/(rearth+alt(:ncol,:)))**2 )

    scaleh(:ncol,:) = avogad * boltz*temp(:ncol,:)/(mbarv(:ncol,:,lchnk)*grvty(:ncol,:)) ! m
    scaleh(:ncol,:) = scaleh(:ncol,:) * 1.0e2_r8 ! m -> cm

    call mee_ap_iprs(ncol, pver, rho(:ncol,:), scaleh(:ncol,:), Ap, ionpairs(:ncol,:), &
                     status=err, maglat=alatm(:ncol,lchnk))
    if (err==mee_ap_error) then
       call endrun('mee_ionpairs: error in Ap based MEE ionization calculation')
    end if

    call outfld( 'APMEEionprs', ionpairs(:ncol,:), ncol, lchnk )

    if (mee_ion_diagonly) then
       ionpairs(:,:) = 0._r8
    end if

  end subroutine mee_ionpairs


end module mee_ionization
