module bam_optics_diags_mod
  use shr_kind_mod, only: r8 => shr_kind_r8
  use cam_history, only: fieldname_len, addfld, outfld, add_default, horiz_only
  use cam_history_support, only : fillvalue
  use rad_constituents, only: rad_cnst_get_info
  use ppgrid, only: pcols, pver
  use phys_control, only: phys_getopts

  implicit none

  private
  public :: bam_optics_diags_active
  public :: bam_optics_diags_init
  public :: bam_optics_diags_out

  integer :: numaerosols=0
  character(len=fieldname_len), pointer :: odv_names(:)  ! outfld names for visible OD
  logical :: bam_optics_diags_active = .false.

contains

  !==============================================================================
  subroutine bam_optics_diags_init()

    integer :: i
    character(len=64), allocatable :: aernames(:)
    logical :: history_aero_optics  ! Output aerosol optics diagnostics

    ! number of bulk aerosols in climate list
    call rad_cnst_get_info(0, naero=numaerosols)

    bam_optics_diags_active = numaerosols>0

    if (.not.bam_optics_diags_active) return

    ! get names of bulk aerosols
    allocate(aernames(numaerosols))
    call rad_cnst_get_info(0, aernames=aernames)

    call phys_getopts( history_aero_optics_out = history_aero_optics )

    ! diagnostic output for bulk aerosols
    ! create outfld names for visible OD
    allocate(odv_names(numaerosols))
    do i = 1, numaerosols
       odv_names(i) = 'ODV_'//trim(aernames(i))
       call addfld (odv_names(i), horiz_only, 'A', '1', &
            trim(aernames(i))//' optical depth in visible band', flag_xyfill=.true.)
       call add_default(odv_names(i), 1, ' ')
    end do

  end subroutine bam_optics_diags_init


  !==============================================================================
  subroutine bam_optics_diags_out(lchnk, ncol, nnite, idxnite, iaer, tau, diag_idx, troplev)

    ! output aerosol optical depth for the visible band

    integer,          intent(in) :: lchnk
    integer,          intent(in) :: ncol           ! number of columns
    integer,          intent(in) :: nnite          ! number of night columns
    integer,          intent(in) :: idxnite(:)     ! local column indices of night columns
    integer,          intent(in) :: iaer           ! aerosol index -- if 0 then tau is a total for all aerosols
    real(r8),         intent(in) :: tau(:,:)       ! aerosol optical depth for the visible band
    integer,          intent(in) :: diag_idx       ! identifies whether the aerosol optics
                                                   ! is for the climate calc or a diagnostic calc
    integer,          intent(in) :: troplev(:)     ! tropopause level

    ! Local variables
    integer  :: i
    real(r8) :: tmp(pcols), tmp2(pcols)
    !-----------------------------------------------------------------------------

!    if (.not.bam_optics_diags_active) return

    ! currently only implemented for climate calc
    if (diag_idx > 0) return

    ! compute total column aerosol optical depth
    tmp(1:ncol) = sum(tau(1:ncol,:), 2)
    ! use fillvalue to indicate night columns
    do i = 1, nnite
       tmp(idxnite(i)) = fillvalue
    end do

    if (iaer > 0) then
       call outfld(odv_names(iaer), tmp, pcols, lchnk)
    else
       call outfld('AEROD_v', tmp, pcols, lchnk)
       do i = 1, ncol
          tmp2(i) = sum(tau(i,:troplev(i)))
       end do
       call outfld('AODvstrt', tmp2, pcols, lchnk)
    end if

  end subroutine bam_optics_diags_out

  !==============================================================================

end module bam_optics_diags_mod
