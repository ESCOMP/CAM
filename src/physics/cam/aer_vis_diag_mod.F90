!-------------------------------------------------------------------------------
! output aerosol optical depths for the visible band
!-------------------------------------------------------------------------------
module aer_vis_diag_mod
  use shr_kind_mod, only: r8 => shr_kind_r8
  use cam_history, only: fieldname_len, addfld, outfld, add_default, horiz_only, hist_fld_active
  use cam_history_support, only : fillvalue
  use aerosol_instances_mod, only: aerosol_instances_get_props, aerosol_instances_get_num_models
  use aerosol_properties_mod, only: aerosol_properties
  use ppgrid, only: pcols, pver
  use phys_control, only: phys_getopts
  use cam_abortutils, only: endrun

  implicit none

  private
  public :: aer_vis_diag_init
  public :: aer_vis_diag_out

  integer :: numaerosols=0
  character(len=fieldname_len), pointer :: odv_names(:)  ! outfld names for visible OD

contains

  !==============================================================================
  subroutine aer_vis_diag_init()

    integer :: i, iaermod, astat
    character(len=64), allocatable :: aernames(:)
    logical :: history_aero_optics  ! Output aerosol optics diagnostics
    class(aerosol_properties), pointer :: aero_props_bam

    ! Find BAM properties object from factory
    aero_props_bam => null()
    do iaermod = 1, aerosol_instances_get_num_models()
       aero_props_bam => aerosol_instances_get_props(iaermod, 0)
       if (associated(aero_props_bam)) then
          if (aero_props_bam%model_is('BAM')) exit
       end if
       aero_props_bam => null()
    end do

    ! number of bulk aerosols in climate list
    if (associated(aero_props_bam)) then
       numaerosols = aero_props_bam%nbins()
    else
       numaerosols = 0
    end if

    if (numaerosols<1) return

    ! get names of bulk aerosols
    allocate(aernames(numaerosols),stat=astat)
    if( astat/= 0 ) call endrun('aer_vis_diag_init: aernames allocate error')

    do i = 1, numaerosols
       aernames(i) = aero_props_bam%bin_name(i)
    end do

    call phys_getopts( history_aero_optics_out = history_aero_optics )

    ! diagnostic output for bulk aerosols
    ! create outfld names for visible OD
    allocate(odv_names(numaerosols),stat=astat)
    if( astat/= 0 ) call endrun('aer_vis_diag_init: odv_names allocate error')

    do i = 1, numaerosols
       odv_names(i) = 'ODV_'//trim(aernames(i))
       call addfld (odv_names(i), horiz_only, 'A', '1', &
            trim(aernames(i))//' optical depth in visible band', flag_xyfill=.true.)
       call add_default(odv_names(i), 1, ' ')
    end do

  end subroutine aer_vis_diag_init

  !==============================================================================
  subroutine aer_vis_diag_out(lchnk, ncol, nnite, idxnite, iaer, tau, diag_idx, troplev)

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
    logical :: do_calc
    !-----------------------------------------------------------------------------

    ! currently only implemented for climate calc
    if (diag_idx > 0) return

    do_calc = .false.
    if (iaer > 0) then
       do_calc = hist_fld_active(odv_names(iaer))
    else
       do_calc = hist_fld_active('AEROD_v')
    end if

    if (do_calc) then
       ! compute total column aerosol optical depth
       tmp(1:ncol) = sum(tau(1:ncol,:), 2)
       ! use fillvalue to indicate night columns
       do i = 1, nnite
          tmp(idxnite(i)) = fillvalue
       end do
    end if

    if (iaer > 0) then
       call outfld(odv_names(iaer), tmp, pcols, lchnk)
    else
       call outfld('AEROD_v', tmp, pcols, lchnk)

       if (hist_fld_active('AODvstrt')) then
          do i = 1, ncol
             tmp2(i) = sum(tau(i,:troplev(i)))
          end do
          call outfld('AODvstrt', tmp2, pcols, lchnk)
       end if
    end if

  end subroutine aer_vis_diag_out

  !==============================================================================

end module aer_vis_diag_mod
