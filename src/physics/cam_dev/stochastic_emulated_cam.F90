module stochastic_emulated_cam
! From Morrison (Lebo, originally TAU bin code)
! Gettelman and Chen 2018
!the subroutines take in air density, air temperature, and the bin mass boundaries, and 
!output the mass and number mixing ratio tendencies in each bin directly.
!this is then wrapped for CAM. 

use shr_kind_mod,      only: cl=>shr_kind_cl
use cam_history,       only: addfld
use cam_logfile,       only: iulog
use cam_abortutils,    only: endrun

implicit none
private
save

! Subroutines
public :: stochastic_emulated_readnl
public :: stochastic_emulated_init_cam

!Module variables
integer, parameter, public  :: ncd = 35
integer, parameter, public  :: ncdp = ncd + 1

character(len=cl) :: stochastic_emulated_filename_quantile = " "
character(len=cl) :: stochastic_emulated_filename_input_scale = " "
character(len=cl) :: stochastic_emulated_filename_output_scale = " "

!===============================================================================
contains
!===============================================================================

subroutine stochastic_emulated_readnl(nlfile)

  use namelist_utils,  only: find_group_name
  use units,           only: getunit, freeunit
  use spmd_utils,      only: mpicom, mstrid=>masterprocid, mpi_character, masterproc
  use string_utils,    only: int2str

  character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input
  integer :: unitn, ierr
  character(len=*), parameter :: sub = 'stochastic_emulated_readnl'

  namelist /stochastic_emulated_nl/ stochastic_emulated_filename_quantile, stochastic_emulated_filename_input_scale, &
                               stochastic_emulated_filename_output_scale

  if (masterproc) then
     unitn = getunit()
     open( unitn, file=trim(nlfile), status='old' )
     call find_group_name(unitn, 'stochastic_emulated_nl', status=ierr)
     if (ierr == 0) then
        read(unitn, stochastic_emulated_nl, iostat=ierr)
        if (ierr /= 0) then
           call endrun(sub // ':: ERROR reading namelist, iostat = ' // int2str(ierr))
        end if
     end if
     close(unitn)
     call freeunit(unitn)
  end if

  ! Broadcast namelist variables
  call mpi_bcast(stochastic_emulated_filename_quantile, cl, mpi_character, mstrid, mpicom, ierr)
  if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: stochastic_emulated_filename_quantile")

  call mpi_bcast(stochastic_emulated_filename_input_scale, cl, mpi_character, mstrid, mpicom, ierr)
  if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: stochastic_emulated_filename_input_scale")

  call mpi_bcast(stochastic_emulated_filename_output_scale, cl, mpi_character, mstrid, mpicom, ierr)
  if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: stochastic_emulated_filename_output_scale")

  write(iulog,*) 'PUMAS stochastic_emulated_readnl, stochastic_emulated_filename_quantile=',&
               stochastic_emulated_filename_quantile

end subroutine stochastic_emulated_readnl

subroutine stochastic_emulated_init_cam(stochastic_emulated_filename_quantile_out, &
                                       stochastic_emulated_filename_input_scale_out, &
                                       stochastic_emulated_filename_output_scale_out)

    use cam_history_support, only:          add_hist_coord

    character(len=cl),intent(out) :: stochastic_emulated_filename_quantile_out
    character(len=cl),intent(out) :: stochastic_emulated_filename_input_scale_out
    character(len=cl),intent(out) :: stochastic_emulated_filename_output_scale_out

    call add_hist_coord('bins_ncd', ncd, 'bins for TAU microphysics')

    call addfld('amk_c',(/'trop_cld_lev','bins_ncd    '/),'A','kg','cloud liquid mass from bins')
    call addfld('ank_c',(/'trop_cld_lev','bins_ncd    '/),'A','1/kg','cloud liquid number concentration from bins')
    call addfld('amk_r',(/'trop_cld_lev','bins_ncd    '/),'A','kg','rain mass from bins')
    call addfld('ank_r',(/'trop_cld_lev','bins_ncd    '/),'A','1/kg','rain number concentration from bins')
    call addfld('amk',(/'trop_cld_lev','bins_ncd    '/),'A','kg','all liquid mass from bins')
    call addfld('ank',(/'trop_cld_lev','bins_ncd    '/),'A','1/kg','all liquid number concentration from bins')
    call addfld('amk_out',(/'trop_cld_lev','bins_ncd    '/),'A','kg','all liquid mass from bins')
    call addfld('ank_out',(/'trop_cld_lev','bins_ncd    '/),'A','1/kg','all liquid number concentration from bins')

    call addfld('scale_nc',(/'trop_cld_lev'/),'A','1','scaling factor for nc') 
    call addfld('scale_nr',(/'trop_cld_lev'/),'A','1','scaling factor for nr') 
    call addfld('scale_qc',(/'trop_cld_lev'/),'A','1','scaling factor for qc') 
    call addfld('scale_qr',(/'trop_cld_lev'/),'A','1','scaling factor for qr')
 
    call addfld('QC_TAU_in',(/'trop_cld_lev'/),'A','kg/kg','qc in TAU')
    call addfld('NC_TAU_in',(/'trop_cld_lev'/),'A','1/kg','nc in TAU')
    call addfld('QR_TAU_in',(/'trop_cld_lev'/),'A','kg/kg','qr in TAU')
    call addfld('NR_TAU_in',(/'trop_cld_lev'/),'A','1/kg','nr in TAU')
    call addfld('QC_TAU_out',(/'trop_cld_lev'/),'A','kg/kg','qc out TAU')
    call addfld('NC_TAU_out',(/'trop_cld_lev'/),'A','1/kg','nc out TAU')
    call addfld('QR_TAU_out',(/'trop_cld_lev'/),'A','kg/kg','qr out TAU')
    call addfld('NR_TAU_out',(/'trop_cld_lev'/),'A','1/kg','nr out TAU')  
 
    call addfld('qctend_TAU',(/'trop_cld_lev'/),'A','kg/kg/s','qc tendency due to TAU bin code') 
    call addfld('nctend_TAU',(/'trop_cld_lev'/),'A','1/kg/s','nc tendency due to TAU bin code')
    call addfld('qrtend_TAU',(/'trop_cld_lev'/),'A','kg/kg/s','qr tendency due to TAU bin code') 
    call addfld('nrtend_TAU',(/'trop_cld_lev'/),'A','1/kg/s','nr tendency due to TAU bin code')
    call addfld('qctend_TAU_diag',(/'trop_cld_lev'/),'A','kg/kg/s','qc tendency due to TAU bin code')
    call addfld('nctend_TAU_diag',(/'trop_cld_lev'/),'A','1/kg/s','nc tendency due to TAU bin code')
    call addfld('qrtend_TAU_diag',(/'trop_cld_lev'/),'A','kg/kg/s','qr tendency due to TAU bin code')
    call addfld('nrtend_TAU_diag',(/'trop_cld_lev'/),'A','1/kg/s','nr tendency due to TAU bin code')

    call addfld('gmnnn_lmnnn_TAU',(/'trop_cld_lev'/),'A','1','sum of mass gain and loss from bin code')
    call addfld('ML_fixer',(/'trop_cld_lev'/),'A','1','frequency that ML fixer is activated')
    call addfld('qc_fixer',(/'trop_cld_lev'/),'A','kg/kg','delta qc due to ML fixer')
    call addfld('nc_fixer',(/'trop_cld_lev'/),'A','kg/kg','delta nc due to ML fixer')
    call addfld('qr_fixer',(/'trop_cld_lev'/),'A','kg/kg','delta qr due to ML fixer')
    call addfld('nr_fixer',(/'trop_cld_lev'/),'A','kg/kg','delta nr due to ML fixer')

    stochastic_emulated_filename_quantile_out     = stochastic_emulated_filename_quantile
    stochastic_emulated_filename_input_scale_out  = stochastic_emulated_filename_input_scale
    stochastic_emulated_filename_output_scale_out = stochastic_emulated_filename_output_scale

end subroutine stochastic_emulated_init_cam
end module stochastic_emulated_cam


