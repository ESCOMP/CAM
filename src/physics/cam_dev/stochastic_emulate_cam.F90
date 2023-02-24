module stochastic_emulate_cam
! From Morrison (Lebo, originally TAU bin code)
! Gettelman and Chen 2018
!the subroutines take in air density, air temperature, and the bin mass boundaries, and 
!output the mass and number mixing ratio tendencies in each bin directly.
!this is then wrapped for CAM. 

! note, this is now coded locally. Want the CAM interface to be over i,k I think.

#ifndef HAVE_GAMMA_INTRINSICS
use shr_spfn_mod, only: gamma => shr_spfn_gamma
#endif

!use statements here
!use
  
use shr_kind_mod,      only: r8=>shr_kind_r8
use cam_history,       only: addfld
use micro_pumas_utils, only: pi, rhow, qsmall
use cam_logfile,       only: iulog

implicit none
private
save

! Subroutines
public :: stochastic_emulate_init_cam



!In the module top, declare the following so that these can be used throughout the module:

integer, parameter, public  :: ncd = 35
integer, parameter, public  :: ncdp = ncd + 1
integer, parameter, public  :: ncdl = ncd
integer, parameter, public  :: ncdpl = ncdl+1

!integer, private :: ncd,ncdp
!integer, private :: ncdl,ncdpl
!PARAMETER(ncd=35,ncdl=ncd) ! set number of ice and liquid bins
!PARAMETER(ncdp=ncd+1,ncdpl=ncdl+1)

! for Zach's collision-coalescence code


real(r8), private :: knn(ncd,ncd)

real(r8), public :: mmean(ncd), diammean(ncd)       ! kg & m at bin mid-points
real(r8), public :: medge(ncdp), diamedge(ncdp)     ! kg & m at bin edges 
integer, private  :: cutoff_id                       ! cutoff between cloud water and rain drop, D = 40 microns

!===============================================================================
contains
!===============================================================================

subroutine stochastic_emulate_init_cam

    use cam_history_support, only:          add_hist_coord

    call add_hist_coord('bins_ncd', ncd, 'bins for TAU microphysics')

    call addfld('amk_c',(/'lev     ','bins_ncd'/),'A','kg','cloud liquid mass from bins')
    call addfld('ank_c',(/'lev     ','bins_ncd'/),'A','1/kg','cloud liquid number concentration from bins')
    call addfld('amk_r',(/'lev     ','bins_ncd'/),'A','kg','cloud liquid mass from bins')
    call addfld('ank_r',(/'lev     ','bins_ncd'/),'A','1/kg','cloud liquid number concentration from bins')
    call addfld('amk',(/'lev     ','bins_ncd'/),'A','kg','all liquid mass from bins')
    call addfld('ank',(/'lev     ','bins_ncd'/),'A','1/kg','all liquid number concentration from bins')
    call addfld('amk_out',(/'lev     ','bins_ncd'/),'A','kg','all liquid mass from bins')
    call addfld('ank_out',(/'lev     ','bins_ncd'/),'A','1/kg','all liquid number concentration from bins')

    call addfld('scale_nc',(/'lev'/),'A','1','scaling factor for nc') 
    call addfld('scale_nr',(/'lev'/),'A','1','scaling factor for nr') 
    call addfld('scale_qc',(/'lev'/),'A','1','scaling factor for qc') 
    call addfld('scale_qr',(/'lev'/),'A','1','scaling factor for qr')
 
    call addfld('QC_TAU_in',(/'lev'/),'A','kg/kg','qc in TAU')
    call addfld('NC_TAU_in',(/'lev'/),'A','1/kg','nc in TAU')
    call addfld('QR_TAU_in',(/'lev'/),'A','kg/kg','qr in TAU')
    call addfld('NR_TAU_in',(/'lev'/),'A','1/kg','nr in TAU')
    call addfld('QC_TAU_out',(/'lev'/),'A','kg/kg','qc out TAU')
    call addfld('NC_TAU_out',(/'lev'/),'A','1/kg','nc out TAU')
    call addfld('QR_TAU_out',(/'lev'/),'A','kg/kg','qr out TAU')
    call addfld('NR_TAU_out',(/'lev'/),'A','1/kg','nr out TAU')  
   
    call addfld('qctend_MG2',(/'lev'/),'A','kg/kg/s','qc tendency due to autoconversion & accretion in MG2') 
    call addfld('nctend_MG2',(/'lev'/),'A','1/kg/s','nc tendency due to autoconversion & accretion in MG2')
    call addfld('qrtend_MG2',(/'lev'/),'A','kg/kg/s','qr tendency due to autoconversion & accretion in MG2') 
    call addfld('nrtend_MG2',(/'lev'/),'A','1/kg/s','nr tendency due to autoconversion & accretion in MG2')
    call addfld('qctend_TAU',(/'lev'/),'A','kg/kg/s','qc tendency due to autoconversion & accretion in MG2') 
    call addfld('nctend_TAU',(/'lev'/),'A','1/kg/s','nc tendency due to autoconversion & accretion in MG2')
    call addfld('qrtend_TAU',(/'lev'/),'A','kg/kg/s','qr tendency due to autoconversion & accretion in MG2') 
    call addfld('nrtend_TAU',(/'lev'/),'A','1/kg/s','nr tendency due to autoconversion & accretion in MG2')
    call addfld('qctend_TAU_diag',(/'lev'/),'A','kg/kg/s','qc tendency due to autoconversion & accretion in MG2')
    call addfld('nctend_TAU_diag',(/'lev'/),'A','1/kg/s','nc tendency due to autoconversion & accretion in MG2')
    call addfld('qrtend_TAU_diag',(/'lev'/),'A','kg/kg/s','qr tendency due to autoconversion & accretion in MG2')
    call addfld('nrtend_TAU_diag',(/'lev'/),'A','1/kg/s','nr tendency due to autoconversion & accretion in MG2')

    call addfld('gmnnn_lmnnn_TAU',(/'lev'/),'A','1','sum of mass gain and loss from bin code')
    call addfld('ML_fixer',(/'lev'/),'A','1','frequency of ML fixer is activated')
    call addfld('qc_fixer',(/'lev'/),'A','kg/kg','delta qc due to ML fixer')
    call addfld('nc_fixer',(/'lev'/),'A','kg/kg','delta nc due to ML fixer')
    call addfld('qr_fixer',(/'lev'/),'A','kg/kg','delta qr due to ML fixer')
    call addfld('nr_fixer',(/'lev'/),'A','kg/kg','delta nr due to ML fixer')

end subroutine stochastic_emulate_init_cam
end module stochastic_emulate_cam


