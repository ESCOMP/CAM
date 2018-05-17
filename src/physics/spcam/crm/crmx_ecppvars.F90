module crmx_ecppvars
#ifdef ECPP
  implicit none

  public 

  integer, public, parameter :: nupdraft_in  = 1  ! Number of updraft class
  integer, public, parameter :: ndndraft_in  = 1  ! Number of dndraft class
  integer, public, parameter :: ncls_ecpp_in = 3  ! Number of total number of ecpp transport class
                                                  ! = nupdraft_in+1+ndndraft_in
  integer, public, parameter :: ncc_in       = 2  ! number of clear/cloudy sub-calsses
  integer, public, parameter :: nprcp_in     = 2  ! Number of non-precipitating/precipitating sub-classes.

  integer, public, parameter :: QUI       = 1, &  !Quiescent class
                                UP1       = 2     !First index for upward classes

  integer, public  :: DN1, & !First index of downward classes
                      NCLASS_TR !Num. of transport classes
                                !Both initialized based on
                                !runtime settings

  integer, public :: NCLASS_CL = ncc_in, &  !Number of cloud classes
                     CLR = 1, &        !Clear sub-class
                     CLD = 2           !Cloudy sub-class

  integer, public :: NCLASS_PR = nprcp_in, &  !Number of precipitaion classes
                     PRN = 1,       &  !Not precipitating sub-class
                     PRY = 2           !Is precipitating sub-class


  real,dimension(:,:,:), allocatable :: qlsink, precr, precsolid, rh, qlsink_bf, prain, qcloud_bf, qvs

  real,dimension(:,:,:),allocatable :: &
       qcloudsum1, qcloud_bfsum1, qrainsum1, qicesum1, qsnowsum1, qgraupsum1, &
       qlsinksum1, qlsink_bfsum1, prainsum1, precrsum1, precsolidsum1, precallsum1, &
       altsum1, rhsum1, cf3dsum1, wwsum1, wwsqsum1, tkesgssum1, qvssum1 

! dim1 = z
  real,dimension(:),allocatable :: &
       xkhvsum, wup_thresh, wdown_thresh, wwqui_cen_sum, wwqui_bnd_sum, wwqui_cloudy_cen_sum, wwqui_cloudy_bnd_sum

! dims = (z, cloud sub-class, transport-class, precip sub-class)
  real, dimension(:,:,:,:), allocatable :: &
       area_bnd_final, area_bnd_sum, area_cen_final, area_cen_sum, &
       mass_bnd_final, mass_bnd_sum, mass_cen_final, mass_cen_sum, &
       ent_bnd_sum, rh_cen_sum, &
       qcloud_cen_sum, qcloud_bf_cen_sum, qrain_cen_sum, &
       qice_cen_sum, qsnow_cen_sum, qgraup_cen_sum, &
       qlsink_cen_sum, precr_cen_sum, precsolid_cen_sum, precall_cen_sum, &
       qlsink_bf_cen_sum, prain_cen_sum, qlsink_avg_cen_sum
#endif /*ECPP*/
end module crmx_ecppvars
