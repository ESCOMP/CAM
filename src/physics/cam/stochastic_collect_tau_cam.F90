module stochastic_collect_tau_cam
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
public :: stochastic_kernel_init, stochastic_collect_tau_tend



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

      
subroutine calc_bins    

real(r8) :: DIAM(ncdp)
real(r8) :: X(ncdp)
real(r8) :: radsl(ncdp)
real(r8) :: radl(ncd)
integer  :: L, lcl  
real(r8) :: kkfac
!Then before doing any calculations you'll need to calculate the bin mass grid 
! (note this code could be cleaned up, I'm just taking it as it's used in our bin scheme). 
! This only needs to be done once, since we'll use the same bin mass grid for all calculations. 

! use mass doubling bins from Graham Feingold (note cgs units)

!      PI=3.14159_r8
!      rho_w=1000._r8                         ! kg/m3
      DIAM(1)=1.5625*2.E-04_r8                ! cm
      X(1)=PI/6._r8*DIAM(1)**3*rhow/1000._r8  ! rhow kg/m3 --> g/cm3 
      radsl(1) = X(1)                         ! grams 
!      radsl(1) = X(1)/1000._r8       

      DO l=2,ncdp
         X(l)=2._r8*X(l-1)
         DIAM(l)=(6._r8/pi*X(l)*1000._r8/rhow)**(1._r8/3._r8)  ! cm
!         radsl(l)=X(l)/1000._r8 ! convert from g to kg
         radsl(l)=X(l)             
      ENDDO

! now get bin mid-points

      do l=1,ncd
         radl(l)=(radsl(l)+radsl(l+1))/2._r8         ! grams   
!         diammean(l) = (DIAM(l)+DIAM(l+1))/2._r8     ! cm
         diammean(l) = (6._r8/pi*radl(l)*1000._r8/rhow)**(1._r8/3._r8) ! cm
      end do

! set bin grid for method of moments

        ! for method of moments

      do lcl = 1,ncd+1
!         medge(lcl) = radsl(lcl)*1000._r8     ! convert to grams
         medge(lcl) = radsl(lcl)               ! grams
         diamedge(lcl) = DIAM(lcl)             ! cm
      enddo

      do lcl = 1,ncd
!         mmean(lcl) = radl(lcl)*1000._r8
         mmean(lcl) = radl(lcl)  
         diammean(lcl) = diammean(lcl)
      enddo

      do lcl = ncdp,1,-1
         if( diamedge(lcl).ge.40.e-4_r8 ) cutoff_id = lcl
      end do  

      write(*,*) 'cutoff_id = ', cutoff_id
         
end subroutine calc_bins

subroutine stochastic_kernel_init

    use cam_history_support, only: add_hist_coord

    integer :: idd, jdd
    real(r8) :: kkfac

    call calc_bins

    call add_hist_coord('bins_ncd', ncd, 'bins for TAU microphysics')

    call addfld('amk_c',(/'lev','bins_ncd'/),'A','kg','cloud liquid mass from bins')
    call addfld('ank_c',(/'lev','bins_ncd'/),'A','1/kg','cloud liquid number concentration from bins')
    call addfld('amk_r',(/'lev','bins_ncd'/),'A','kg','cloud liquid mass from bins')
    call addfld('ank_r',(/'lev','bins_ncd'/),'A','1/kg','cloud liquid number concentration from bins')
    call addfld('amk',(/'lev','bins_ncd'/),'A','kg','all liquid mass from bins')
    call addfld('ank',(/'lev','bins_ncd'/),'A','1/kg','all liquid number concentration from bins')
    call addfld('amk_out',(/'lev','bins_ncd'/),'A','kg','all liquid mass from bins')
    call addfld('ank_out',(/'lev','bins_ncd'/),'A','1/kg','all liquid number concentration from bins')

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

! Read in the collection kernel code from a lookup table. Again, this only needs to be done once.
! use kernel from Zach (who got it from Jerry)

     KNN(:,:)=0._r8 ! initialize values
     kkfac=1.5_r8   ! from Zach
!CACNOTE  - Need to fix the opening and reading of this file
     open(unit=40,file='/glade/u/home/cchen/TAU/input/KBARF',status='old')

 941 FORMAT(2X,E12.5)

     do idd=1,ncd
        do jdd=1,idd
   	  READ(40,941) KNN(IDD,JDD)
!     KNN(IDD,JDD)=(XK_GR(IDD)*kkfac+XK_GR(JDD)*kkfac)*KNN(IDD,JDD)
	  KNN(IDD,JDD)=(mmean(IDD)*kkfac+mmean(JDD)*kkfac)*KNN(IDD,JDD)

      	  if (knn(idd,jdd).lt.0._r8) knn(idd,jdd)=0._r8

     	end do
     end do

end subroutine stochastic_kernel_init

!main driver routine
!needs to pull in i,k fields (so might need dimensions here too)

subroutine stochastic_collect_tau_tend(deltatin, t,rho, qc, qr, qcin,ncin,qrin,nrin, lcldm, precip_frac, &
                                       mu_c, lambda_c, n0r, lambda_r, &
                                       qcin_new,ncin_new,qrin_new,nrin_new, &
                                       qctend,nctend,qrtend,nrtend,qctend_TAU,nctend_TAU,qrtend_TAU,nrtend_TAU, &
                                       scale_qc,scale_nc,scale_qr,scale_nr, &
                                       amk_c, ank_c, amk_r, ank_r, amk, ank, amk_out, ank_out, gmnnn_lmnnn_TAU, mgncol)


!use micro_pumas_utils, only: &
!       mg_liq_props, &
!       mg_rain_props

!use micro_pumasmg_utils, only: &
!       size_dist_param_liq, &
!       size_dist_param_basic

!inputs: mgncol,nlev,t,rho,qcin,ncin,qrin,nrin
!outputs: qctend,nctend,qrtend,nrtend
!not sure if we want to output bins (extra dimension). Good for testing?  

integer, intent(in) :: mgncol

real(r8), intent(in) :: deltatin
real(r8), intent(in) :: t(mgncol)
real(r8), intent(in) :: rho(mgncol)
real(r8), intent(in) :: qc(mgncol)
real(r8), intent(in) :: qr(mgncol)
real(r8), intent(in) :: qcin(mgncol)
real(r8), intent(in) :: ncin(mgncol)
real(r8), intent(in) :: qrin(mgncol)
real(r8), intent(in) :: nrin(mgncol)
real(r8), intent(in) :: lcldm(mgncol)
real(r8), intent(in) :: precip_frac(mgncol)
real(r8), intent(inout) :: qctend(mgncol)
real(r8), intent(inout) :: nctend(mgncol)
real(r8), intent(inout) :: qrtend(mgncol)
real(r8), intent(inout) :: nrtend(mgncol)
real(r8), intent(out) :: qctend_TAU(mgncol)
real(r8), intent(out) :: nctend_TAU(mgncol)
real(r8), intent(out) :: qrtend_TAU(mgncol)
real(r8), intent(out) :: nrtend_TAU(mgncol)

real(r8), intent(out) :: scale_qc(mgncol)
real(r8), intent(out) :: scale_nc(mgncol)
real(r8), intent(out) :: scale_qr(mgncol)
real(r8), intent(out) :: scale_nr(mgncol)

real(r8), intent(out) :: amk_c(mgncol,ncd)
real(r8), intent(out) :: ank_c(mgncol,ncd)
real(r8), intent(out) :: amk_r(mgncol,ncd)
real(r8), intent(out) :: ank_r(mgncol,ncd)
real(r8), intent(out) :: amk(mgncol,ncd)
real(r8), intent(out) :: ank(mgncol,ncd)
real(r8), intent(out) :: amk_out(mgncol,ncd)
real(r8), intent(out) :: ank_out(mgncol,ncd)

real(r8), intent(out) :: mu_c(mgncol)
real(r8), intent(out) :: lambda_c(mgncol)
real(r8), intent(out) :: lambda_r(mgncol)
real(r8), intent(out) :: n0r(mgncol)

real(r8) :: amk0(mgncol,ncd)
real(r8) :: ank0(mgncol,ncd)
real(r8) :: gnnnn(ncd)
real(r8) :: gmnnn(ncd)
real(r8) :: lnnnn(ncd)
real(r8) :: lmnnn(ncd)
real(r8) :: gnnnn0(ncd)
real(r8) :: gmnnn0(ncd)
real(r8) :: lnnnn0(ncd)
real(r8) :: lmnnn0(ncd)

real(r8), intent(out) :: qcin_new(mgncol)
real(r8), intent(out) :: ncin_new(mgncol)
real(r8), intent(out) :: qrin_new(mgncol)
real(r8), intent(out) :: nrin_new(mgncol)
real(r8), intent(out) :: gmnnn_lmnnn_TAU(mgncol)

real(r8) :: qcin_old(mgncol)
real(r8) :: ncin_old(mgncol)
real(r8) :: qrin_old(mgncol)
real(r8) :: nrin_old(mgncol)

integer :: i, n, lcl, cutoff_amk, cutoff(mgncol)

real(r8) :: all_gmnnn, all_lmnnn

integer, parameter :: sub_step = 60

cutoff = cutoff_id-1


!write(iulog,*) 'TAU time step = ', deltatin

!do k = 1,nlev
!call size_dist_param_liq(mg_liq_props, qcin, ncin, rho, mu_c, lambda_c, mgncol)
!call size_dist_param_basic(mg_rain_props, qrin, nrin, lambda_r, mgncol, n0=n0r)
!end do

! First make bins from cam size distribution (bins are diagnostic)
  
!call cam_bin_distribute(qcin,ncin,qrin,nrin,medge,amk,ank)  
do i=1,mgncol
!do k=1,nlev
   call cam_bin_distribute(qc(i), qr(i), qcin(i),ncin(i),qrin(i),nrin(i), &
                           mu_c(i),lambda_c(i),lambda_r(i),n0r(i), lcldm(i), precip_frac(i), &
                           scale_qc(i), scale_nc(i), scale_qr(i), scale_nr(i), &
                           amk_c(i,1:ncd),ank_c(i,1:ncd), amk_r(i,1:ncd), ank_r(i,1:ncd), amk(i,1:ncd), ank(i,1:ncd), cutoff_amk)
!end do
   if( cutoff_amk.gt.0 ) then
      cutoff(i) = cutoff_amk
   end if
!   cutoff(i) = cutoff_id-1
end do
  
!Then call the subroutines that actually do the calculations. The inputs/ouputs are described in comments below. 

!This part of the code needs to be called each time for each process rate calculation 
! (i.e., for each sampled cloud/rain gamma distribution):

! note: variables passed to compute_column_params are all inputs,
! outputs from this subroutine are stored as global variables

! inputs: t --> input air temperature (K)
!         rho --> input air density (kg/m^3)
!         medge --> bin mass boundary (g) 
!         amk --> array of bin mass mixing ratio, i.e., the input drop mass distribution (kg/kg)
!         ank --> array of bin number mixing ratio, i.e., the input drop number distribution (kg^-1)

! inputs: medge --> bin mass boundary (g), same as above

! outputs: gnnnn --> bin number mass mixing tendency gain, array in bins (#/cm^3/s)
!          lnnnn --> bin number mass mixing tendency loss, array in bins (#/cm^3/s)
!          gmnnn --> bin mass mixing ratio tendency gain, array in bins (g/cm^3/s) 
!          lmnnn --> bin mass mixing ratio tendency loss, array in bins (g/cm^3/s)


! Call Kernel

!do i=1,mgncol
!do k=1,nlev
!   call do_nn_n(gnnnn(i,:),gmnnn(i,:),lnnnn(i,:),lmnnn(i,:),medge)
!end do
!end do

qcin_new = 0._r8
ncin_new = 0._r8
qrin_new = 0._r8
nrin_new = 0._r8

qcin_old = 0._r8
ncin_old = 0._r8
qrin_old = 0._r8
nrin_old = 0._r8

qctend_TAU = 0._r8
nctend_TAU = 0._r8
qrtend_TAU = 0._r8
nrtend_TAU = 0._r8

amk0 = amk
ank0 = ank

! update qc, nc, qr, nr
do i=1,mgncol

gnnnn = 0._r8
gmnnn = 0._r8
lnnnn = 0._r8
lmnnn = 0._r8

! substep bin code
do n=1,sub_step
   call compute_coll_params(rho(i),medge,amk0(i,1:ncd),ank0(i,1:ncd),gnnnn0,gmnnn0,lnnnn0,lmnnn0)

   all_gmnnn=0._r8
   all_lmnnn=0._r8
! scaling gmnnn, lmnnn
   do lcl=1,ncd
      all_gmnnn = all_gmnnn+gmnnn0(lcl)
      all_lmnnn = all_lmnnn+lmnnn0(lcl)
   end do
 
   if( (all_gmnnn.eq.0._r8).or.(all_lmnnn.eq.0._r8) ) then
     gmnnn0(:) = 0._r8
     lmnnn0(:) = 0._r8
   else
     lmnnn0 = lmnnn0*(all_gmnnn/all_lmnnn)
   end if

   do lcl=1,ncd
      amk0(i,lcl) = amk0(i,lcl)+(gmnnn0(lcl)-lmnnn0(lcl))*1.e3_r8/rho(i)*deltatin/dble(sub_step)
      ank0(i,lcl) = ank0(i,lcl)+(gnnnn0(lcl)-lnnnn0(lcl))*1.e6_r8/rho(i)*deltatin/dble(sub_step)

      gmnnn(lcl) = gmnnn(lcl)+gmnnn0(lcl)/sub_step
      gnnnn(lcl) = gnnnn(lcl)+gnnnn0(lcl)/sub_step
      lmnnn(lcl) = lmnnn(lcl)+lmnnn0(lcl)/sub_step
      lnnnn(lcl) = lnnnn(lcl)+lnnnn0(lcl)/sub_step
   end do

end do

   ! cloud water
   do lcl = 1,cutoff(i)
      qcin_old(i) = qcin_old(i)+amk(i,lcl)
      ncin_old(i) = ncin_old(i)+ank(i,lcl)
      qcin_new(i) = qcin_new(i)+(gmnnn(lcl)-lmnnn(lcl))*1.e3_r8/rho(i)*deltatin
      ncin_new(i) = ncin_new(i)+(gnnnn(lcl)-lnnnn(lcl))*1.e6_r8/rho(i)*deltatin

!      qctend_TAU(i) = qctend_TAU(i)+(gmnnn(lcl)-lmnnn(lcl))*1.e3_r8/rho(i)
!      nctend_TAU(i) = nctend_TAU(i)+(gnnnn(lcl)-lnnnn(lcl))*1.e6_r8/rho(i)
      qctend_TAU(i) = qctend_TAU(i)+(amk0(i,lcl)-amk(i,lcl))/deltatin
      nctend_TAU(i) = nctend_TAU(i)+(ank0(i,lcl)-ank(i,lcl))/deltatin

      gmnnn_lmnnn_TAU(i) = gmnnn_lmnnn_TAU(i)+gmnnn(lcl)-lmnnn(lcl)
   end do
   ! rain 
   do lcl = cutoff(i)+1, ncd
      qrin_old(i) = qrin_old(i)+amk(i,lcl)
      nrin_old(i) = nrin_old(i)+ank(i,lcl)
      qrin_new(i) = qrin_new(i)+(gmnnn(lcl)-lmnnn(lcl))*1.e3_r8/rho(i)*deltatin
      nrin_new(i) = nrin_new(i)+(gnnnn(lcl)-lnnnn(lcl))*1.e6_r8/rho(i)*deltatin

!      qrtend_TAU(i) = qrtend_TAU(i)+(gmnnn(lcl)-lmnnn(lcl))*1.e3_r8/rho(i)
!      nrtend_TAU(i) = nrtend_TAU(i)+(gnnnn(lcl)-lnnnn(lcl))*1.e6_r8/rho(i)
      qrtend_TAU(i) = qrtend_TAU(i)+(amk0(i,lcl)-amk(i,lcl))/deltatin
      nrtend_TAU(i) = nrtend_TAU(i)+(ank0(i,lcl)-ank(i,lcl))/deltatin
      gmnnn_lmnnn_TAU(i) = gmnnn_lmnnn_TAU(i)+gmnnn(lcl)-lmnnn(lcl)
   end do

   do lcl = 1,ncd
      amk_out(i,lcl) = amk(i,lcl) + (gmnnn(lcl)-lmnnn(lcl))*1.e3_r8/rho(i)*deltatin
      ank_out(i,lcl) = ank(i,lcl) + (gnnnn(lcl)-lnnnn(lcl))*1.e6_r8/rho(i)*deltatin
   end do

   qcin_new(i) = qcin_new(i)+qcin_old(i)
   ncin_new(i) = ncin_new(i)+ncin_old(i)
   qrin_new(i) = qrin_new(i)+qrin_old(i)
   nrin_new(i) = nrin_new(i)+nrin_old(i)
end do


end subroutine stochastic_collect_tau_tend


subroutine cam_bin_distribute(qc_all, qr_all, qc,nc,qr,nr,mu_c,lambda_c,lambda_r,n0r, &
                              lcldm, precip_frac, scale_qc, scale_nc, scale_qr, scale_nr, &
                              amk_c,ank_c, amk_r, ank_r, amk, ank, cutoff_amk)

real(r8) :: qc_all, qr_all, qc, nc, qr, nr, mu_c, lambda_c, lambda_r, n0r, lcldm, precip_frac
real(r8), dimension(ncd) :: amk_c, ank_c, amk_r, ank_r, amk, ank
integer  :: i
real(r8) :: phi
real(r8) :: scale_nc, scale_qc, scale_nr, scale_qr 

integer  :: id_max_qc, id_max_qr, cutoff_amk
real(r8) :: max_qc, max_qr, min_amk 

ank_c = 0._r8
amk_c = 0._r8
ank_r = 0._r8
amk_r = 0._r8
ank = 0._r8
amk = 0._r8

scale_nc = 0._r8
scale_qc = 0._r8
scale_nr = 0._r8
scale_qr = 0._r8

id_max_qc = 0
id_max_qr = 0
cutoff_amk = 0
max_qc = 0._r8
max_qr = 0._r8

! cloud water, nc in #/m3 --> #/cm3
if( (qc_all.gt.qsmall).and.(qc.gt.qsmall) ) then
do i = 1,ncd
   phi = nc*lambda_c**(mu_c+1._r8)/gamma(mu_c+1._r8)*(diammean(i)*1.e-2_r8)**mu_c*exp(-lambda_c*diammean(i)*1.e-2_r8) ! D cm --> m
   ank_c(i) = phi*(diamedge(i+1)-diamedge(i))*1.e-2_r8                     ! D cm --> m
   amk_c(i) = phi*(diamedge(i+1)-diamedge(i))*1.e-2_r8*mmean(i)*1.e-3_r8   ! mass in bin g --> kg 

   scale_nc = scale_nc+ank_c(i)
   scale_qc = scale_qc+amk_c(i) 
end do

scale_nc = scale_nc/nc
scale_qc = scale_qc/qc

ank_c = ank_c/scale_nc*lcldm
amk_c = amk_c/scale_qc*lcldm
!ank_c = ank_c*lcldm
!amk_c = amk_c*lcldm

do i=1,ncd
   if( amk_c(i).gt.max_qc ) then
      id_max_qc = i
      max_qc = amk_c(i)
   end if
end do

!else

!do i=1,ncd
!   ank_c(i) = 0._r8
!   amk_c(i) = 0._r8
!end do

end if

! rain drop
if( (qr_all.gt.qsmall).and.(qr.gt.qsmall) ) then
do i = 1,ncd
   phi = n0r*exp(-lambda_r*diammean(i)*1.e-2_r8)                   ! D cm --> m
   ank_r(i) = phi*(diamedge(i+1)-diamedge(i))*1.e-2_r8    ! D cm --> m  
   amk_r(i) = phi*(diamedge(i+1)-diamedge(i))*1.e-2_r8*mmean(i)*1.e-3_r8

   scale_nr = scale_nr + ank_r(i)
   scale_qr = scale_qr + amk_r(i)
end do

scale_nr = scale_nr/nr
scale_qr = scale_qr/qr

ank_r = ank_r/scale_nr*precip_frac
amk_r = amk_r/scale_qr*precip_frac
!ank_r = ank_r*precip_frac
!amk_r = amk_r*precip_frac

!else

!do i=1,ncd
!   ank_r(i) = 0._r8
!   amk_r(i) = 0._r8
!end do

do i=1,ncd
   if( amk_r(i).gt.max_qr ) then
      id_max_qr = i
      max_qr = amk_r(i)
   end if
end do

end if

amk = amk_c + amk_r
ank = ank_c + ank_r


if( (id_max_qc.gt.0).and.(id_max_qr.gt.0) ) then
   if( (max_qc/max_qr.lt.10._r8).or.(max_qc/max_qr.gt.0.1_r8) )then
      min_amk = amk(id_max_qc)

      do i=id_max_qc,id_max_qr
         if( amk(i).le.min_amk ) then
           cutoff_amk = i
           min_amk = amk(i)
         end if
      end do
   end if
end if


!if( qc_all.gt.qsmall.OR.qr_all.gt.qsmall ) then
!   do i=1,ncd
!      ank(i) = ank_c(i) + ank_r(i)
!      amk(i) = amk_c(i) + amk_r(i)
!   end do
!else
!   do i=1,ncd
!      amk(i) = 0._r8
!      ank(i) = 0._r8
!   end do
!end if
!input: qc,nc,qr,nr, medge (bin edges). May also need # bins?
!output: amk, ank (mixing ratio and number in each bin)

!this part will take a bit of thinking about.
!use size distribution parameters (mu, lambda) to generate the values at discrete size points
!need to also ensure mass conservation  
  
  
end subroutine cam_bin_distribute

         
! here are the subroutines called above that actually do the collision-coalescence calculations:

! The Kernel is from Jerry from many moons ago (included)

! I read in the file data and multiply by the summed mass of the individual bins 
! (with a factor of 1.5 so that the values represent the middle of the bin

! 941 FORMAT(2X,E12.5)
!     READ(40,941) KNN(IDD,JDD)
!     KNN(IDD,JDD)=(XK_GR(IDD)*kkfac+XK_GR(JDD)*kkfac)*KNN(IDD,JDD)

!where idd and jdd are the indexes for the bins and xk_gr is the mass of drops in a bin in grams
!

!************************************************************************************
! Setup variables needed for collection
! Either pass in or define globally the following variables
! tbase(height) - temperature in K as a function of height
! rhon(height) - air density as a function of height in kg/m^3
! xk_gr(bins) - mass of single drop in each bin in grams
! lsmall - small number
! QC - mass mixing ratio in kg/kg
! QN - number mixing ratio in #/kg
! All parameters are defined to be global in my version so that they are readily available throughout the code:
! SMN0,SNN0,SMCN,APN,AMN2,AMN3,PSIN,FN,FPSIN,XPSIN,HPSIN,FN2,XXPSIN (all arrays of drop bins)
!************************************************************************************

!AG: Global arrays need to be passed around I think? Right now at the module level. Is that okay?

SUBROUTINE COMPUTE_COLL_PARAMS(rhon,xk_gr,qc,qn,gnnnn,gmnnn,lnnnn,lmnnn)
  IMPLICIT NONE

! variable declarations (added by hm, 020118)
! note: vertical array looping is stripped out, this subroutine operates
! only on LOCAL values

  real(r8), dimension(ncd) :: qc,qn
  real(r8), dimension(ncdp) :: xk_gr
  real(r8) :: tbase,rhon
!  real(r8) :: TAIRC,UMMS,UMMS2
  integer :: lk
  integer :: l
  real(r8), parameter :: lsmall = 1.e-12_r8
  real(r8), dimension(ncd) :: smn0,snn0,smcn,amn2,amn3,psin,fn,fpsin, &
                               xpsin,hpsin,fn2,xxpsin
  real(r8) :: apn

  real(r8), dimension(ncd) :: gnnnn,gmnnn,lnnnn,lmnnn
  integer :: lm1,ll

  lk=ncd



!....................................................................................
!  TAIRC=TBASE(K)-273.15
!  TAIRC=TBASE-273.15_r8
!  UMMS=UMM(TAIRC)
!!  UMMS2=UMMS*4.66/(RHON(K)/1.E3)
!!  UMMS=UMMS/(RHON(K)/1.E3)
!  UMMS2=UMMS*4.66_r8/(RHON/1.E3_r8)
!  UMMS=UMMS/(RHON/1.E3_r8)

  DO L=1,LK
!     SMN0(L)=QC(L,K)*RHON(K)/1.E3
!     SNN0(L)=QN(L,K)*RHON(K)/1.E6
     SMN0(L)=QC(L)*RHON/1.E3_r8
     SNN0(L)=QN(L)*RHON/1.E6_r8

     IF(SMN0(L).LT.lsmall.OR.SNN0(L).LT.lsmall)THEN
        SMN0(L)=0.0_r8
        SNN0(L)=0.0_r8
     ENDIF
  ENDDO

  DO L=1,LK
     IF(SMN0(L) .gt. 0._r8.AND.SNN0(L) .gt. 0._r8)THEN
        SMCN(L)=SMN0(L)/SNN0(L)
        IF((SMCN(L) .GT. 2._r8*XK_GR(L)))THEN
!           SMCN(L) = (2*XK_GR(L))
           SMCN(L) = (2._r8*XK_GR(L))
        ENDIF
        IF((SMCN(L) .LT. XK_GR(L)))THEN
           SMCN(L) = XK_GR(L)
        ENDIF
     ELSE
        SMCN(L)=0._r8
     ENDIF
     IF (SMCN(L).LT.XK_GR(L).OR.SMCN(L).GT.(2._r8*XK_GR(L)).OR.SMCN(L).EQ.0.0_r8)THEN
        APN=1.0_r8
     ELSE
!        APN=0.5*(1.+3.*(XK_GR(L)/SMCN(L))-2*((XK_GR(L)/SMCN(L))**2.))
        APN=0.5_r8*(1._r8+3._r8*(XK_GR(L)/SMCN(L))-2._r8*((XK_GR(L)/SMCN(L))**2._r8))
     ENDIF

     IF(SNN0(L) .GT. LSMALL)THEN
        AMN2(L)=APN*SMN0(L)*SMN0(L)/SNN0(L)
        AMN3(L)=APN*APN*APN*SMN0(L)*SMN0(L)*SMN0(L)/(SNN0(L)*SNN0(L))
     ELSE
        AMN2(L)=0._r8
        AMN3(L)=0._r8
     ENDIF

     IF(SMCN(L).LT.XK_GR(L))THEN
        PSIN(L)=0.0_r8
        FN(L)=2._r8*SNN0(L)/XK_GR(L)
     ELSE
        IF(SMCN(L).GT.(2._r8*XK_GR(L)))THEN
           FN(L)=0.0_r8
           PSIN(L)=2._r8*SNN0(L)/XK_GR(L)
        ELSE
           PSIN(L)=2._r8/XK_GR(L)*(SMN0(L)/XK_GR(L)-SNN0(L))
           FN(L)=2._r8/XK_GR(L)*(2._r8*SNN0(L)-SMN0(L)/XK_GR(L))
        ENDIF
     ENDIF

     IF(SNN0(L).LT.LSMALL.OR.SMN0(L).LT.LSMALL)THEN
        PSIN(L)=0.0_r8
        FN(L)=0.0_r8
     ENDIF

     FPSIN(L)=0.5_r8/XK_GR(L)*(PSIN(L)-FN(L))
     XPSIN(L)=2._r8*XK_GR(L)*PSIN(L)
     HPSIN(L)=PSIN(L)-0.5_r8*FN(L)
     FN2(L)=FN(L)/2._r8

     IF(L.GT.1)THEN
        XXPSIN(L)=XK_GR(L)*PSIN(L-1)
     ENDIF
  ENDDO

!************************************************************************************
! Compute collision coalescence
! Either pass in or define globally the following variables
! Gain terms begin with G, loss terms begin with L
! Second letter defines mass (M) or number (N)
! Third and fourth letters define the types of particles colling, i.e., NN means drops colliding with drops
! Last letter defines the category the new particles go into, in this case just N for liquid drops
! The resulting rates are in units of #/cm^3/s and g/cm^3/s
! Relies on predefined kernel array KNN(bins,bins) - see top of this file
!************************************************************************************

   GMNNN = 0._r8
   GNNNN = 0._r8
   LMNNN = 0._r8
   LNNNN = 0._r8
! remove verical array index, calculate gain/loss terms locally

  DO L=3,LK-1
     LM1=L-1
     DO LL=1,L-2
!        GNNNN(L,K)=GNNNN(L,K)+(PSIN(LM1)*SMN0(LL)-FPSIN(LM1)*AMN2(LL))*KNN(LM1,LL)
!        GMNNN(L,K)=GMNNN(L,K)+(XK_GR(L)*PSIN(LM1)*SMN0(LL)+FN2(LM1)*AMN2(LL)-FPSIN(LM1)*AMN3(LL))*KNN(LM1,LL)
        GNNNN(L)=GNNNN(L)+(PSIN(LM1)*SMN0(LL)-FPSIN(LM1)*AMN2(LL))*KNN(LM1,LL)
        GMNNN(L)=GMNNN(L)+(XK_GR(L)*PSIN(LM1)*SMN0(LL)+FN2(LM1)*AMN2(LL)-FPSIN(LM1)*AMN3(LL))*KNN(LM1,LL)
     ENDDO
  ENDDO

  DO L=2,LK-1
     LM1=L-1
     GNNNN(L)=GNNNN(L)+0.5_r8*SNN0(LM1)*SNN0(LM1)*KNN(LM1,LM1)
     GMNNN(L)=GMNNN(L)+0.5_r8*(SNN0(LM1)*SMN0(LM1)+SMN0(LM1)*SNN0(LM1))*KNN(LM1,LM1)
     DO LL=1,L-1
!        LNNNN(L,K)=LNNNN(L,K)+(PSIN(L)*SMN0(LL)-FPSIN(L)*AMN2(LL))*KNN(L,LL)
!        GMNNN(L,K)=GMNNN(L,K)+(SMN0(LL)*SNN0(L)-PSIN(L)*AMN2(LL)+FPSIN(L)*AMN3(LL))*KNN(L,LL)
!        LMNNN(L,K)=LMNNN(L,K)+(XPSIN(L)*SMN0(LL)-HPSIN(L)*AMN2(LL))*KNN(L,LL)
        LNNNN(L)=LNNNN(L)+(PSIN(L)*SMN0(LL)-FPSIN(L)*AMN2(LL))*KNN(L,LL)
        GMNNN(L)=GMNNN(L)+(SMN0(LL)*SNN0(L)-PSIN(L)*AMN2(LL)+FPSIN(L)*AMN3(LL))*KNN(L,LL)
        LMNNN(L)=LMNNN(L)+(XPSIN(L)*SMN0(LL)-HPSIN(L)*AMN2(LL))*KNN(L,LL)
     ENDDO
  ENDDO

  DO L=1,LK-1
     DO LL=L,LK-1
!        LNNNN(L,K)=LNNNN(L,K)+(SNN0(LL)*SNN0(L))*KNN(LL,L)
!        LMNNN(L,K)=LMNNN(L,K)+(SNN0(LL)*SMN0(L))*KNN(LL,L)
        LNNNN(L)=LNNNN(L)+(SNN0(LL)*SNN0(L))*KNN(LL,L)
        LMNNN(L)=LMNNN(L)+(SNN0(LL)*SMN0(L))*KNN(LL,L)
     ENDDO
  ENDDO




END SUBROUTINE COMPUTE_COLL_PARAMS




end module stochastic_collect_tau_cam


