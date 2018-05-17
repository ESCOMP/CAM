module crmx_sgs

! module for original SAM subgrid-scale SGS closure (Smagorinsky or 1st-order TKE)
! Marat Khairoutdinov, 2012

use crmx_grid, only: nx,nxp1,ny,nyp1,YES3D,nzm,nz,dimx1_s,dimx2_s,dimy1_s,dimy2_s 
use crmx_params, only: dosgs
use crmx_vars, only: tke2, tk2
implicit none

!----------------------------------------------------------------------
! Required definitions:

!!! prognostic scalar (need to be advected arround the grid):

integer, parameter :: nsgs_fields = 1   ! total number of prognostic sgs vars

real sgs_field(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm, nsgs_fields)

!!! sgs diagnostic variables that need to exchange boundary information (via MPI):

integer, parameter :: nsgs_fields_diag = 2   ! total number of diagnostic sgs vars

! diagnostic fields' boundaries:
integer, parameter :: dimx1_d=0, dimx2_d=nxp1, dimy1_d=1-YES3D, dimy2_d=nyp1

real sgs_field_diag(dimx1_d:dimx2_d, dimy1_d:dimy2_d, nzm, nsgs_fields_diag)

logical:: advect_sgs = .false. ! advect prognostics or not, default - not (Smagorinsky)
logical, parameter:: do_sgsdiag_bound = .true.  ! exchange boundaries for diagnostics fields

! SGS fields that output by default (if =1).
integer, parameter :: flag_sgs3Dout(nsgs_fields) = (/0/)
integer, parameter :: flag_sgsdiag3Dout(nsgs_fields_diag) = (/0,0/)

real fluxbsgs (nx, ny, 1:nsgs_fields) ! surface fluxes 
real fluxtsgs (nx, ny, 1:nsgs_fields) ! top boundary fluxes 

!!! these arrays may be needed for output statistics:

real sgswle(nz,1:nsgs_fields)  ! resolved vertical flux
real sgswsb(nz,1:nsgs_fields)  ! SGS vertical flux
real sgsadv(nz,1:nsgs_fields)  ! tendency due to vertical advection
real sgslsadv(nz,1:nsgs_fields)  ! tendency due to large-scale vertical advection
real sgsdiff(nz,1:nsgs_fields)  ! tendency due to vertical diffusion

!------------------------------------------------------------------
! internal (optional) definitions:

! make aliases for prognostic variables:

real tke(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)   ! SGS TKE
equivalence (tke(dimx1_s,dimy1_s,1),sgs_field(dimx1_s,dimy1_s,1,1))

! make aliases for diagnostic variables:

real tk  (dimx1_d:dimx2_d, dimy1_d:dimy2_d, nzm) ! SGS eddy viscosity
real tkh (dimx1_d:dimx2_d, dimy1_d:dimy2_d, nzm) ! SGS eddy conductivity
equivalence (tk(dimx1_d,dimy1_d,1), sgs_field_diag(dimx1_d, dimy1_d,1,1))
equivalence (tkh(dimx1_d,dimy1_d,1), sgs_field_diag(dimx1_d, dimy1_d,1,2))


real grdf_x(nzm)! grid factor for eddy diffusion in x
real grdf_y(nzm)! grid factor for eddy diffusion in y
real grdf_z(nzm)! grid factor for eddy diffusion in z

logical:: dosmagor   ! if true, then use Smagorinsky closure

! Local diagnostics:

real tkesbbuoy(nz), tkesbshear(nz),tkesbdiss(nz), tkesbdiff(nz)

CONTAINS

! required microphysics subroutines and function:
!----------------------------------------------------------------------
!!! Read microphysics options from prm (namelist) file

subroutine sgs_setparm()

  use crmx_grid, only: case
  implicit none

  integer ierr, ios, ios_missing_namelist, place_holder

  !======================================================================
  NAMELIST /SGS_TKE/ &
       dosmagor ! Diagnostic Smagorinsky closure

  NAMELIST /BNCUIODSBJCB/ place_holder

  dosmagor = .true. ! default 

  !----------------------------------
  !  Read namelist for microphysics options from prm file:
  !------------
  !open(55,file='./'//trim(case)//'/prm', status='old',form='formatted')

  !read (UNIT=55,NML=BNCUIODSBJCB,IOSTAT=ios_missing_namelist)
  !rewind(55) !note that one must rewind before searching for new namelists

  !read (55,SGS_TKE,IOSTAT=ios)

  advect_sgs = .not.dosmagor

  !if (ios.ne.0) then
  !   !namelist error checking
  !   if(ios.ne.ios_missing_namelist) then
  !      write(*,*) '****** ERROR: bad specification in SGS_TKE namelist'
  !      call task_abort()
  !   end if
  !end if
  !close(55)

  ! END UW ADDITION
  !======================================================================

end subroutine sgs_setparm

!----------------------------------------------------------------------
!!! Initialize sgs:


subroutine sgs_init()

  use crmx_grid, only: nrestart, dx, dy, dz, adz, masterproc
  use crmx_params, only: LES
  integer k

  if(nrestart.eq.0) then

     sgs_field = 0.
     sgs_field_diag = 0.

     fluxbsgs = 0.
     fluxtsgs = 0.

  end if

!  if(masterproc) then
!     if(dosmagor) then
!        write(*,*) 'Smagorinsky SGS Closure'
!     else
!        write(*,*) 'Prognostic TKE 1.5-order SGS Closure'
!     end if
!  end if

  if(LES) then
    do k=1,nzm
       grdf_x(k) = dx**2/(adz(k)*dz)**2
       grdf_y(k) = dy**2/(adz(k)*dz)**2
       grdf_z(k) = 1.
    end do
  else
    do k=1,nzm
       grdf_x(k) = min(16.,dx**2/(adz(k)*dz)**2)
       grdf_y(k) = min(16.,dy**2/(adz(k)*dz)**2)
       grdf_z(k) = 1.
    end do
  end if

  sgswle = 0.
  sgswsb = 0.
  sgsadv = 0.
  sgsdiff = 0.
  sgslsadv = 0.


end subroutine sgs_init

!----------------------------------------------------------------------
!!! make some initial noise in sgs:
!
subroutine setperturb_sgs(ptype)

use crmx_vars, only: q0, z
integer, intent(in) :: ptype
integer i,j,k

select case (ptype)

  case(0)

     do k=1,nzm
      do j=1,ny
       do i=1,nx
         if(k.le.4.and..not.dosmagor) then
            tke(i,j,k)=0.04*(5-k)
         endif
       end do
      end do
     end do

  case(1)

     do k=1,nzm
      do j=1,ny
       do i=1,nx
         if(q0(k).gt.6.e-3.and..not.dosmagor) then
            tke(i,j,k)=1.
         endif
       end do
      end do
     end do

  case(2)

  case(3)   ! gcss wg1 smoke-cloud case

     do k=1,nzm
      do j=1,ny
       do i=1,nx
         if(q0(k).gt.0.5e-3.and..not.dosmagor) then
            tke(i,j,k)=1.
         endif
       end do
      end do
     end do


  case(4)  ! gcss wg1 arm case

     do k=1,nzm
      do j=1,ny
       do i=1,nx
         if(z(k).le.150..and..not.dosmagor) then
            tke(i,j,k)=0.15*(1.-z(k)/150.)
         endif
       end do
      end do
     end do


  case(5)  ! gcss wg1 BOMEX case

     do k=1,nzm
      do j=1,ny
       do i=1,nx
         if(z(k).le.3000..and..not.dosmagor) then
            tke(i,j,k)=1.-z(k)/3000.
         endif
       end do
      end do
     end do

  case(6)  ! GCSS Lagragngian ASTEX


     do k=1,nzm
      do j=1,ny
       do i=1,nx
         if(q0(k).gt.6.e-3.and..not.dosmagor) then
            tke(i,j,k)=1.
         endif
       end do
      end do
     end do


  case default

end select

end subroutine setperturb_sgs

!----------------------------------------------------------------------
!!! Estimate Courant number limit for SGS
!

subroutine kurant_sgs(cfl)

use crmx_grid, only: dt, dx, dy, dz, adz, adzw
implicit none

real, intent(out) :: cfl

integer k
real tkhmax(nz)

do k = 1,nzm
 tkhmax(k) = maxval(tkh(1:nx,1:ny,k))
end do

cfl = 0.
do k=1,nzm
  cfl = max(cfl,        &
     0.5*tkhmax(k)*grdf_z(k)*dt/(dz*adzw(k))**2, &
     0.5*tkhmax(k)*grdf_x(k)*dt/dx**2, &
     YES3D*0.5*tkhmax(k)*grdf_y(k)*dt/dy**2)
end do

end subroutine kurant_sgs


!----------------------------------------------------------------------
!!! compute sgs diffusion of momentum:
!
subroutine sgs_mom()

   call diffuse_mom()

end subroutine sgs_mom

!----------------------------------------------------------------------
!!! compute sgs diffusion of scalars:
!
subroutine sgs_scalars()

  use crmx_vars
  use crmx_microphysics
  use crmx_crmtracers
  use crmx_params, only: dotracers
  implicit none

    real dummy(nz)
    real fluxbtmp(nx,ny), fluxttmp(nx,ny) !bloss
    integer k


      call diffuse_scalar(t,fluxbt,fluxtt,tdiff,twsb, &
                           t2lediff,t2lediss,twlediff,.true.)
    
      if(advect_sgs) then
         call diffuse_scalar(tke,fzero,fzero,dummy,sgswsb, &
                                    dummy,dummy,dummy,.false.)
      end if


!
!    diffusion of microphysics prognostics:
!
      call micro_flux()

      total_water_evap = total_water_evap - total_water()

      do k = 1,nmicro_fields
        if(   k.eq.index_water_vapor             &! transport water-vapor variable no metter what
         .or. docloud.and.flag_precip(k).ne.1    & ! transport non-precipitation vars
         .or. doprecip.and.flag_precip(k).eq.1 ) then
           fluxbtmp(1:nx,1:ny) = fluxbmk(1:nx,1:ny,k)
           fluxttmp(1:nx,1:ny) = fluxtmk(1:nx,1:ny,k)
           call diffuse_scalar(micro_field(:,:,:,k),fluxbtmp,fluxttmp, &
                mkdiff(:,k),mkwsb(:,k), dummy,dummy,dummy,.false.)
       end if
      end do

      total_water_evap = total_water_evap + total_water()

 ! diffusion of tracers:

      if(dotracers) then

        call tracers_flux()

        do k = 1,ntracers

          fluxbtmp = fluxbtr(:,:,k)
          fluxttmp = fluxttr(:,:,k)
          call diffuse_scalar(tracer(:,:,:,k),fluxbtmp,fluxttmp, &
               trdiff(:,k),trwsb(:,k), &
               dummy,dummy,dummy,.false.)
!!$          call diffuse_scalar(tracer(:,:,:,k),fluxbtr(:,:,k),fluxttr(:,:,k),trdiff(:,k),trwsb(:,k), &
!!$                           dummy,dummy,dummy,.false.)

        end do

      end if



end subroutine sgs_scalars

!----------------------------------------------------------------------
!!! compute sgs processes (beyond advection):
!
subroutine sgs_proc()

   use crmx_grid, only: nstep,dt,icycle
   use crmx_params, only: dosmoke

!    SGS TKE equation:

     if(dosgs) call tke_full()

     tke2 = tke
     tk2 = tk

end subroutine sgs_proc

!----------------------------------------------------------------------
!!! Diagnose arrays nessesary for dynamical core and statistics:
!
subroutine sgs_diagnose()
! None 

end subroutine sgs_diagnose

!----------------------------------------------------------------------
! called when stepout() called

subroutine sgs_print()

 call fminmax_print('tke:',tke,dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm)
 call fminmax_print('tk:',tk,0,nxp1,1-YES3D,nyp1,nzm)
 call fminmax_print('tkh:',tkh,0,nxp1,1-YES3D,nyp1,nzm)

end subroutine sgs_print

!----------------------------------------------------------------------
!!! Initialize the list of sgs statistics 
!
subroutine sgs_hbuf_init(namelist,deflist,unitlist,status,average_type,count,sgscount)
character(*) namelist(*), deflist(*), unitlist(*)
integer status(*),average_type(*),count,sgscount

end subroutine sgs_hbuf_init


end module crmx_sgs



