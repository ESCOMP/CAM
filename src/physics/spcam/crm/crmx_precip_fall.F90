subroutine precip_fall(qp, term_vel, hydro_type, omega, ind)

!     positively definite monotonic advection with non-oscillatory option
!     and gravitational sedimentation 

use crmx_vars
use crmx_params
implicit none



real qp(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm) ! falling hydrometeor
integer hydro_type   ! 0 - all liquid, 1 - all ice, 2 - mixed
real omega(nx,ny,nzm)   !  = 1: liquid, = 0: ice;  = 0-1: mixed : used only when hydro_type=2
integer ind

! Terminal velocity fnction 	

real, external :: term_vel  ! terminal velocity function


! Local:

real mx(nzm),mn(nzm), lfac(nz)
real www(nz),fz(nz)
real df(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)
real f0(nzm),df0(nzm)
real eps
integer i,j,k,kc,kb
logical nonos

real y,pp,pn
pp(y)= max(0.,y)
pn(y)=-min(0.,y)

real lat_heat, wmax

real wp(nzm), tmp_qp(nzm), irhoadz(nzm), iwmax(nzm), rhofac(nzm), prec_cfl
integer nprec, iprec
real  flagstat

!--------------------------------------------------------

!call t_startf ('precip_fall')

eps = 1.e-10
nonos = .true.
  
 do k = 1,nzm
    rhofac(k) = sqrt(1.29/rho(k))
    irhoadz(k) = 1./(rho(k)*adz(k)) ! Useful factor
    kb = max(1,k-1)
    wmax       = dz*adz(kb)/dtn   ! Velocity equivalent to a cfl of 1.0.
    iwmax(k)   = 1./wmax
 end do

! 	Add sedimentation of precipitation field to the vert. vel.

do j=1,ny
   do i=1,nx

      ! Compute precipitation velocity and flux column-by-column
      
      prec_cfl = 0.

      do k=1,nzm

         select case (hydro_type)
         case(0) 
            lfac(k) = fac_cond 
            flagstat = 1.
         case(1) 
            lfac(k) = fac_sub
            flagstat = 1.
         case(2)
            lfac(k) = fac_cond + (1-omega(i,j,k))*fac_fus
            flagstat = 1.
         case(3)
            lfac(k) = 0.
            flagstat = 0.
         case default
           if(masterproc) then
             print*, 'unknown hydro_type in precip_fall. exitting ...'
             call task_abort
           end if
         end select

         wp(k)=rhofac(k)*term_vel(i,j,k,ind)
         prec_cfl = max(prec_cfl,wp(k)*iwmax(k)) ! Keep column maximum CFL
         wp(k) = -wp(k)*rhow(k)*dtn/dz 

      end do  ! k

      fz(nz)=0.
      www(nz)=0.
      lfac(nz)=0

      ! If maximum CFL due to precipitation velocity is greater than 0.9,
      ! take more than one advection step to maintain stability.
      if (prec_cfl.gt.0.9) then
         nprec = CEILING(prec_cfl/0.9)
         do k = 1,nzm
            ! wp already includes factor of dt, so reduce it by a
            ! factor equal to the number of precipitation steps.
            wp(k) = wp(k)/float(nprec) 
         end do
      else
         nprec = 1
      end if

      do iprec = 1,nprec

         do k = 1,nzm
            tmp_qp(k) = qp(i,j,k) ! Temporary array for qp in this column
         end do

         !-----------------------------------------

         if(nonos) then

            do k=1,nzm
               kc=min(nzm,k+1)
               kb=max(1,k-1)
               mx(k)=max(tmp_qp(kb),tmp_qp(kc),tmp_qp(k))
               mn(k)=min(tmp_qp(kb),tmp_qp(kc),tmp_qp(k))	  
            end do

         end if  ! nonos

         !  loop over iterations

         do k=1,nzm
            ! Define upwind precipitation flux
            fz(k)=tmp_qp(k)*wp(k)
         end do

         do k=1,nzm
            kc=k+1
            tmp_qp(k)=tmp_qp(k)-(fz(kc)-fz(k))*irhoadz(k) !Update temporary qp
         end do

         do k=1,nzm
            ! Also, compute anti-diffusive correction to previous
            ! (upwind) approximation to the flux
            kb=max(1,k-1)
            ! The precipitation velocity is a cell-centered quantity,
            ! since it is computed from the cell-centered
            ! precipitation mass fraction.  Therefore, a reformulated
            ! anti-diffusive flux is used here which accounts for
            ! this and results in reduced numerical diffusion.
            www(k) = 0.5*(1.+wp(k)*irhoadz(k)) &
                 *(tmp_qp(kb)*wp(kb) - tmp_qp(k)*wp(k)) ! works for wp(k)<0
         end do

         !---------- non-osscilatory option ---------------

         if(nonos) then

            do k=1,nzm
               kc=min(nzm,k+1)
               kb=max(1,k-1)
               mx(k)=max(tmp_qp(kb),tmp_qp(kc),tmp_qp(k),mx(k))
               mn(k)=min(tmp_qp(kb),tmp_qp(kc),tmp_qp(k),mn(k))	  
            end do

            do k=1,nzm
               kc=min(nzm,k+1)
               mx(k)=rho(k)*adz(k)*(mx(k)-tmp_qp(k))/(pn(www(kc)) + pp(www(k))+eps)
               mn(k)=rho(k)*adz(k)*(tmp_qp(k)-mn(k))/(pp(www(kc)) + pn(www(k))+eps)
            end do

            do k=1,nzm
               kb=max(1,k-1)
               ! Add limited flux correction to fz(k).
               fz(k) = fz(k) &                        ! Upwind flux
                    + pp(www(k))*min(1.,mx(k), mn(kb)) &
                    - pn(www(k))*min(1.,mx(kb),mn(k)) ! Anti-diffusive flux
            end do

         endif ! nonos

         ! Update precipitation mass fraction and liquid-ice static
         ! energy using precipitation fluxes computed in this column.
         do k=1,nzm
            kc=k+1
            ! Update precipitation mass fraction.
            ! Note that fz is the total flux, including both the
            ! upwind flux and the anti-diffusive correction.
            qp(i,j,k)=qp(i,j,k)-(fz(kc)-fz(k))*irhoadz(k)
            qpfall(k)=qpfall(k)-(fz(kc)-fz(k))*irhoadz(k)*flagstat  ! For qp budget
            lat_heat = -(lfac(kc)*fz(kc)-lfac(k)*fz(k))*irhoadz(k)
            t(i,j,k)=t(i,j,k)-lat_heat
            tlat(k)=tlat(k)-lat_heat            ! For energy budget
            precflux(k) = precflux(k) - fz(k)*flagstat   ! For statistics
         end do
         precsfc(i,j) = precsfc(i,j) - fz(1)*flagstat ! For statistics
         precssfc(i,j) = precssfc(i,j) - fz(1)*(1.-omega(i,j,1))*flagstat ! For statistics
         prec_xy(i,j) = prec_xy(i,j) - fz(1)*flagstat ! For 2D output

         if (iprec.lt.nprec) then

            ! Re-compute precipitation velocity using new value of qp.
            do k=1,nzm
                  wp(k) = rhofac(k)*term_vel(i,j,k,ind)
                  ! Decrease precipitation velocity by factor of nprec
                  wp(k) = -wp(k)*rhow(k)*dtn/dz/float(nprec)
                  ! Note: Don't bother checking CFL condition at each
                  ! substep since it's unlikely that the CFL will
                  ! increase very much between substeps when using
                  ! monotonic advection schemes.
            end do

            fz(nz)=0.
            www(nz)=0.
            lfac(nz)=0.

         end if

      end do !iprec

  end do
end do	
	 

!call t_stopf ('precip_fall')

end subroutine precip_fall


