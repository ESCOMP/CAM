	subroutine crmsurface(bflx)
	
	
        use crmx_vars
        use crmx_params

	implicit none
	
	real, intent (in) :: bflx
	real u_h0, tau00, tauxm, tauym
	real diag_ustar
	integer i,j

!--------------------------------------------------------


        if(SFC_FLX_FXD.and..not.SFC_TAU_FXD) then

          uhl = uhl + dtn*utend(1)
          vhl = vhl + dtn*vtend(1)

	  tauxm = 0.
	  tauym = 0.

          do j=1,ny
           do i=1,nx
             u_h0 = max(1.,sqrt((0.5*(u(i+1,j,1)+u(i,j,1))+ug)**2+ &
                                   (0.5*(v(i,j+YES3D,1)+v(i,j,1))+vg)**2))
             tau00 = rho(1) * diag_ustar(z(1),bflx,u_h0,z0)**2 
             fluxbu(i,j) = -(0.5*(u(i+1,j,1)+u(i,j,1))+ug-uhl)/u_h0*tau00
             fluxbv(i,j) = -(0.5*(v(i,j+YES3D,1)+v(i,j,1))+vg-vhl)/u_h0*tau00
             tauxm = tauxm + fluxbu(i,j)
             tauym = tauym + fluxbv(i,j)
           end do
          end do

	  taux0 = taux0 + tauxm/dble(nx*ny)
	  tauy0 = tauy0 + tauym/dble(nx*ny)

        end if ! SFC_FLX_FXD

        return
        end





! ----------------------------------------------------------------------
!
! DISCLAIMER : this code appears to be correct but has not been
!              very thouroughly tested. If you do notice any
!              anomalous behaviour then please contact Andy and/or
!              Bjorn
!
! Function diag_ustar:  returns value of ustar using the below 
! similarity functions and a specified buoyancy flux (bflx) given in
! kinematic units
!
! phi_m (zeta > 0) =  (1 + am * zeta)
! phi_m (zeta < 0) =  (1 - bm * zeta)^(-1/4)
!
! where zeta = z/lmo and lmo = (theta_rev/g*vonk) * (ustar^2/tstar)
!
! Ref: Businger, 1973, Turbulent Transfer in the Atmospheric Surface 
! Layer, in Workshop on Micormeteorology, pages 67-100.
!
! Code writen March, 1999 by Bjorn Stevens
!
! Code corrected 8th June 1999 (obukhov length was wrong way up,
! so now used as reciprocal of obukhov length)

      real function diag_ustar(z,bflx,wnd,z0)

      implicit none
      real, parameter      :: vonk =  0.4   ! von Karmans constant
      real, parameter      :: g    = 9.81   ! gravitational acceleration
      real, parameter      :: am   =  4.8   !   "          "         "
      real, parameter      :: bm   = 19.3   !   "          "         "
      real, parameter      :: eps  = 1.e-10 ! non-zero, small number

      real, intent (in)    :: z             ! height where u locates
      real, intent (in)    :: bflx          ! surface buoyancy flux (m^2/s^3)
      real, intent (in)    :: wnd           ! wind speed at z
      real, intent (in)    :: z0            ! momentum roughness height

      integer :: iterate
      real    :: lnz, klnz, c1, x, psi1, zeta, rlmo, ustar

      lnz   = log(z/z0) 
      klnz  = vonk/lnz              
      c1    = 3.14159/2. - 3.*log(2.)

      ustar =  wnd*klnz
      if (bflx /= 0.0) then 
        do iterate=1,8
          rlmo   = -bflx * vonk/(ustar**3 + eps)   !reciprocal of
                                                   !obukhov length
          zeta  = min(1.,z*rlmo)
          if (zeta > 0.) then
            ustar =  vonk*wnd  /(lnz + am*zeta)
          else
            x     = sqrt( sqrt( 1.0 - bm*zeta ) )
            psi1  = 2.*log(1.0+x) + log(1.0+x*x) - 2.*atan(x) + c1
            ustar = wnd*vonk/(lnz - psi1)
          end if
        end do
      end if

      diag_ustar = ustar

      return
      end function diag_ustar
! ----------------------------------------------------------------------



      real function z0_est(z,bflx,wnd,ustar)

!
! Compute z0 from buoyancy flux, wind, and friction velocity
!
! 2004, Marat Khairoutdinov
!

      implicit none
      real, parameter      :: vonk =  0.4   ! von Karmans constant
      real, parameter      :: g    = 9.81   ! gravitational acceleration
      real, parameter      :: am   =  4.8   !   "          "         "
      real, parameter      :: bm   = 19.3   !   "          "         "
      real, parameter      :: eps  = 1.e-10 ! non-zero, small number

      real, intent (in)    :: z             ! height where u locates
      real, intent (in)    :: bflx          ! surface buoyancy flux (m^2/s^3)
      real, intent (in)    :: wnd           ! wind speed at z
      real, intent (in)    :: ustar         ! friction velocity

      real    :: lnz, klnz, c1, x, psi1, zeta, rlmo

      c1    = 3.14159/2. - 3.*log(2.)
      rlmo   = -bflx*vonk/(ustar**3+eps)   !reciprocal of
      zeta   = min(1.,z*rlmo)
      if (zeta >= 0.) then
            psi1 = -am*zeta
      else
            x     = sqrt( sqrt( 1.0 - bm*zeta ) )
            psi1  = 2.*log(1.0+x) + log(1.0+x*x) - 2.*atan(x) + c1
      end if
      lnz = max(0.,vonk*wnd/(ustar + eps) + psi1)
      z0_est = z*exp(-lnz)

      return
      end function z0_est
! ----------------------------------------------------------------------

