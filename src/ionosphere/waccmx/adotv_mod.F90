module adotv_mod
  use shr_kind_mod,only: r8 => shr_kind_r8 ! 8-byte reals

  implicit none

contains

  subroutine calc_adotv(z, un, vn, wn, adotv1, adotv2, adota1, adota2, &
                        a1dta2, be3, sini, lev0, lev1, lon0, lon1, lat0, lat1)
     !
     ! Calculate adotv1,2, adota1,2, a1dta2 and be3.
     ! All fields should be on O+ grid
     !
     use edyn_params,  only: r0,h0
     use edyn_geogrid, only: jspole, jnpole
     use getapex,      only: &
          zb,                & ! downward component of magnetic field
          bmod,              & ! magnitude of magnetic field (gauss)
          dvec,              & ! (nlonp1,nlat,3,2)
          dddarr,            & ! (nlonp1,nlat)
          be3arr,            & ! (nlonp1,nlat)
          alatm                ! (nlonp1,0:nlatp1)
     !
     ! Args:
     integer,intent(in) :: lev0, lev1, lon0, lon1, lat0, lat1
     real(r8), dimension(lev0:lev1,lon0:lon1,lat0:lat1), intent(in) :: &
          z,    & ! geopotential height (cm)
          un,   & ! neutral zonal velocity (cm/s)
          vn      ! neutral meridional velocity (cm/s)
     real(r8), dimension(lev0:lev1,lon0:lon1,lat0:lat1), intent(in) :: &
          wn    ! vertical velocity (cm/s)

     real(r8), dimension(lon0:lon1,lat0:lat1,lev0:lev1), intent(out) :: &
          adotv1, adotv2
     real(r8), dimension(lon0:lon1,lat0:lat1), intent(out) :: &
          adota1, adota2, a1dta2, be3, sini
     !
     ! Local:
     integer  :: k, i, j
     real(r8) :: r0or, rat, sinalat
     real(r8) :: clm2(lon0:lon1,lat0:lat1)
     !
     adotv1 = 0.0_r8
     adotv2 = 0.0_r8
     adota1 = 0.0_r8
     adota2 = 0.0_r8
     a1dta2 = 0.0_r8
     be3 = 0.0_r8
     sini = 0.0_r8

     do j = lat0, lat1
        do i = lon0, lon1
           sinalat = sin(alatm(i,j))               ! sin(lam)
           clm2(i,j) = 1._r8 - (sinalat * sinalat) ! cos^2(lam)
           be3(i,j) = 1.e-9_r8*be3arr(i,j)         ! be3 is in T (be3arr in nT)
           sini(i,j) = zb(i,j)/bmod(i,j) ! sin(I_m)

           do k=lev0,lev1-1
              !
              ! d_1 = (R_0/R)^1.5
              r0or = r0/(r0 + 0.5_r8* (z(k,i,j) + z(k+1,i,j)) - h0)
              rat = 1.e-2_r8*r0or**1.5_r8 ! 1/100 conversion in cm
              !
              ! A_1 dot V = fac( d_1(1) u + d_1(2) v + d_1(3) w
              adotv1(i,j,k) = rat*(            &
                   dvec(i,j,1,1) * un(k,i,j) + &
                   dvec(i,j,2,1) * vn(k,i,j) + &
                   dvec(i,j,3,1) * wn(k,i,j))

              !
              ! Note: clm2 is being used here to represent the squared cosine
              !       of the quasi-dipole latitude, not of the M(90) latitude,
              !       since the wind values are aligned vertically,
              !       not along the field line.
              !
              rat = rat * sqrt((4._r8 - (3._r8 * clm2(i,j))) / &
                   (4._r8 - (3._r8 * r0or * clm2(i,j))))
              !
              ! A_2 dot V = fac( d_2(1) u + d_2(2) v + d_2(3) w
              adotv2(i,j,k) = rat * (          &
                   dvec(i,j,1,2) * un(k,i,j) + &
                   dvec(i,j,2,2) * vn(k,i,j) + &
                   dvec(i,j,3,2) * wn(k,i,j))
           end do ! k=lev0,lev1-1
        end do
     end do

     do j = lat0, lat1
        if (j==jspole .or. j==jnpole) cycle
        do i = lon0, lon1
           !
           ! Calculation of adota(n) = d(n)**2/D
           !   a1dta2   = (d(1) dot d(2)) /D
           !
           adota1(i,j) = (dvec(i,j,1,1)**2 + dvec(i,j,2,1)**2 + &
                dvec(i,j,3,1)**2) / dddarr(i,j)
           adota2(i,j) = (dvec(i,j,1,2)**2 + dvec(i,j,2,2)**2 + &
                dvec(i,j,3,2)**2) / dddarr(i,j)
           a1dta2(i,j) = (dvec(i,j,1,1) * dvec(i,j,1,2) +       &
                dvec(i,j,2,1) * dvec(i,j,2,2) +                 &
                dvec(i,j,3,1) * dvec(i,j,3,2)) / dddarr(i,j)
        end do ! i=lon0,lon1

     end do ! j=lat0,lat1

  end subroutine calc_adotv


end module adotv_mod
