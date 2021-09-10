module wv_sat_methods

! This portable module contains all CAM methods for estimating
! the saturation vapor pressure of water.
!
! wv_saturation provides CAM-specific interfaces and utilities
! based on these formulae.
!
! Typical usage of this module:
!
! Init:
! call wv_sat_methods_init(r8, <constants>, errstring)
!
! Get scheme index from a name string:
! scheme_idx = wv_sat_get_scheme_idx(scheme_name)
! if (.not. wv_sat_valid_idx(scheme_idx)) <throw some error>
!
! Get pressures:
! es = wv_sat_svp_water(t, scheme_idx)
! es = wv_sat_svp_ice(t, scheme_idx)
!
! Use ice/water transition range:
! es = wv_sat_svp_trice(t, ttrice, scheme_idx)
!
! Note that elemental functions cannot be pointed to, nor passed
! as arguments. If you need to do either, it is recommended to
! wrap the function so that it can be given an explicit (non-
! elemental) interface.

implicit none
private
save

integer, parameter :: r8 = selected_real_kind(12) ! 8 byte real

integer, parameter :: VLENS = 128 ! vector length for a GPU kernel

real(r8) :: tmelt   ! Melting point of water at 1 atm (K)
real(r8) :: h2otrip ! Triple point temperature of water (K)
real(r8) :: tboil   ! Boiling point of water at 1 atm (K)

real(r8) :: ttrice  ! Ice-water transition range

real(r8) :: epsilo  ! Ice-water transition range
real(r8) :: omeps   ! 1._r8 - epsilo

! Indices representing individual schemes
integer, parameter :: Invalid_idx = -1
integer, parameter :: GoffGratch_idx = 1
integer, parameter :: MurphyKoop_idx = 2
integer, parameter :: Bolton_idx = 3

! Index representing the current default scheme.
integer, parameter :: initial_default_idx = GoffGratch_idx
integer :: default_idx = initial_default_idx

!$acc declare create (epsilo, tmelt, tboil, omeps, h2otrip, ttrice)

public wv_sat_methods_init
public wv_sat_get_scheme_idx
public wv_sat_valid_idx

public wv_sat_set_default
public wv_sat_reset_default

public wv_sat_qsat_water, wv_sat_qsat_water_vect
public wv_sat_qsat_ice, wv_sat_qsat_ice_vect

public wv_sat_svp_trans, wv_sat_svp_trans_vect

! pressure -> humidity conversion
public wv_sat_svp_to_qsat, wv_sat_svp_to_qsat_vect

! Combined qsat operations
public wv_sat_qsat_trans

public wv_sat_svp_water, wv_sat_svp_water_vect
public wv_sat_svp_ice, wv_sat_svp_ice_vect

contains

!---------------------------------------------------------------------
! ADMINISTRATIVE FUNCTIONS
!---------------------------------------------------------------------

! Get physical constants
subroutine wv_sat_methods_init(kind, tmelt_in, h2otrip_in, tboil_in, &
     ttrice_in, epsilo_in, errstring)
  integer, intent(in) :: kind
  real(r8), intent(in) :: tmelt_in
  real(r8), intent(in) :: h2otrip_in
  real(r8), intent(in) :: tboil_in
  real(r8), intent(in) :: ttrice_in
  real(r8), intent(in) :: epsilo_in
  character(len=*), intent(out)  :: errstring

  errstring = ' '

  if (kind /= r8) then
     write(errstring,*) 'wv_sat_methods_init: ERROR: ', &
          kind,' was input kind but ',r8,' is internal kind.'
     return
  end if

  if (ttrice_in < 0._r8) then
     write(errstring,*) 'wv_sat_methods_init: ERROR: ', &
          ttrice_in,' was input for ttrice, but negative range is invalid.'
     return
  end if

  tmelt = tmelt_in
  h2otrip = h2otrip_in
  tboil = tboil_in
  ttrice = ttrice_in
  epsilo = epsilo_in

  omeps = 1._r8 - epsilo

  !$acc update device (tmelt,h2otrip,tboil,ttrice,epsilo,omeps)

end subroutine wv_sat_methods_init

! Look up index by name.
pure function wv_sat_get_scheme_idx(name) result(idx)
  character(len=*), intent(in) :: name
  integer :: idx
  
  select case (name)
  case("GoffGratch")
     idx = GoffGratch_idx
  case("MurphyKoop")
     idx = MurphyKoop_idx
  case("Bolton")
     idx = Bolton_idx
  case default
     idx = Invalid_idx
  end select

end function wv_sat_get_scheme_idx

! Check validity of an index from the above routine.
pure function wv_sat_valid_idx(idx) result(status)
  integer, intent(in) :: idx
  logical :: status

  status = (idx /= Invalid_idx)

end function wv_sat_valid_idx

! Set default scheme (otherwise, Goff & Gratch is default)
! Returns a logical representing success (.true.) or
! failure (.false.).
function wv_sat_set_default(name) result(status)
  character(len=*), intent(in) :: name
  logical :: status

  ! Don't want to overwrite valid default with invalid,
  ! so assign to temporary and check it first.
  integer :: tmp_idx

  tmp_idx = wv_sat_get_scheme_idx(name)

  status = wv_sat_valid_idx(tmp_idx)

  if (status) default_idx = tmp_idx

end function wv_sat_set_default

! Reset default scheme to initial value.
! The same thing can be accomplished with wv_sat_set_default;
! the real reason to provide this routine is to reset the
! module for testing purposes.
subroutine wv_sat_reset_default()

  default_idx = initial_default_idx

end subroutine wv_sat_reset_default

!---------------------------------------------------------------------
! UTILITIES
!---------------------------------------------------------------------

! Get saturation specific humidity given pressure and SVP.
! Specific humidity is limited to range 0-1.
function wv_sat_svp_to_qsat(es, p) result(qs)
  real(r8), intent(in) :: es  ! SVP
  real(r8), intent(in) :: p   ! Current pressure.
  real(r8)             :: qs

  ! If pressure is less than SVP, set qs to maximum of 1.
  if ( (p - es) <= 0._r8 ) then
     qs = 1.0_r8
  else
     qs = epsilo*es / (p - omeps*es)
  end if

end function wv_sat_svp_to_qsat

! Get saturation specific humidity given pressure and SVP.
! Specific humidity is limited to range 0-1.
subroutine  wv_sat_svp_to_qsat_vect(es, p, qs, vlen)

  integer,  intent(in) :: vlen
  real(r8), intent(in)  :: es(vlen)  ! SVP
  real(r8), intent(in)  :: p(vlen)   ! Current pressure.
  real(r8), intent(out) :: qs(vlen)
  integer :: i

  ! If pressure is less than SVP, set qs to maximum of 1.

  !$acc data present (es,p,qs)

  !$acc parallel vector_length(VLENS) default(present)
  !$acc loop gang vector
  do i=1,vlen
     if ( (p(i) - es(i)) <= 0._r8 ) then
        qs(i) = 1.0_r8
     else
        qs(i) = epsilo*es(i) / (p(i) - omeps*es(i))
     end if
  end do
  !$acc end parallel

  !$acc end data
end subroutine wv_sat_svp_to_qsat_vect

subroutine wv_sat_qsat_water(t, p, es, qs, idx)
  !------------------------------------------------------------------!
  ! Purpose:                                                         !
  !   Calculate SVP over water at a given temperature, and then      !
  !   calculate and return saturation specific humidity.             !
  !------------------------------------------------------------------!

  ! Inputs
  real(r8), intent(in) :: t    ! Temperature
  real(r8), intent(in) :: p    ! Pressure
  ! Outputs
  real(r8), intent(out) :: es  ! Saturation vapor pressure
  real(r8), intent(out) :: qs  ! Saturation specific humidity

  integer,  intent(in), optional :: idx ! Scheme index

  es = wv_sat_svp_water(t, idx)

  qs = wv_sat_svp_to_qsat(es, p)

  ! Ensures returned es is consistent with limiters on qs.
  es = min(es, p)

end subroutine wv_sat_qsat_water

subroutine wv_sat_qsat_water_vect(t, p, es, qs, vlen, idx)
  !------------------------------------------------------------------!
  ! Purpose:                                                         !
  !   Calculate SVP over water at a given temperature, and then      !
  !   calculate and return saturation specific humidity.             !
  !------------------------------------------------------------------!
  ! Inputs

  integer,  intent(in) :: vlen
  real(r8), intent(in) :: t(vlen)    ! Temperature
  real(r8), intent(in) :: p(vlen)    ! Pressure
  ! Outputs
  real(r8), intent(out) :: es(vlen)  ! Saturation vapor pressure
  real(r8), intent(out) :: qs(vlen)  ! Saturation specific humidity

  integer,  intent(in), optional :: idx ! Scheme index
  integer :: i

  !$acc data present (t,p,es,qs)

  call wv_sat_svp_water_vect(t, es, vlen, idx)
  call wv_sat_svp_to_qsat_vect(es, p, qs, vlen)

  !$acc parallel vector_length(VLENS) default(present)
  !$acc loop gang vector
  do i=1,vlen
     ! Ensures returned es is consistent with limiters on qs.
     es(i) = min(es(i), p(i))
  enddo
  !$acc end parallel

  !$acc end data
end subroutine wv_sat_qsat_water_vect

subroutine wv_sat_qsat_ice(t, p, es, qs, idx)
  !------------------------------------------------------------------!
  ! Purpose:                                                         !
  !   Calculate SVP over ice at a given temperature, and then        !
  !   calculate and return saturation specific humidity.             !
  !------------------------------------------------------------------!

  ! Inputs
  real(r8), intent(in) :: t    ! Temperature
  real(r8), intent(in) :: p    ! Pressure
  ! Outputs
  real(r8), intent(out) :: es  ! Saturation vapor pressure
  real(r8), intent(out) :: qs  ! Saturation specific humidity

  integer,  intent(in), optional :: idx ! Scheme index

  es = wv_sat_svp_ice(t, idx)

  qs = wv_sat_svp_to_qsat(es, p)

  ! Ensures returned es is consistent with limiters on qs.
  es = min(es, p)

end subroutine wv_sat_qsat_ice

subroutine wv_sat_qsat_ice_vect(t, p, es, qs, vlen, idx)
  !------------------------------------------------------------------!
  ! Purpose:                                                         !
  !   Calculate SVP over ice at a given temperature, and then        !
  !   calculate and return saturation specific humidity.             !
  !------------------------------------------------------------------!
  ! Inputs

  integer,  intent(in) :: vlen
  real(r8), intent(in) :: t(vlen)    ! Temperature
  real(r8), intent(in) :: p(vlen)    ! Pressure
  ! Outputs
  real(r8), intent(out) :: es(vlen)  ! Saturation vapor pressure
  real(r8), intent(out) :: qs(vlen)  ! Saturation specific humidity

  integer,  intent(in), optional :: idx ! Scheme index
  integer :: i

  !$acc data present (t,p,es,qs)

  call wv_sat_svp_ice_vect(t, es, vlen, idx)
  call wv_sat_svp_to_qsat_vect(es, p, qs, vlen)

  !$acc parallel vector_length(VLENS) default(present)
  !$acc loop gang vector
  do i=1,vlen
     ! Ensures returned es is consistent with limiters on qs.
     es(i) = min(es(i), p(i))
  enddo
  !$acc end parallel

  !$acc end data
end subroutine wv_sat_qsat_ice_vect

subroutine wv_sat_qsat_trans(t, p, es, qs, idx)
  !------------------------------------------------------------------!
  ! Purpose:                                                         !
  !   Calculate SVP over ice at a given temperature, and then        !
  !   calculate and return saturation specific humidity.             !
  !------------------------------------------------------------------!

  ! Inputs
  real(r8), intent(in) :: t    ! Temperature
  real(r8), intent(in) :: p    ! Pressure
  ! Outputs
  real(r8), intent(out) :: es  ! Saturation vapor pressure
  real(r8), intent(out) :: qs  ! Saturation specific humidity

  integer,  intent(in), optional :: idx ! Scheme index

  es = wv_sat_svp_trans(t, idx)

  qs = wv_sat_svp_to_qsat(es, p)

  ! Ensures returned es is consistent with limiters on qs.
  es = min(es, p)

end subroutine wv_sat_qsat_trans

!---------------------------------------------------------------------
! SVP INTERFACE FUNCTIONS
!---------------------------------------------------------------------

function wv_sat_svp_water(t, idx) result(es)
  real(r8), intent(in) :: t
  integer,  intent(in), optional :: idx
  real(r8)             :: es

  integer :: use_idx

  if (present(idx)) then
     use_idx = idx
  else
     use_idx = default_idx
  end if

  select case (use_idx)
  case(GoffGratch_idx)
     es = GoffGratch_svp_water(t)
  case(MurphyKoop_idx)
     es = MurphyKoop_svp_water(t)
  case(Bolton_idx)
     es = Bolton_svp_water(t)
  end select

end function wv_sat_svp_water

subroutine  wv_sat_svp_water_vect(t, es, vlen, idx)
  integer,  intent(in) :: vlen
  real(r8), intent(in) :: t(vlen)
  integer,  intent(in), optional :: idx
  real(r8), intent(out) :: es(vlen)
  integer :: i
  integer :: use_idx

  !$acc data present (t,es)

  if (present(idx)) then
     use_idx = idx
  else
     use_idx = default_idx
  end if

  select case (use_idx)
  case(GoffGratch_idx)
     call GoffGratch_svp_water_vect(t,es,vlen)
  case(MurphyKoop_idx)
     call MurphyKoop_svp_water_vect(t,es,vlen)
  case(Bolton_idx)
     call Bolton_svp_water_vect(t,es,vlen)
  end select

  !$acc end data
end subroutine wv_sat_svp_water_vect

function wv_sat_svp_ice(t, idx) result(es)
  real(r8), intent(in)  :: t
  integer,  intent(in), optional :: idx
  real(r8)              :: es

  integer :: use_idx

  if (present(idx)) then
     use_idx = idx
  else
     use_idx = default_idx
  end if

  select case (use_idx)
  case(GoffGratch_idx)
     es = GoffGratch_svp_ice(t)
  case(MurphyKoop_idx)
     es = MurphyKoop_svp_ice(t)
  case(Bolton_idx)
     es = Bolton_svp_water(t)
  end select

end function wv_sat_svp_ice

subroutine wv_sat_svp_ice_vect(t, es, vlen, idx)
  integer,  intent(in) :: vlen
  real(r8), intent(in) :: t(vlen)
  integer,  intent(in), optional :: idx
  real(r8), intent(out) :: es(vlen)
  integer :: i

  integer :: use_idx

  !$acc data present (t,es)

  if (present(idx)) then
     use_idx = idx
  else
     use_idx = default_idx
  end if

  select case (use_idx)
  case(GoffGratch_idx)
     call GoffGratch_svp_ice_vect(t,es,vlen)
  case(MurphyKoop_idx)
     call MurphyKoop_svp_ice_vect(t,es,vlen)
  case(Bolton_idx)
     call Bolton_svp_water_vect(t,es,vlen)
  end select

  !$acc end data
end subroutine wv_sat_svp_ice_vect

function wv_sat_svp_trans(t, idx) result(es)

  real(r8), intent(in) :: t
  integer,  intent(in), optional :: idx
  real(r8) :: es

  real(r8) :: esice      ! Saturation vapor pressure over ice
  real(r8) :: weight     ! Intermediate scratch variable for es transition

!
! Water
!
  if (t >= (tmelt - ttrice)) then
     es = wv_sat_svp_water(t,idx)
  else
     es = 0.0_r8
  end if

!
! Ice
!
  if (t < tmelt) then

     esice = wv_sat_svp_ice(t,idx)

     if ( (tmelt - t) > ttrice ) then
        weight = 1.0_r8
     else
        weight = (tmelt - t)/ttrice
     end if

     es = weight*esice + (1.0_r8 - weight)*es
  end if

end function wv_sat_svp_trans

subroutine wv_sat_svp_trans_vect(t, es, vlen, idx)

  integer,  intent(in)  :: vlen
  real(r8), intent(in)  :: t(vlen)
  integer,  intent(in), optional :: idx
  real(r8), intent(out) :: es(vlen)

  real(r8) :: esice(vlen)      ! Saturation vapor pressure over ice
  real(r8) :: weight           ! Intermediate scratch variable for es transition
  integer  :: i

  !$acc data present (t,es) &
  !$acc      create  (esice)

!
! Water
!
  call wv_sat_svp_water_vect(t,es,vlen,idx)
  !$acc parallel vector_length(VLENS) default(present)
  !$acc loop gang vector
  do i = 1, vlen
     if (t(i) < (tmelt - ttrice)) then
        es(i) = 0.0_r8
     end if
  end do
  !$acc end parallel
!
! Ice
!
  call wv_sat_svp_ice_vect(t,esice,vlen,idx)
  !$acc parallel vector_length(VLENS) default(present)
  !$acc loop gang vector
  do i = 1, vlen
     if (t(i) < tmelt) then
        if ( (tmelt - t(i)) > ttrice ) then
           weight = 1.0_r8
        else
           weight = (tmelt - t(i))/ttrice
        end if
   
        es(i) = weight*esice(i) + (1.0_r8 - weight)*es(i)
     end if
  end do
  !$acc end parallel

  !$acc end data
end subroutine wv_sat_svp_trans_vect

!---------------------------------------------------------------------
! SVP METHODS
!---------------------------------------------------------------------

! Goff & Gratch (1946)

function GoffGratch_svp_water(t) result(es)
  real(r8), intent(in) :: t  ! Temperature in Kelvin
  real(r8) :: es             ! SVP in Pa

  ! uncertain below -70 C
  es = 10._r8**(-7.90298_r8*(tboil/t-1._r8)+ &
       5.02808_r8*log10(tboil/t)- &
       1.3816e-7_r8*(10._r8**(11.344_r8*(1._r8-t/tboil))-1._r8)+ &
       8.1328e-3_r8*(10._r8**(-3.49149_r8*(tboil/t-1._r8))-1._r8)+ &
       log10(1013.246_r8))*100._r8

end function GoffGratch_svp_water

subroutine GoffGratch_svp_water_vect(t, es, vlen)
  integer, intent(in) :: vlen
  real(r8), intent(in)  :: t(vlen)  ! Temperature in Kelvin
  real(r8), intent(out) :: es(vlen) ! SVP in Pa
  real(r8) :: log_tboil
  integer :: i

  !$acc data present (t,es)

  ! Goff, J. A., and S. Gratch. “Low-Pressure Properties of Water from -160F
  ! to 212F.” Trans. Am. Soc. Heat. Vent. Eng. 52 (1946): 95–121.
  ! uncertain below -70 C

  log_tboil = log10(tboil)

  !$acc parallel vector_length(VLENS) default(present)
  !$acc loop gang vector
  do i=1,vlen
     es(i) = 10._r8**(-7.90298_r8*(tboil/t(i)-1._r8)+ &
       5.02808_r8*(log_tboil-log10(t(i)))- &
       1.3816e-7_r8*(10._r8**(11.344_r8*(1._r8-t(i)/tboil))-1._r8)+ &
       8.1328e-3_r8*(10._r8**(-3.49149_r8*(tboil/t(i)-1._r8))-1._r8)+ &
       log10(1013.246_r8))*100._r8
  enddo
  !$acc end parallel

  !$acc end data
end subroutine GoffGratch_svp_water_vect

function GoffGratch_svp_ice(t) result(es)
  real(r8), intent(in) :: t  ! Temperature in Kelvin
  real(r8) :: es             ! SVP in Pa

  ! good down to -100 C
  es = 10._r8**(-9.09718_r8*(h2otrip/t-1._r8)-3.56654_r8* &
       log10(h2otrip/t)+0.876793_r8*(1._r8-t/h2otrip)+ &
       log10(6.1071_r8))*100._r8

end function GoffGratch_svp_ice

subroutine GoffGratch_svp_ice_vect(t, es, vlen)
  integer, intent(in)   :: vlen
  real(r8), intent(in)  :: t(vlen)  ! Temperature in Kelvin
  real(r8), intent(out) :: es(vlen) ! SVP in Pa
  real(r8), parameter :: log_param = log10(6.1071_r8)
  integer :: i
  ! good down to -100 C

  !$acc data present (t,es)

  !$acc parallel vector_length(VLENS) default(present)
  !$acc loop gang vector
  do i=1,vlen
     es(i) = 10._r8**(-9.09718_r8*(h2otrip/t(i)-1._r8)-3.56654_r8* &
          log10(h2otrip/t(i))+0.876793_r8*(1._r8-t(i)/h2otrip)+ &
          log_param)*100._r8
  enddo
  !$acc end parallel

  !$acc end data
end subroutine GoffGratch_svp_ice_vect

! Murphy & Koop (2005)

function MurphyKoop_svp_water(t) result(es)
  real(r8), intent(in) :: t  ! Temperature in Kelvin
  real(r8) :: es             ! SVP in Pa

  ! (good for 123 < T < 332 K)
  es = exp(54.842763_r8 - (6763.22_r8 / t) - (4.210_r8 * log(t)) + &
       (0.000367_r8 * t) + (tanh(0.0415_r8 * (t - 218.8_r8)) * &
       (53.878_r8 - (1331.22_r8 / t) - (9.44523_r8 * log(t)) + &
       0.014025_r8 * t)))

end function MurphyKoop_svp_water

subroutine MurphyKoop_svp_water_vect(t, es, vlen)
  integer, intent(in)   :: vlen
  real(r8), intent(in)  :: t(vlen)  ! Temperature in Kelvin
  real(r8), intent(out) :: es(vlen)             ! SVP in Pa
  
  integer :: i
  ! Murphy, D. M., and T. Koop. “Review of the Vapour Pressure of Ice and
  ! Supercooled Water for Atmospheric Applications.” Q. J. R. Meteorol.
  ! Soc. 131, no. 608 (2005): 1539–65. 10.1256/qj.04.94
  ! (good for 123 < T < 332 K)

  !$acc data present (t,es)

  !$acc parallel vector_length(VLENS) default(present)
  !$acc loop gang vector
  do i = 1, vlen
     es(i) = exp(54.842763_r8 - (6763.22_r8 / t(i)) - (4.210_r8 * log(t(i))) + &
          (0.000367_r8 * t(i)) + (tanh(0.0415_r8 * (t(i) - 218.8_r8)) * &
          (53.878_r8 - (1331.22_r8 / t(i)) - (9.44523_r8 * log(t(i))) + &
          0.014025_r8 * t(i))))
  end do
  !$acc end parallel

  !$acc end data
end subroutine MurphyKoop_svp_water_vect

function MurphyKoop_svp_ice(t) result(es)
  real(r8), intent(in) :: t  ! Temperature in Kelvin
  real(r8) :: es             ! SVP in Pa

  ! (good down to 110 K)
  es = exp(9.550426_r8 - (5723.265_r8 / t) + (3.53068_r8 * log(t)) &
       - (0.00728332_r8 * t))

end function MurphyKoop_svp_ice

subroutine MurphyKoop_svp_ice_vect(t, es, vlen)
  integer, intent(in)   :: vlen
  real(r8), intent(in) :: t(vlen)  ! Temperature in Kelvin
  real(r8), intent(out) :: es(vlen)             ! SVP in Pa
  
  integer :: i
  ! (good down to 110 K)

  !$acc data present (t,es)

  !$acc parallel vector_length(VLENS) default(present)
  !$acc loop gang vector
  do i = 1, vlen
     es(i) = exp(9.550426_r8 - (5723.265_r8 / t(i)) + (3.53068_r8 * log(t(i))) &
             - (0.00728332_r8 * t(i)))
  end do
  !$acc end parallel

  !$acc end data
end subroutine MurphyKoop_svp_ice_vect

! Bolton (1980)
! zm_conv deep convection scheme contained this SVP calculation.
! It appears to be from D. Bolton, 1980, Monthly Weather Review.
! Unlike the other schemes, no distinct ice formula is associated
! with it. (However, a Bolton ice formula exists in CLUBB.)

! The original formula used degrees C, but this function
! takes Kelvin and internally converts.

function Bolton_svp_water(t) result(es)
  real(r8),parameter :: c1 = 611.2_r8
  real(r8),parameter :: c2 = 17.67_r8
  real(r8),parameter :: c3 = 243.5_r8

  real(r8), intent(in) :: t  ! Temperature in Kelvin
  real(r8) :: es             ! SVP in Pa

  es = c1*exp( (c2*(t - tmelt))/((t - tmelt)+c3) )

end function Bolton_svp_water

subroutine Bolton_svp_water_vect(t, es,vlen)
  real(r8),parameter :: c1 = 611.2_r8
  real(r8),parameter :: c2 = 17.67_r8
  real(r8),parameter :: c3 = 243.5_r8

  integer, intent(in)   :: vlen
  real(r8), intent(in)  :: t(vlen)  ! Temperature in Kelvin
  real(r8), intent(out) :: es(vlen)             ! SVP in Pa

  integer :: i

  !$acc data present (t,es)

  !$acc parallel vector_length(VLENS) default(present)
  !$acc loop gang vector
  do i = 1, vlen
     es(i) = c1*exp( (c2*(t(i) - tmelt))/((t(i) - tmelt)+c3) )
  end do
  !$acc end parallel

  !$acc end data
end subroutine Bolton_svp_water_vect

end module wv_sat_methods
