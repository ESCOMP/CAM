!----------------------------------------------------------------------------
! For calculating background ionization due to star light, geo-corona radiation
! Applicable to high altitudes of WACCM and WACCMX 
! Module created by Francis Vitt 14 Feb 2013
! Background ionization algorithm provided by Stan Solomn
!----------------------------------------------------------------------------
module photo_bkgrnd

  use shr_kind_mod,   only : r8 => shr_kind_r8
  use mo_chem_utls,   only : get_rxt_ndx
  use ppgrid,         only : pver

  implicit none

  private
  public :: photo_bkgrnd_calc
  public :: photo_bkgrnd_init

  integer :: jo_ndx, jo2_ndx, jn2_ndx, jn_ndx, jno_ndx

contains

  !----------------------------------------------------------------------------
  ! look up corresponding reaction rate indices
  !----------------------------------------------------------------------------
  subroutine photo_bkgrnd_init()
    jo_ndx  = get_rxt_ndx( 'jeuv_1' ) ! O + hv  -> Op + e
    jo2_ndx = get_rxt_ndx( 'jeuv_5' ) ! O2 + hv -> O2p + e 
    jn2_ndx = get_rxt_ndx( 'jeuv_6' ) ! N2 + hv -> N2p + e 
    jn_ndx  = get_rxt_ndx( 'jeuv_10') ! N2 + hv -> N + Np + e
    jno_ndx = get_rxt_ndx( 'jno_i' )  ! NO + hv -> NOp + e
  end subroutine photo_bkgrnd_init

  !----------------------------------------------------------------------------
  ! Adds background ionization rates to WACCM's photolysis rates
  !----------------------------------------------------------------------------
  subroutine photo_bkgrnd_calc(f107, o_den, o2_den, n2_den, no_den, zint, rates, & 
                               qbko1_out, qbko2_out, qbkn2_out, qbkn1_out, qbkno_out )

   ! arguments
    real(r8), intent(in) :: f107
    real(r8), intent(in) :: o_den(:)           ! N density (molecules/cm^3)
    real(r8), intent(in) :: o2_den(:)          ! O2 density (molecules/cm^3)
    real(r8), intent(in) :: n2_den(:)          ! N2 density (molecules/cm^3)
    real(r8), intent(in) :: no_den(:)          ! NO density (molecules/cm^3)
    real(r8), intent(in) :: zint(:)            ! interface height (km)

    real(r8), intent(inout) :: rates(:,:)      ! photolysis rates (sec-1)

    real(r8), intent(out), optional :: qbko1_out(:) ! rate (cm-3 sec-1) of O  + hv ->  Op  + e 
    real(r8), intent(out), optional :: qbko2_out(:) ! rate (cm-3 sec-1) of O2 + hv ->  O2p + e 
    real(r8), intent(out), optional :: qbkn2_out(:) ! rate (cm-3 sec-1) of N2 + hv ->  N2p + e 
    real(r8), intent(out), optional :: qbkn1_out(:) ! rate (cm-3 sec-1) of N  + hv ->  Np  + e   
    real(r8), intent(out), optional :: qbkno_out(:) ! rate (cm-3 sec-1) of N2 + hv ->  NOp + e 

    ! local vars
    real(r8), parameter :: km2cm = 1.e5_r8
    integer, parameter :: nmaj = 3

    real(r8) :: zmaj(nmaj,pver)
    real(r8) :: zno(pver)
    real(r8) :: zvcd(nmaj,pver)
    real(r8) :: delz(pver)
    real(r8) :: qbko1(pver)
    real(r8) :: qbko2(pver)
    real(r8) :: qbkn2(pver)
    real(r8) :: qbkn1(pver)
    real(r8) :: qbkno(pver)

    integer :: k

    zmaj(1,:) = o_den(:) 
    zmaj(2,:) = o2_den(:) 
    zmaj(3,:) = n2_den(:) 
    zno(:)    = no_den(:)

    ! thickness of each layer (cm)
    delz(1:pver-1) = km2cm*(zint(1:pver-1) - zint(2:pver))
    delz(pver) = delz(pver-1)

    zvcd(:,:) = 0._r8

    ! intergate column above each layer
    do k = 2,pver
       zvcd(1,k) = zvcd(1,k-1) + delz(k) * o_den(k)
       zvcd(2,k) = zvcd(2,k-1) + delz(k) * o2_den(k)
       zvcd(3,k) = zvcd(3,k-1) + delz(k) * n2_den(k)
    enddo

    ! invoke Stan's background ionization method -- returns rates (cm-3 sec-1)
    call qback (zmaj,zno,zvcd,f107,nmaj,pver,qbko1,qbko2,qbkn2,qbkn1,qbkno)

    ! divide by densities to get photolysis rates (sec-1)
    if (jo_ndx>0)  rates(:,jo_ndx)  = rates(:,jo_ndx)  + qbko1(:)/o_den(:)
    if (jo2_ndx>0) rates(:,jo2_ndx) = rates(:,jo2_ndx) + qbko2(:)/o2_den(:)
    if (jn2_ndx>0) rates(:,jn2_ndx) = rates(:,jn2_ndx) + qbkn2(:)/n2_den(:)
    if (jn_ndx >0) rates(:,jn_ndx)  = rates(:,jn_ndx)  + qbkn1(:)/n2_den(:)
    if (jno_ndx>0) rates(:,jno_ndx) = rates(:,jno_ndx) + qbkno(:)/no_den(:)

    ! for diagnostics
    if (present(qbko1_out)) qbko1_out(:) = qbko1(:)
    if (present(qbko2_out)) qbko2_out(:) = qbko2(:)
    if (present(qbkn2_out)) qbkn2_out(:) = qbkn2(:)
    if (present(qbkn1_out)) qbkn1_out(:) = qbkn1(:)
    if (present(qbkno_out)) qbkno_out(:) = qbkno(:)

  endsubroutine photo_bkgrnd_calc

!----------------------------------------------------------------------------
! Private Method
!----------------------------------------------------------------------------
!
! Stan Solomon, 11/88, 11/92
! Comment updated 3/05
! New version uses updated TIE-GCM and TIME-GCM qinite.F formulation, 1/13
!
! Estimates background ("nighttime") ionization rates.
! Four components are used:
!    Geocoronal Lyman-beta 102.6 nm (ionizes O2 only)
!    Geocoronal He I 58.4 nm
!    Geocoronal He II 30.4 nm
!    Geocoronal Lyman-alpha 121.6 nm (ionizes NO only)
!
! Definitions:
!
! zmaj    major species O, O2, N2 at each altitude
! zno     nitric oxide at each altitude
! zvcd    vertical column density for each major species above each altitude
! photoi  photoionization rates for each state, species, altitude
! f107    solar 10.7 cm radio flux activity index
! jm      number of altitude levels
! nmaj    number of major species (3)
! nst     number of states
!
! al      photon flux at 102.6 nm, 58.4 nm, 30.4 nm
! flyan   photon flux at 121.6 nm
! sa      absorption cross sections for O, O2, N2 at each wavelength
! si      ionization cross sections for O, O2, N2 at each wavelength
! flyan   photon flux at 121.6 nm
! salyao2 absorption cross section for O2 at 121.6 nm
! silyano ionization cross section for NO at 121.6 nm
! bn2p    branching ratio for N2+ from ionization of N2
! bn1p    branching ratio for N+ from ionization of N2
! tau     optical depth
! qbko1    production rate of O+
! qbko2    production rate of O2+
! qbkn2    production rate of N2+
! qbkn1    production rate of N+
! qbkno    production rate of NO+
!
! All units cgs.
!
!
subroutine qback (zmaj,zno,zvcd,f107,nmaj,jm,qbko1,qbko2,qbkn2,qbkn1,qbkno)

 ! args:
  integer,  intent(in)  :: nmaj,jm
  real(r8), intent(in)  :: f107
  real(r8), intent(in)  :: zmaj(nmaj,jm), zno(jm), zvcd(nmaj,jm)
  real(r8), intent(out) :: qbko1(jm),qbko2(jm),qbkn2(jm),qbkn1(jm),qbkno(jm)

 ! local vars:
  real(r8) :: al(3), sa(3,3), si(3,3)
  real(r8) :: salyao2, silyano, bn2p, bn1p
  real(r8) :: flyan
  real(r8) :: tau
  integer :: j,l

  data al /1.5e7_r8, 1.5e6_r8, 1.5e6_r8/

  data sa /      0._r8,  1.6e-18_r8,       0._r8, &
           10.2e-18_r8, 22.0e-18_r8, 23.1e-18_r8, &
            8.4e-18_r8, 16.0e-18_r8, 11.6e-18_r8/

  data si /      0._r8,  1.0e-18_r8,       0._r8, &
           10.2e-18_r8, 22.0e-18_r8, 23.1e-18_r8, &
            8.4e-18_r8, 16.0e-18_r8, 11.6e-18_r8/

  data salyao2/8.0e-21_r8/
  data silyano/2.0e-18_r8/
  data bn2p/0.86_r8/
  data bn1p/0.14_r8/

!
! Calculate Lyman-alpha 121.6 nm geocoronal flux as a function of F10.7:
!
  flyan = 5.E9_r8*(1._r8+0.002_r8*(f107-65._r8))
!
! Loop over altitudes:
!
  do j=1,jm
!
! Calculate optical depth and ionization rates for major species:
!
     qbko1(j)=0._r8
     qbko2(j)=0._r8
     qbkn2(j)=0._r8
     qbkn1(j)=0._r8
     do l=1,3
        tau=(sa(1,l)*zvcd(1,j)+sa(2,l)*zvcd(2,j)+sa(3,l)*zvcd(3,j))
        qbko1(j) = qbko1(j) +       al(l)*si(1,l)*zmaj(1,j)*exp(-tau)
        qbko2(j) = qbko2(j) +       al(l)*si(2,l)*zmaj(2,j)*exp(-tau)
        qbkn2(j) = qbkn2(j) + bn2p*(al(l)*si(3,l)*zmaj(3,j)*exp(-tau))
        qbkn1(j) = qbkn1(j) + bn1p*(al(l)*si(3,l)*zmaj(3,j)*exp(-tau))
     enddo
!
! Calculate optical depth of Ly-alpha, and ionization rate of NO:
!
     tau = salyao2*zvcd(2,j)
     qbkno(j) = flyan*silyano*zno(j)*exp(-tau)

  enddo

  return
end subroutine qback

end module photo_bkgrnd
