module ionos_state_mod
  use shr_kind_mod, only : r8 => shr_kind_r8
  use ppgrid,       only : pcols, pver, pverp
  use mo_jeuv,      only : nIonRates ! Number of ionization rates in mo_photo

  implicit none

  type ionos_state

     real(r8), dimension(pcols)      :: cosZenAngR ! cosine of zenith angle (radians)
     real(r8), dimension(pcols)      :: zenAngD    ! zenith angle (degrees)

     real(r8), dimension(pcols,pver) :: bNorth3d   ! northward component of magnetic field units?
     real(r8), dimension(pcols,pver) :: bEast3d    ! eastward component of magnetic field
     real(r8), dimension(pcols,pver) :: bDown3d    ! downward component of magnetic field

     real(r8), dimension(pcols,pver,nIonRates) :: ionPRates    ! ionization rates temporary array (s-1 cm-3)
     real(r8), dimension(pcols,pver)           :: sumIonPRates ! Sum of ionization rates for O+,O2+,N+,N2+,NO+ (s-2 cm-3)
     ! ?? should be (s-1 cm-3) ??

     real(r8), dimension(pcols,pver)  :: dipMag   ! dip angle for each column (radians)

     real(r8), dimension(pcols,pverp) :: tNInt    ! Interface Temperature (K)

     real(r8), dimension(pcols,pver)  :: ndensN2  ! N2 number density (cm-3)
     real(r8), dimension(pcols,pver)  :: ndensO2  ! O2 number density (cm-3)
     real(r8), dimension(pcols,pver)  :: ndensO1  ! O number density (cm-3)
     real(r8), dimension(pcols,pver)  :: ndensNO  ! NO number density (cm-3)
     real(r8), dimension(pcols,pver)  :: ndensN1  ! N number density  (cm-3)
     real(r8), dimension(pcols,pver)  :: ndensE   ! E electron number density (cm-3)
     real(r8), dimension(pcols,pver)  :: ndensOp  ! O plus number density (cm-3)
     real(r8), dimension(pcols,pver)  :: ndensO2p ! O2 plus ion number density (cm-3)
     real(r8), dimension(pcols,pver)  :: ndensNOp ! NO plus ion number density  (cm-3)

     real(r8), dimension(pcols,pver)  :: sourceg4 ! g4 source term for electron/ion temperature update

     real(r8), dimension(pcols,pverp) :: rairvi   ! Constituent dependent gas constant on interface levels

     real(r8) :: n2_mmr(pcols,pver)

  end type ionos_state

end module ionos_state_mod
