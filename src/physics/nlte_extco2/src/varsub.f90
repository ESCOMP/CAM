!=======================================================================!
! MODULE VARSUB ---     module with some useful routines                !
!-----------------------------------------------------------------------!
! STATUS: bf 23-Feb-2024                      CREATED: mh 23-Feb-1997   !
!-----------------------------------------------------------------------!
! COPYRIGHT: (C) 2000-2023 Karlsruhe Institute of Technology (KIT)      ! 
!            (C) 2000-2023 Instituto de Astrofisica de Andalucia (CSIC) !
!-----------------------------------------------------------------------!
! LICENSE: GNU Lesser General Public License (LGPL) version 2.1         !
!          (see file LICENSE_LGPLv2.1)                                  !
!-----------------------------------------------------------------------!
! CONTAINS:                                                             !
! public:                                                               !
!             errunit, error, seconds                                   !
!-----------------------------------------------------------------------!
! MAINTENANCE HISTORY:                                                  !
!=======================================================================!

module varsub

integer, parameter :: errunit = 6 

private
public :: errunit, error, seconds
contains

!___________________________________________________________________________
!                             H E A D E R
!___________________________________________________________________________
!SUBROUTINE:                  error
!Created by:                  Michael Hoepfner
!Date of creation:            12-1-99
!___________________________________________________________________________
!MAINTENANCE HISTORY:
!          grabow  30-Sep-99: output to logunit only if different from errunit
!___________________________________________________________________________
!DESCRIPTION:                 stop program
!___________________________________________________________________________
!VARIABLES:
!
!CHARACTER:
!
!argument list:
! chardum                     in change integer in character
! mess1                       in first error message
! mess2                       in second error message
! mess3                       in third error message
! mess4                       in fourth error message
! mess5                       in fifth error message
! place                       in subroutine where error occurs
!___________________________________________________________________________
!FUNCTIONS AND SUBROUTINES:
!___________________________________________________________________________
!___________________________________________________________________________
subroutine error(place,mess1,mess2,mess3,mess4,mess5,mess6,mess7)

  character(*),intent(in)             :: place, mess1
  character(*),intent(in),optional    :: mess2, mess3, mess4, mess5, mess6, mess7

  write(errunit,*)
  write(errunit,*) '>>>>>>>>              STOP              <<<<<<<<'
  write(errunit,*)
  write(errunit,*) 'Error occurred in subroutine: '//place
  write(errunit,*) mess1
  if (present(mess2)) write(errunit,*) mess2
  if (present(mess3)) write(errunit,*) mess3
  if (present(mess4)) write(errunit,*) mess4
  if (present(mess5)) write(errunit,*) mess5
  if (present(mess6)) write(errunit,*) mess6
  if (present(mess7)) write(errunit,*) mess7
  write(errunit,*) 'Stop program!'
  write(errunit,*)
  write(errunit,*) '>>>>>>>>              STOP              <<<<<<<<'

  close(errunit)
  stop
end subroutine error


!___________________________________________________________________________
!                             H E A D E R
!___________________________________________________________________________
!FUNCTION:                    seconds
!Created by:                  Michael Hoepfner
!Date of creation:            23-7-98
!Date of last modification:
!Modified by:
!Modification:
!___________________________________________________________________________
!DESCRIPTION:                 @@@!!!@@@
!___________________________________________________________________________
!VARIABLES:
!
!INTEGER:
!
!argument list:
! i1                          in @@@!!!@@@
! i2                          in
! ifirst                      in
! ilast                       in
! istart                      in
!
!COMPLEX:
!
!argument list:
! seconds                     in
!___________________________________________________________________________
!FUNCTIONS AND SUBROUTINES:
!                             call system_clock
!___________________________________________________________________________
!___________________________________________________________________________
function seconds()
  integer,save :: i1,i2,istart,ifirst,ilast
  complex :: seconds
  data ifirst /0/
  if (ifirst == 0) then
     call system_clock(istart,i2)
     ilast = istart
     seconds = (0.,0.)
     ifirst = 1
  else
     call system_clock(i1,i2)
     seconds = cmplx((i1-ilast)/real(i2),(i1-istart)/real(i2))
     ilast = i1
  endif
end function seconds


end module varsub
