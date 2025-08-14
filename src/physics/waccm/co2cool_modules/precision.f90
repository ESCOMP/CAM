!==========================================================================!
! MODULE PRECISION --- Defines the system definitions of representation of !
!                      numbers so that desired accuracy can be determined  !
!                      automatically by the chosen computer system         !
!--------------------------------------------------------------------------!
! STATUS : gra 25-Oct-05                            CREATED: afg 20-Aug-96 !
!--------------------------------------------------------------------------!
! COPYRIGHT: (C) 2020-2021 Karlsruhe Institute of Technology (KIT)         !
!            (C) 1996-2021 Dr.Udo Grabowski                                !
!--------------------------------------------------------------------------!
! LICENSE: GNU Lesser General Public License (LGPL) version 2.1            !
!          (see file LICENSE_LGPLv2.1)                                     !
!--------------------------------------------------------------------------!
! ACKNOWLEDGEMENTS:                                                        !
!          Parts of this code have been developed during the 'DarkStar'    !
!          Project at the Astronomical Institute Basel,supported by the    !
!          Swiss National Science Foundation (SNSF) project 2100-045568    !
!--------------------------------------------------------------------------!
! MAINTENANCE HISTORY:                                                     !
! gra 25-Oct-05: bytes_xx moved to os_specific, no standard way possible   !
! gra 26-Feb-99: added db=dp precision parameter                           !
! gra 08-Feb-99: precision now specified with selected_real_kind           !
! gra 11-Sep-96: changed name to precision, dg to dp, added sp             !
! afg 20-Aug-96: First specified                                           !
!==========================================================================!
module PRECISION

  use shr_kind_mod, only: SHR_KIND_R4, SHR_KIND_R8

  implicit NONE
!
! --> single and double precision kind parameters
!
!  integer, parameter :: bytes_sp = 4
!  integer, parameter :: bytes_dp = 8
!  integer, parameter :: bytes_qp = 16
!  integer, parameter :: bytes_db = bytes_dp, bytes_sg = bytes_sp
  integer, parameter :: sp = SHR_KIND_R4 ! bytes_sp
  integer, parameter :: dp = SHR_KIND_R8 ! bytes_dp
!  integer, parameter :: qp = bytes_qp
  integer, parameter :: db = dp
!  integer, parameter :: sg = sp
!
! --> shortcuts for double precision small numbers
!
  real (KIND=dp), parameter :: zero  = 0.0_dp, one  = 1.0_dp, two  = 2.0_dp
  real (KIND=dp), parameter :: three = 3.0_dp, four = 4.0_dp, five = 5.0_dp

private
!public :: sp,sg,dp,db,qp, bytes_sp, bytes_dp, bytes_qp, bytes_sg, bytes_db, zero,one,two,three,four,five

public :: sp, dp, db
public :: zero, one, two, three, four, five

end module PRECISION
