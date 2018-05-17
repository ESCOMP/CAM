      module mo_lin_matrix
      private
      public :: linmat
      contains
      subroutine linmat01( mat, y, rxt, het_rates )
!----------------------------------------------
! ... linear matrix entries for implicit species
!----------------------------------------------
      use chem_mods, only : gas_pcnst, rxntot, nzcnt
      use shr_kind_mod, only : r8 => shr_kind_r8
      implicit none
!----------------------------------------------
! ... dummy arguments
!----------------------------------------------
      real(r8), intent(in) :: y(gas_pcnst)
      real(r8), intent(in) :: rxt(rxntot)
      real(r8), intent(in) :: het_rates(max(1,gas_pcnst))
      real(r8), intent(inout) :: mat(nzcnt)
         mat(423) = -( rxt(3) + rxt(4) + het_rates(1) )
         mat(735) = -( rxt(66) + rxt(67) + rxt(68) + rxt(79) + rxt(80) + rxt(81) &
                 + het_rates(2) )
         mat(602) = rxt(1) + 2.000_r8*rxt(2) + rxt(72) + rxt(73) + rxt(74) &
                      + 2.000_r8*rxt(77) + rxt(84) + rxt(85) + rxt(86) + 2.000_r8*rxt(89)
         mat(434) = rxt(4)
         mat(457) = rxt(6)
         mat(484) = rxt(8)
         mat(63) = rxt(10)
         mat(546) = rxt(12)
         mat(571) = rxt(21)
         mat(390) = rxt(24)
         mat(68) = rxt(25)
         mat(322) = rxt(32)
         mat(232) = rxt(62)
         mat(51) = rxt(63)
         mat(271) = rxt(65)
         mat(694) = rxt(105)
         mat(693) = -( rxt(105) + rxt(109)*y(4) + rxt(110)*y(4) + rxt(112)*y(36) &
                      + rxt(113)*y(37) + rxt(114)*y(38) + rxt(115)*y(46) + rxt(116)*y(47) &
                      + rxt(117)*y(39) + rxt(118)*y(44) + rxt(119)*y(45) + rxt(120)*y(40) &
                      + rxt(121)*y(35) + rxt(122)*y(43) + rxt(123)*y(42) + rxt(124)*y(48) &
                      + rxt(125)*y(49) + rxt(126)*y(50) + rxt(127)*y(51) + rxt(130)*y(12) &
                      + rxt(131)*y(12) + rxt(132)*y(12) + het_rates(97) )
         mat(601) = rxt(1)
         mat(433) = rxt(3)
         mat(570) = rxt(20)
         mat(597) = -( rxt(1) + rxt(2) + rxt(70) + rxt(72) + rxt(73) + rxt(74) + rxt(77) &
                      + rxt(82) + rxt(84) + rxt(85) + rxt(86) + rxt(89) + het_rates(3) )
         mat(429) = rxt(4)
         mat(541) = rxt(13)
         mat(34) = rxt(100)
         mat(31) = rxt(103) + rxt(104)
         mat(689) = rxt(110)*y(4)
         mat(33) = -( rxt(97) + rxt(100) + rxt(99)*y(56) + het_rates(95) )
         mat(30) = -( rxt(103) + rxt(104) + het_rates(96) )
         mat(412) = rxt(3)
         mat(32) = rxt(97) + rxt(99)*y(56)
         mat(327) = -( het_rates(17) )
         mat(285) = rxt(18)
         mat(556) = rxt(20)
         mat(679) = rxt(132)*y(12)
         mat(102) = -( het_rates(16) )
         mat(280) = rxt(17) + rxt(18)
         mat(70) = rxt(64)
         mat(341) = rxt(225)*y(34)
         mat(139) = rxt(288)*y(56)
         mat(209) = -( rxt(69) + het_rates(5) )
         mat(438) = rxt(6)
         mat(143) = rxt(285)
         mat(447) = -( rxt(6) + rxt(7) + het_rates(6) )
         mat(474) = rxt(8) + .500_r8*rxt(248)
         mat(60) = rxt(10)
         mat(536) = rxt(13)
         mat(162) = rxt(295)
         mat(684) = 2.000_r8*rxt(109)*y(4)
         mat(475) = -( rxt(8) + rxt(248) + het_rates(7) )
         mat(61) = rxt(9) + rxt(168)
         mat(250) = rxt(11)
         mat(537) = rxt(12)
         mat(91) = rxt(15) + rxt(177)
         mat(239) = rxt(30)
         mat(111) = rxt(36)
         mat(513) = -( rxt(226)*y(34) + rxt(227)*y(41) + rxt(228)*y(39) + rxt(229)*y(35) &
                      + rxt(231)*y(44) + rxt(232)*y(45) + rxt(233)*y(51) + rxt(234)*y(50) &
                      + rxt(237)*y(12) + het_rates(86) )
         mat(251) = rxt(11)
         mat(92) = rxt(14)
         mat(80) = rxt(16)
         mat(563) = rxt(19)
         mat(116) = 2.000_r8*rxt(22)
         mat(222) = rxt(27)
         mat(193) = rxt(33)
         mat(476) = .500_r8*rxt(248)
         mat(686) = rxt(130)*y(12)
         mat(539) = -( rxt(12) + rxt(13) + rxt(247) + het_rates(8) )
         mat(62) = rxt(9) + rxt(10) + rxt(168)
         mat(93) = rxt(14)
         mat(241) = rxt(29)
         mat(112) = rxt(35)
         mat(247) = -( rxt(11) + het_rates(9) )
         mat(59) = 2.000_r8*rxt(246) + 2.000_r8*rxt(267) + 2.000_r8*rxt(273) &
                      + 2.000_r8*rxt(278)
         mat(526) = rxt(247)
         mat(465) = .500_r8*rxt(248)
         mat(236) = rxt(268) + rxt(274) + rxt(279)
         mat(108) = rxt(269) + rxt(277) + rxt(280)
         mat(90) = -( rxt(14) + rxt(15) + rxt(177) + het_rates(10) )
         mat(58) = -( rxt(9) + rxt(10) + rxt(168) + rxt(246) + rxt(267) + rxt(273) &
                      + rxt(278) + het_rates(11) )
         mat(643) = -( het_rates(13) )
         mat(691) = rxt(130)*y(12)
         mat(359) = rxt(184)*y(12)
         mat(157) = rxt(223)*y(12)
         mat(518) = rxt(237)*y(12)
         mat(77) = -( rxt(16) + het_rates(14) )
         mat(284) = -( rxt(17) + rxt(18) + het_rates(15) )
         mat(79) = rxt(16)
         mat(678) = rxt(131)*y(12) + rxt(132)*y(12)
         mat(272) = -( het_rates(18) )
         mat(78) = rxt(16)
         mat(283) = 2.000_r8*rxt(17)
         mat(554) = rxt(19) + 2.000_r8*rxt(21)
         mat(653) = rxt(28)
         mat(198) = rxt(34)
         mat(46) = rxt(57)
         mat(677) = rxt(131)*y(12)
         mat(623) = -( rxt(249) + het_rates(87) )
         mat(96) = rxt(15) + rxt(177)
         mat(690) = rxt(131)*y(12)
         mat(358) = rxt(225)*y(34) + rxt(230)*y(35)
         mat(517) = rxt(226)*y(34) + rxt(229)*y(35)
         mat(114) = -( rxt(22) + het_rates(19) )
         mat(605) = .500_r8*rxt(249)
         mat(565) = -( rxt(19) + rxt(20) + rxt(21) + het_rates(98) )
         mat(29) = rxt(61)
         mat(515) = rxt(226)*y(34) + rxt(227)*y(41) + rxt(228)*y(39) + rxt(229)*y(35) &
                      + rxt(233)*y(51) + rxt(237)*y(12)
         mat(349) = -( rxt(184)*y(12) + rxt(225)*y(34) + rxt(230)*y(35) + rxt(235)*y(51) &
                      + rxt(236)*y(50) + het_rates(84) )
         mat(36) = 2.000_r8*rxt(23)
         mat(376) = rxt(24)
         mat(22) = 2.000_r8*rxt(26)
         mat(220) = rxt(27)
         mat(656) = rxt(28)
         mat(237) = rxt(29)
         mat(42) = rxt(31)
         mat(39) = rxt(56)
         mat(680) = 2.000_r8*rxt(112)*y(36) + 2.000_r8*rxt(113)*y(37) &
                      + 2.000_r8*rxt(114)*y(38) + 2.000_r8*rxt(115)*y(46) + rxt(116)*y(47) &
                      + rxt(117)*y(39) + rxt(118)*y(44) + rxt(119)*y(45) &
                      + 4.000_r8*rxt(120)*y(40) + rxt(122)*y(43)
         mat(507) = rxt(226)*y(34) + 3.000_r8*rxt(227)*y(41) + rxt(228)*y(39) &
                      + rxt(231)*y(44) + rxt(232)*y(45)
         mat(35) = -( rxt(23) + het_rates(22) )
         mat(377) = -( rxt(24) + het_rates(23) )
         mat(67) = rxt(25)
         mat(238) = rxt(30)
         mat(23) = 2.000_r8*rxt(196)
         mat(64) = -( rxt(25) + het_rates(24) )
         mat(21) = -( rxt(26) + rxt(196) + het_rates(25) )
         mat(668) = -( rxt(28) + het_rates(26) )
         mat(360) = rxt(184)*y(12) + 2.000_r8*rxt(225)*y(34) + rxt(230)*y(35) &
                      + rxt(235)*y(51) + rxt(236)*y(50)
         mat(219) = -( rxt(27) + het_rates(27) )
         mat(234) = rxt(268) + rxt(274) + rxt(279)
         mat(235) = -( rxt(29) + rxt(30) + rxt(268) + rxt(274) + rxt(279) + het_rates(28) &
       )
         mat(41) = -( rxt(31) + het_rates(29) )
         mat(399) = -( het_rates(85) )
         mat(43) = rxt(31)
         mat(311) = rxt(32)
         mat(192) = rxt(33)
         mat(199) = rxt(34)
         mat(110) = rxt(35)
         mat(682) = rxt(121)*y(35) + rxt(122)*y(43) + rxt(123)*y(42) &
                      + 2.000_r8*rxt(124)*y(48) + 2.000_r8*rxt(125)*y(49) &
                      + 3.000_r8*rxt(126)*y(50) + 2.000_r8*rxt(127)*y(51)
         mat(509) = rxt(229)*y(35) + 2.000_r8*rxt(233)*y(51) + 3.000_r8*rxt(234)*y(50)
         mat(351) = rxt(230)*y(35) + 2.000_r8*rxt(235)*y(51) + 3.000_r8*rxt(236)*y(50)
         mat(307) = -( rxt(32) + het_rates(30) )
         mat(109) = rxt(36)
         mat(197) = -( rxt(34) + het_rates(31) )
         mat(189) = -( rxt(33) + het_rates(32) )
         mat(107) = rxt(269) + rxt(277) + rxt(280)
         mat(106) = -( rxt(35) + rxt(36) + rxt(269) + rxt(277) + rxt(280) + het_rates(33) &
       )
         mat(122) = -( het_rates(88) )
         mat(158) = -( rxt(295) + het_rates(89) )
         mat(579) = rxt(70) + rxt(82)
         mat(141) = rxt(288)*y(56)
         mat(83) = -( het_rates(90) )
         mat(204) = rxt(69)
         mat(140) = -( rxt(285) + rxt(288)*y(56) + het_rates(91) )
         mat(705) = rxt(66) + rxt(67) + rxt(68) + rxt(79) + rxt(80) + rxt(81)
         mat(578) = rxt(72) + rxt(73) + rxt(74) + rxt(84) + rxt(85) + rxt(86)
         mat(167) = -( het_rates(92) )
         mat(436) = rxt(7)
         mat(142) = rxt(285)
         mat(159) = rxt(295)
         mat(97) = -( het_rates(94) )
         mat(179) = -( het_rates(93) )
         mat(437) = rxt(7)
         mat(708) = rxt(66) + rxt(67) + rxt(68) + rxt(79) + rxt(80) + rxt(81)
         mat(208) = rxt(69)
         mat(581) = rxt(70) + rxt(72) + rxt(73) + rxt(74) + rxt(82) + rxt(84) + rxt(85) &
                      + rxt(86)
         mat(24) = -( rxt(55) + het_rates(52) )
         mat(671) = rxt(113)*y(37) + rxt(114)*y(38) + 2.000_r8*rxt(115)*y(46) &
                      + 2.000_r8*rxt(116)*y(47) + rxt(117)*y(39) + rxt(119)*y(45) &
                      + rxt(122)*y(43) + rxt(123)*y(42) + rxt(124)*y(48) &
                      + 2.000_r8*rxt(125)*y(49)
         mat(485) = rxt(228)*y(39) + rxt(232)*y(45)
         mat(37) = -( rxt(56) + het_rates(53) )
         mat(673) = rxt(112)*y(36) + rxt(114)*y(38) + rxt(118)*y(44)
         mat(486) = rxt(231)*y(44)
         mat(44) = -( rxt(57) + het_rates(54) )
         mat(149) = rxt(223)*y(12)
         mat(150) = -( rxt(223)*y(12) + het_rates(55) )
         mat(25) = 2.000_r8*rxt(55)
         mat(38) = rxt(56)
         mat(45) = rxt(57)
         mat(674) = rxt(116)*y(47) + rxt(123)*y(42)
      end subroutine linmat01
      subroutine linmat02( mat, y, rxt, het_rates )
!----------------------------------------------
! ... linear matrix entries for implicit species
!----------------------------------------------
      use chem_mods, only : gas_pcnst, rxntot, nzcnt
      use shr_kind_mod, only : r8 => shr_kind_r8
      implicit none
!----------------------------------------------
! ... dummy arguments
!----------------------------------------------
      real(r8), intent(in) :: y(gas_pcnst)
      real(r8), intent(in) :: rxt(rxntot)
      real(r8), intent(in) :: het_rates(max(1,gas_pcnst))
      real(r8), intent(inout) :: mat(nzcnt)
         mat(69) = -( rxt(64) + het_rates(57) )
         mat(132) = -( het_rates(58) )
         mat(71) = rxt(64)
         mat(256) = rxt(65)
         mat(258) = -( rxt(65) + het_rates(59) )
         mat(228) = rxt(62)
         mat(227) = -( rxt(62) + het_rates(60) )
         mat(49) = rxt(63)
         mat(48) = -( rxt(63) + het_rates(61) )
         mat(28) = rxt(61)
         mat(27) = -( rxt(61) + het_rates(62) )
         mat(52) = -( het_rates(63) )
         mat(1) = -( het_rates(64) )
         mat(2) = -( het_rates(65) )
         mat(3) = -( het_rates(66) )
         mat(4) = -( het_rates(67) )
         mat(5) = -( het_rates(68) )
         mat(6) = -( het_rates(69) )
         mat(7) = -( het_rates(70) )
         mat(8) = -( het_rates(71) )
         mat(9) = -( het_rates(72) )
         mat(10) = -( het_rates(73) )
         mat(11) = -( het_rates(74) )
         mat(12) = -( het_rates(75) )
         mat(13) = -( het_rates(76) )
         mat(14) = -( het_rates(77) )
         mat(15) = -( het_rates(78) )
         mat(16) = -( het_rates(79) )
         mat(17) = -( het_rates(80) )
         mat(18) = -( het_rates(81) )
         mat(19) = -( het_rates(82) )
         mat(20) = -( het_rates(83) )
      end subroutine linmat02
      subroutine linmat( mat, y, rxt, het_rates )
!----------------------------------------------
! ... linear matrix entries for implicit species
!----------------------------------------------
      use chem_mods, only : gas_pcnst, rxntot, nzcnt
      use shr_kind_mod, only : r8 => shr_kind_r8
      implicit none
!----------------------------------------------
! ... dummy arguments
!----------------------------------------------
      real(r8), intent(in) :: y(gas_pcnst)
      real(r8), intent(in) :: rxt(rxntot)
      real(r8), intent(in) :: het_rates(max(1,gas_pcnst))
      real(r8), intent(inout) :: mat(nzcnt)
      call linmat01( mat, y, rxt, het_rates )
      call linmat02( mat, y, rxt, het_rates )
      end subroutine linmat
      end module mo_lin_matrix
