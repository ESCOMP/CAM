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
         mat(610) = -( rxt(3) + rxt(4) + het_rates(1) )
         mat(499) = -( rxt(62) + rxt(63) + rxt(64) + rxt(75) + rxt(76) + rxt(77) &
                 + het_rates(2) )
         mat(230) = rxt(1) + 2.000_r8*rxt(2) + rxt(68) + rxt(69) + rxt(70) &
                      + 2.000_r8*rxt(73) + rxt(80) + rxt(81) + rxt(82) + 2.000_r8*rxt(85)
         mat(605) = rxt(4)
         mat(587) = rxt(6)
         mat(380) = rxt(8)
         mat(30) = rxt(10)
         mat(519) = rxt(12)
         mat(242) = rxt(21)
         mat(319) = rxt(24)
         mat(11) = rxt(25)
         mat(293) = rxt(32)
         mat(565) = rxt(101)
         mat(51) = rxt(278)
         mat(58) = rxt(285)
         mat(568) = -( rxt(101) + rxt(105)*y(7) + rxt(106)*y(7) + rxt(108)*y(44) &
                      + rxt(109)*y(45) + rxt(110)*y(46) + rxt(111)*y(54) + rxt(112)*y(55) &
                      + rxt(113)*y(47) + rxt(114)*y(52) + rxt(115)*y(53) + rxt(116)*y(48) &
                      + rxt(117)*y(43) + rxt(118)*y(51) + rxt(119)*y(50) + rxt(120)*y(56) &
                      + rxt(121)*y(57) + rxt(122)*y(58) + rxt(123)*y(59) + rxt(126)*y(15) &
                      + rxt(127)*y(15) + rxt(128)*y(15) + het_rates(3) )
         mat(231) = rxt(1)
         mat(608) = rxt(3)
         mat(244) = rxt(20)
         mat(225) = -( rxt(1) + rxt(2) + rxt(66) + rxt(68) + rxt(69) + rxt(70) + rxt(73) &
                      + rxt(78) + rxt(80) + rxt(81) + rxt(82) + rxt(85) + het_rates(4) )
         mat(596) = rxt(4)
         mat(509) = rxt(13)
         mat(8) = rxt(96)
         mat(5) = rxt(99) + rxt(100)
         mat(555) = rxt(106)*y(7)
         mat(7) = -( rxt(93) + rxt(96) + rxt(95)*y(64) + het_rates(5) )
         mat(4) = -( rxt(99) + rxt(100) + het_rates(6) )
         mat(594) = rxt(3)
         mat(6) = rxt(93) + rxt(95)*y(64)
         mat(267) = -( het_rates(20) )
         mat(251) = rxt(18)
         mat(238) = rxt(20)
         mat(558) = rxt(128)*y(15)
         mat(64) = -( het_rates(19) )
         mat(246) = rxt(17) + rxt(18)
         mat(413) = rxt(221)*y(42)
         mat(101) = rxt(267)*y(64)
         mat(174) = -( rxt(65) + het_rates(8) )
         mat(575) = rxt(6)
         mat(105) = rxt(264)
         mat(591) = -( rxt(6) + rxt(7) + het_rates(9) )
         mat(384) = rxt(8) + .500_r8*rxt(244)
         mat(32) = rxt(10)
         mat(523) = rxt(13)
         mat(127) = rxt(274)
         mat(52) = rxt(279)
         mat(569) = 2.000_r8*rxt(105)*y(7)
         mat(376) = -( rxt(8) + rxt(244) + het_rates(10) )
         mat(29) = rxt(9) + rxt(164)
         mat(209) = rxt(11)
         mat(515) = rxt(12)
         mat(42) = rxt(15) + rxt(173)
         mat(197) = rxt(30)
         mat(73) = rxt(36)
         mat(459) = -( rxt(222)*y(42) + rxt(223)*y(49) + rxt(224)*y(47) + rxt(225)*y(43) &
                      + rxt(227)*y(52) + rxt(228)*y(53) + rxt(229)*y(59) + rxt(230)*y(58) &
                      + rxt(233)*y(15) + het_rates(22) )
         mat(210) = rxt(11)
         mat(44) = rxt(14)
         mat(37) = rxt(16)
         mat(241) = rxt(19)
         mat(80) = 2.000_r8*rxt(22)
         mat(188) = rxt(27)
         mat(117) = rxt(33)
         mat(379) = .500_r8*rxt(244)
         mat(564) = rxt(126)*y(15)
         mat(520) = -( rxt(12) + rxt(13) + rxt(243) + het_rates(11) )
         mat(31) = rxt(9) + rxt(10) + rxt(164)
         mat(45) = rxt(14)
         mat(201) = rxt(29)
         mat(75) = rxt(35)
         mat(205) = -( rxt(11) + het_rates(12) )
         mat(28) = 2.000_r8*rxt(242) + 2.000_r8*rxt(246) + 2.000_r8*rxt(252) &
                      + 2.000_r8*rxt(257)
         mat(508) = rxt(243)
         mat(368) = .500_r8*rxt(244)
         mat(194) = rxt(247) + rxt(253) + rxt(258)
         mat(70) = rxt(248) + rxt(256) + rxt(259)
         mat(39) = -( rxt(14) + rxt(15) + rxt(173) + het_rates(13) )
         mat(27) = -( rxt(9) + rxt(10) + rxt(164) + rxt(242) + rxt(246) + rxt(252) &
                      + rxt(257) + het_rates(14) )
         mat(630) = -( het_rates(16) )
         mat(571) = rxt(126)*y(15)
         mat(435) = rxt(180)*y(15)
         mat(161) = rxt(219)*y(15)
         mat(466) = rxt(233)*y(15)
         mat(33) = -( rxt(16) + het_rates(17) )
         mat(250) = -( rxt(17) + rxt(18) + het_rates(18) )
         mat(35) = rxt(16)
         mat(557) = rxt(127)*y(15) + rxt(128)*y(15)
         mat(329) = -( het_rates(21) )
         mat(36) = rxt(16)
         mat(252) = 2.000_r8*rxt(17)
         mat(239) = rxt(19) + 2.000_r8*rxt(21)
         mat(536) = rxt(28)
         mat(164) = rxt(34)
         mat(26) = rxt(57)
         mat(559) = rxt(127)*y(15)
         mat(402) = -( rxt(245) + het_rates(23) )
         mat(43) = rxt(15) + rxt(173)
         mat(562) = rxt(127)*y(15)
         mat(426) = rxt(221)*y(42) + rxt(226)*y(43)
         mat(457) = rxt(222)*y(42) + rxt(225)*y(43)
         mat(76) = -( rxt(22) + het_rates(24) )
         mat(389) = .500_r8*rxt(245)
         mat(237) = -( rxt(19) + rxt(20) + rxt(21) + het_rates(74) )
         mat(449) = rxt(222)*y(42) + rxt(223)*y(49) + rxt(224)*y(47) + rxt(225)*y(43) &
                      + rxt(229)*y(59) + rxt(233)*y(15)
         mat(427) = -( rxt(180)*y(15) + rxt(221)*y(42) + rxt(226)*y(43) + rxt(231)*y(59) &
                      + rxt(232)*y(58) + het_rates(28) )
         mat(16) = 2.000_r8*rxt(23)
         mat(317) = rxt(24)
         mat(3) = 2.000_r8*rxt(26)
         mat(187) = rxt(27)
         mat(540) = rxt(28)
         mat(198) = rxt(29)
         mat(23) = rxt(31)
         mat(19) = rxt(56)
         mat(563) = 2.000_r8*rxt(108)*y(44) + 2.000_r8*rxt(109)*y(45) &
                      + 2.000_r8*rxt(110)*y(46) + 2.000_r8*rxt(111)*y(54) + rxt(112)*y(55) &
                      + rxt(113)*y(47) + rxt(114)*y(52) + rxt(115)*y(53) &
                      + 4.000_r8*rxt(116)*y(48) + rxt(118)*y(51)
         mat(458) = rxt(222)*y(42) + 3.000_r8*rxt(223)*y(49) + rxt(224)*y(47) &
                      + rxt(227)*y(52) + rxt(228)*y(53)
         mat(15) = -( rxt(23) + het_rates(29) )
         mat(312) = -( rxt(24) + het_rates(30) )
         mat(10) = rxt(25)
         mat(196) = rxt(30)
         mat(2) = 2.000_r8*rxt(192)
         mat(9) = -( rxt(25) + het_rates(31) )
         mat(1) = -( rxt(26) + rxt(192) + het_rates(32) )
         mat(544) = -( rxt(28) + het_rates(33) )
         mat(431) = rxt(180)*y(15) + 2.000_r8*rxt(221)*y(42) + rxt(226)*y(43) &
                      + rxt(231)*y(59) + rxt(232)*y(58)
         mat(184) = -( rxt(27) + het_rates(34) )
         mat(192) = rxt(247) + rxt(253) + rxt(258)
         mat(193) = -( rxt(29) + rxt(30) + rxt(247) + rxt(253) + rxt(258) + het_rates(35) &
       )
         mat(21) = -( rxt(31) + het_rates(36) )
         mat(350) = -( het_rates(37) )
         mat(22) = rxt(31)
         mat(288) = rxt(32)
         mat(115) = rxt(33)
         mat(165) = rxt(34)
         mat(72) = rxt(35)
         mat(560) = rxt(117)*y(43) + rxt(118)*y(51) + rxt(119)*y(50) &
                      + 2.000_r8*rxt(120)*y(56) + 2.000_r8*rxt(121)*y(57) &
                      + 3.000_r8*rxt(122)*y(58) + 2.000_r8*rxt(123)*y(59)
         mat(455) = rxt(225)*y(43) + 2.000_r8*rxt(229)*y(59) + 3.000_r8*rxt(230)*y(58)
         mat(424) = rxt(226)*y(43) + 2.000_r8*rxt(231)*y(59) + 3.000_r8*rxt(232)*y(58)
         mat(285) = -( rxt(32) + het_rates(38) )
         mat(71) = rxt(36)
         mat(162) = -( rxt(34) + het_rates(39) )
         mat(112) = -( rxt(33) + het_rates(40) )
         mat(69) = rxt(248) + rxt(256) + rxt(259)
         mat(68) = -( rxt(35) + rxt(36) + rxt(248) + rxt(256) + rxt(259) + het_rates(41) &
       )
         mat(91) = -( het_rates(65) )
         mat(49) = rxt(278)
         mat(54) = rxt(285)
         mat(120) = -( rxt(274) + het_rates(66) )
         mat(221) = rxt(66) + rxt(78)
         mat(103) = rxt(267)*y(64)
         mat(83) = -( het_rates(67) )
         mat(169) = rxt(65)
         mat(48) = rxt(279)
         mat(102) = -( rxt(264) + rxt(267)*y(64) + het_rates(68) )
         mat(477) = rxt(62) + rxt(75)
         mat(220) = rxt(68) + rxt(80)
         mat(50) = rxt(287)
         mat(55) = rxt(289)
         mat(129) = -( het_rates(69) )
         mat(573) = rxt(7)
         mat(104) = rxt(264)
         mat(121) = rxt(274)
         mat(59) = -( het_rates(71) )
         mat(144) = -( het_rates(70) )
         mat(574) = rxt(7)
         mat(481) = rxt(62) + rxt(63) + rxt(64) + rxt(75) + rxt(76) + rxt(77)
         mat(173) = rxt(65)
         mat(223) = rxt(66) + rxt(68) + rxt(69) + rxt(70) + rxt(78) + rxt(80) + rxt(81) &
                      + rxt(82)
         mat(46) = -( rxt(278) + rxt(279) + rxt(287) + rxt(288) + het_rates(72) )
         mat(469) = rxt(64) + rxt(77)
         mat(215) = rxt(70) + rxt(82)
         mat(53) = -( rxt(285) + rxt(289) + het_rates(73) )
         mat(470) = rxt(63) + rxt(76)
         mat(216) = rxt(69) + rxt(81)
         mat(47) = rxt(288)
         mat(12) = -( rxt(55) + het_rates(60) )
         mat(550) = rxt(109)*y(45) + rxt(110)*y(46) + 2.000_r8*rxt(111)*y(54) &
                      + 2.000_r8*rxt(112)*y(55) + rxt(113)*y(47) + rxt(115)*y(53) &
                      + rxt(118)*y(51) + rxt(119)*y(50) + rxt(120)*y(56) &
                      + 2.000_r8*rxt(121)*y(57)
         mat(436) = rxt(224)*y(47) + rxt(228)*y(53)
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
         mat(17) = -( rxt(56) + het_rates(61) )
         mat(551) = rxt(108)*y(44) + rxt(110)*y(46) + rxt(114)*y(52)
         mat(437) = rxt(227)*y(52)
         mat(24) = -( rxt(57) + het_rates(62) )
         mat(153) = rxt(219)*y(15)
         mat(154) = -( rxt(219)*y(15) + het_rates(63) )
         mat(13) = 2.000_r8*rxt(55)
         mat(18) = rxt(56)
         mat(25) = rxt(57)
         mat(552) = rxt(112)*y(55) + rxt(119)*y(50)
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
