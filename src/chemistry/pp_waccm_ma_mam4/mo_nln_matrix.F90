      module mo_nln_matrix
      use shr_kind_mod, only : r8 => shr_kind_r8
      private
      public :: nlnmat
      contains
      subroutine nlnmat01( mat, y, rxt )
      use chem_mods, only : gas_pcnst, rxntot, nzcnt
      implicit none
!----------------------------------------------
! ... dummy arguments
!----------------------------------------------
      real(r8), intent(in) :: y(gas_pcnst)
      real(r8), intent(in) :: rxt(rxntot)
      real(r8), intent(inout) :: mat(nzcnt)
!----------------------------------------------
! ... local variables
!----------------------------------------------
!----------------------------------------------
! ... complete matrix entries implicit species
!----------------------------------------------
         mat(423) = -(rxt(93)*y(2) + rxt(111)*y(97) + rxt(137)*y(18) + rxt(142)*y(86) &
                      + rxt(150)*y(87) + rxt(163)*y(6) + rxt(166)*y(7) + rxt(178) &
                      *y(84) + rxt(205)*y(85) + rxt(254)*y(58) + rxt(257)*y(59))
         mat(724) = -rxt(93)*y(1)
         mat(683) = -rxt(111)*y(1)
         mat(274) = -rxt(137)*y(1)
         mat(510) = -rxt(142)*y(1)
         mat(616) = -rxt(150)*y(1)
         mat(446) = -rxt(163)*y(1)
         mat(473) = -rxt(166)*y(1)
         mat(352) = -rxt(178)*y(1)
         mat(400) = -rxt(205)*y(1)
         mat(135) = -rxt(254)*y(1)
         mat(264) = -rxt(257)*y(1)
         mat(724) = mat(724) + rxt(92)*y(3)
         mat(591) = rxt(92)*y(2)
         mat(735) = -(rxt(92)*y(3) + rxt(93)*y(1) + 4._r8*rxt(94)*y(2) + rxt(141) &
                      *y(86) + rxt(148)*y(17) + rxt(149)*y(87) + rxt(152)*y(19) &
                      + rxt(161)*y(6) + (rxt(164) + rxt(165)) * y(7) + rxt(172)*y(8) &
                      + rxt(185)*y(23) + rxt(198)*y(26) + rxt(199)*y(27) + rxt(202) &
                      *y(28) + rxt(208)*y(30) + rxt(218)*y(31) + rxt(219)*y(32) &
                      + rxt(220)*y(33) + rxt(242)*y(15) + rxt(250)*y(57) + (rxt(286) &
                      + rxt(287)) * y(88) + rxt(293)*y(90))
         mat(602) = -rxt(92)*y(2)
         mat(434) = -rxt(93)*y(2)
         mat(521) = -rxt(141)*y(2)
         mat(339) = -rxt(148)*y(2)
         mat(627) = -rxt(149)*y(2)
         mat(120) = -rxt(152)*y(2)
         mat(457) = -rxt(161)*y(2)
         mat(484) = -(rxt(164) + rxt(165)) * y(2)
         mat(546) = -rxt(172)*y(2)
         mat(390) = -rxt(185)*y(2)
         mat(670) = -rxt(198)*y(2)
         mat(225) = -rxt(199)*y(2)
         mat(244) = -rxt(202)*y(2)
         mat(322) = -rxt(208)*y(2)
         mat(203) = -rxt(218)*y(2)
         mat(196) = -rxt(219)*y(2)
         mat(113) = -rxt(220)*y(2)
         mat(298) = -rxt(242)*y(2)
         mat(76) = -rxt(250)*y(2)
         mat(131) = -(rxt(286) + rxt(287)) * y(2)
         mat(89) = -rxt(293)*y(2)
         mat(694) = (rxt(106)+rxt(107))*y(3)
         mat(602) = mat(602) + (rxt(106)+rxt(107))*y(97) + rxt(156)*y(5) + rxt(292) &
                      *y(90) + rxt(284)*y(91) + rxt(253)*y(58) + rxt(256)*y(59)
         mat(217) = rxt(156)*y(3) + rxt(157)*y(6) + rxt(158)*y(7) + rxt(289)*y(89)
         mat(457) = mat(457) + rxt(157)*y(5)
         mat(484) = mat(484) + rxt(158)*y(5)
         mat(521) = mat(521) + 2.000_r8*rxt(144)*y(86)
         mat(279) = rxt(140)*y(87)
         mat(627) = mat(627) + rxt(140)*y(18)
         mat(165) = rxt(289)*y(5) + 1.150_r8*rxt(297)*y(93)
         mat(89) = mat(89) + rxt(292)*y(3)
         mat(148) = rxt(284)*y(3)
         mat(173) = rxt(296)*y(93)
         mat(187) = 1.150_r8*rxt(297)*y(89) + rxt(296)*y(92)
         mat(138) = rxt(253)*y(3)
         mat(271) = rxt(256)*y(3)
         mat(693) = -((rxt(106) + rxt(107)) * y(3) + rxt(108)*y(98) + rxt(111)*y(1) &
                      + rxt(128)*y(52) + rxt(129)*y(53) + rxt(133)*y(17) + rxt(134) &
                      *y(26) + rxt(135)*y(31))
         mat(601) = -(rxt(106) + rxt(107)) * y(97)
         mat(570) = -rxt(108)*y(97)
         mat(433) = -rxt(111)*y(97)
         mat(26) = -rxt(128)*y(97)
         mat(40) = -rxt(129)*y(97)
         mat(338) = -rxt(133)*y(97)
         mat(669) = -rxt(134)*y(97)
         mat(202) = -rxt(135)*y(97)
         mat(601) = mat(601) + rxt(153)*y(94)
         mat(164) = .850_r8*rxt(297)*y(93)
         mat(101) = rxt(153)*y(3)
         mat(186) = .850_r8*rxt(297)*y(89)
         mat(597) = -(rxt(92)*y(2) + rxt(102)*y(96) + rxt(106)*y(97) + rxt(136)*y(18) &
                      + rxt(153)*y(94) + rxt(156)*y(5) + rxt(253)*y(58) + rxt(256) &
                      *y(59) + rxt(284)*y(91) + (rxt(291) + rxt(292)) * y(90) + rxt(294) &
                      *y(88))
         mat(730) = -rxt(92)*y(3)
         mat(31) = -rxt(102)*y(3)
         mat(689) = -rxt(106)*y(3)
         mat(277) = -rxt(136)*y(3)
         mat(100) = -rxt(153)*y(3)
         mat(214) = -rxt(156)*y(3)
         mat(137) = -rxt(253)*y(3)
         mat(269) = -rxt(256)*y(3)
         mat(146) = -rxt(284)*y(3)
         mat(88) = -(rxt(291) + rxt(292)) * y(3)
         mat(129) = -rxt(294)*y(3)
         mat(429) = 2.000_r8*rxt(93)*y(2) + 2.000_r8*rxt(111)*y(97) + rxt(163)*y(6) &
                      + rxt(166)*y(7) + rxt(142)*y(86) + rxt(137)*y(18) &
                      + 2.000_r8*rxt(150)*y(87) + rxt(178)*y(84) + rxt(205)*y(85) &
                      + rxt(254)*y(58) + rxt(257)*y(59)
         mat(730) = mat(730) + 2.000_r8*rxt(93)*y(1) + 2.000_r8*rxt(94)*y(2) &
                      + rxt(101)*y(96) + rxt(164)*y(7) + rxt(141)*y(86) + rxt(172) &
                      *y(8) + rxt(149)*y(87) + rxt(185)*y(23) + rxt(208)*y(30)
         mat(689) = mat(689) + 2.000_r8*rxt(111)*y(1)
         mat(597) = mat(597) + 2.000_r8*rxt(102)*y(96)
         mat(31) = mat(31) + rxt(101)*y(2) + 2.000_r8*rxt(102)*y(3)
         mat(214) = mat(214) + rxt(160)*y(7)
         mat(452) = rxt(163)*y(1) + rxt(290)*y(89)
         mat(479) = rxt(166)*y(1) + rxt(164)*y(2) + rxt(160)*y(5)
         mat(516) = rxt(142)*y(1) + rxt(141)*y(2) + rxt(176)*y(10) + rxt(143)*y(87) &
                      + rxt(187)*y(23)
         mat(541) = rxt(172)*y(2) + rxt(174)*y(87)
         mat(95) = rxt(176)*y(86)
         mat(641) = rxt(244)*y(87)
         mat(277) = mat(277) + rxt(137)*y(1) + rxt(139)*y(87)
         mat(622) = 2.000_r8*rxt(150)*y(1) + rxt(149)*y(2) + rxt(143)*y(86) + rxt(174) &
                      *y(8) + rxt(244)*y(13) + rxt(139)*y(18) + 2.000_r8*rxt(151) &
                      *y(87) + rxt(181)*y(84) + rxt(188)*y(23) + rxt(206)*y(85) &
                      + rxt(210)*y(30)
         mat(357) = rxt(178)*y(1) + rxt(181)*y(87)
         mat(385) = rxt(185)*y(2) + rxt(187)*y(86) + rxt(188)*y(87) + ( &
                      + 2.000_r8*rxt(192)+2.000_r8*rxt(193))*y(23) + (rxt(214) &
                       +rxt(215))*y(30)
         mat(406) = rxt(205)*y(1) + rxt(206)*y(87)
         mat(318) = rxt(208)*y(2) + rxt(210)*y(87) + (rxt(214)+rxt(215))*y(23) &
                      + 2.000_r8*rxt(216)*y(30)
         mat(163) = rxt(290)*y(6)
         mat(137) = mat(137) + rxt(254)*y(1)
         mat(269) = mat(269) + rxt(257)*y(1)
         mat(33) = -(rxt(95)*y(2) + rxt(96)*y(3) + rxt(98)*y(1))
         mat(696) = -rxt(95)*y(95)
         mat(573) = -rxt(96)*y(95)
         mat(413) = -rxt(98)*y(95)
         mat(672) = rxt(106)*y(3)
         mat(573) = mat(573) + rxt(106)*y(97)
         mat(30) = -(rxt(101)*y(2) + rxt(102)*y(3))
         mat(695) = -rxt(101)*y(96)
         mat(572) = -rxt(102)*y(96)
         mat(412) = rxt(98)*y(95)
         mat(695) = mat(695) + rxt(95)*y(95)
         mat(572) = mat(572) + rxt(96)*y(95)
         mat(32) = rxt(98)*y(1) + rxt(95)*y(2) + rxt(96)*y(3)
         mat(327) = -(rxt(133)*y(97) + rxt(146)*y(86) + rxt(148)*y(2) + rxt(179)*y(84) &
                      + rxt(222)*y(55))
         mat(679) = -rxt(133)*y(17)
         mat(506) = -rxt(146)*y(17)
         mat(720) = -rxt(148)*y(17)
         mat(348) = -rxt(179)*y(17)
         mat(153) = -rxt(222)*y(17)
         mat(273) = rxt(139)*y(87)
         mat(612) = rxt(139)*y(18)
         mat(102) = -((rxt(238) + rxt(239)) * y(86))
         mat(492) = -(rxt(238) + rxt(239)) * y(16)
         mat(700) = rxt(242)*y(15) + rxt(250)*y(57)
         mat(492) = mat(492) + rxt(241)*y(15) + rxt(251)*y(57)
         mat(524) = rxt(240)*y(15)
         mat(280) = rxt(242)*y(2) + rxt(241)*y(86) + rxt(240)*y(8) + rxt(183)*y(84) &
                      + rxt(207)*y(85)
         mat(341) = rxt(183)*y(15)
         mat(391) = rxt(207)*y(15)
         mat(70) = rxt(250)*y(2) + rxt(251)*y(86)
         mat(209) = -(rxt(155)*y(86) + rxt(156)*y(3) + rxt(157)*y(6) + (rxt(158) &
                      + rxt(159) + rxt(160)) * y(7) + rxt(289)*y(89))
         mat(497) = -rxt(155)*y(5)
         mat(582) = -rxt(156)*y(5)
         mat(438) = -rxt(157)*y(5)
         mat(462) = -(rxt(158) + rxt(159) + rxt(160)) * y(5)
         mat(161) = -rxt(289)*y(5)
         mat(711) = rxt(293)*y(90) + rxt(154)*y(94)
         mat(582) = mat(582) + rxt(291)*y(90)
         mat(127) = 1.100_r8*rxt(298)*y(93)
         mat(87) = rxt(293)*y(2) + rxt(291)*y(3)
         mat(169) = .200_r8*rxt(296)*y(93)
         mat(98) = rxt(154)*y(2)
         mat(180) = 1.100_r8*rxt(298)*y(88) + .200_r8*rxt(296)*y(92)
         mat(447) = -(rxt(157)*y(5) + rxt(161)*y(2) + rxt(162)*y(87) + rxt(163)*y(1) &
                      + rxt(171)*y(8) + rxt(190)*y(23) + rxt(211)*y(30) + rxt(243) &
                      *y(13) + rxt(290)*y(89))
         mat(211) = -rxt(157)*y(6)
         mat(725) = -rxt(161)*y(6)
         mat(617) = -rxt(162)*y(6)
         mat(424) = -rxt(163)*y(6)
         mat(536) = -rxt(171)*y(6)
         mat(380) = -rxt(190)*y(6)
         mat(313) = -rxt(211)*y(6)
         mat(636) = -rxt(243)*y(6)
         mat(162) = -rxt(290)*y(6)
         mat(725) = mat(725) + rxt(164)*y(7)
         mat(592) = rxt(156)*y(5) + rxt(153)*y(94)
         mat(211) = mat(211) + rxt(156)*y(3) + 2.000_r8*rxt(159)*y(7) + rxt(155)*y(86)
         mat(474) = rxt(164)*y(2) + 2.000_r8*rxt(159)*y(5) + rxt(258)*y(59)
         mat(511) = rxt(155)*y(5)
         mat(99) = rxt(153)*y(3)
         mat(265) = rxt(258)*y(7)
         mat(475) = -((rxt(158) + rxt(159) + rxt(160)) * y(5) + (rxt(164) + rxt(165) &
                      ) * y(2) + rxt(166)*y(1) + rxt(167)*y(8) + rxt(169)*y(86) &
                      + rxt(175)*y(87) + rxt(191)*y(23) + rxt(212)*y(30) + rxt(258) &
                      *y(59))
         mat(212) = -(rxt(158) + rxt(159) + rxt(160)) * y(7)
         mat(726) = -(rxt(164) + rxt(165)) * y(7)
         mat(425) = -rxt(166)*y(7)
         mat(537) = -rxt(167)*y(7)
         mat(512) = -rxt(169)*y(7)
         mat(618) = -rxt(175)*y(7)
         mat(381) = -rxt(191)*y(7)
         mat(314) = -rxt(212)*y(7)
         mat(266) = -rxt(258)*y(7)
         mat(425) = mat(425) + rxt(163)*y(6)
         mat(726) = mat(726) + rxt(161)*y(6) + rxt(172)*y(8)
         mat(448) = rxt(163)*y(1) + rxt(161)*y(2) + 2.000_r8*rxt(171)*y(8) + rxt(243) &
                      *y(13) + rxt(162)*y(87) + rxt(190)*y(23) + rxt(211)*y(30)
         mat(512) = mat(512) + rxt(173)*y(8) + rxt(176)*y(10)
         mat(537) = mat(537) + rxt(172)*y(2) + 2.000_r8*rxt(171)*y(6) + rxt(173)*y(86) &
                      + rxt(174)*y(87)
         mat(91) = rxt(176)*y(86)
         mat(637) = rxt(243)*y(6)
         mat(618) = mat(618) + rxt(162)*y(6) + rxt(174)*y(8)
         mat(381) = mat(381) + rxt(190)*y(6)
         mat(314) = mat(314) + rxt(211)*y(6)
      end subroutine nlnmat01
      subroutine nlnmat02( mat, y, rxt )
      use chem_mods, only : gas_pcnst, rxntot, nzcnt
      implicit none
!----------------------------------------------
! ... dummy arguments
!----------------------------------------------
      real(r8), intent(in) :: y(gas_pcnst)
      real(r8), intent(in) :: rxt(rxntot)
      real(r8), intent(inout) :: mat(nzcnt)
!----------------------------------------------
! ... local variables
!----------------------------------------------
!----------------------------------------------
! ... complete matrix entries implicit species
!----------------------------------------------
         mat(513) = -(rxt(141)*y(2) + rxt(142)*y(1) + rxt(143)*y(87) + (4._r8*rxt(144) &
                      + 4._r8*rxt(145)) * y(86) + rxt(146)*y(17) + rxt(147)*y(19) &
                      + rxt(155)*y(5) + rxt(169)*y(7) + rxt(170)*y(9) + rxt(173)*y(8) &
                      + rxt(176)*y(10) + (rxt(186) + rxt(187)) * y(23) + rxt(197) &
                      *y(26) + rxt(201)*y(27) + rxt(203)*y(28) + rxt(209)*y(30) &
                      + rxt(217)*y(31) + (rxt(238) + rxt(239)) * y(16) + rxt(241) &
                      *y(15) + rxt(245)*y(14) + rxt(251)*y(57) + rxt(252)*y(58) &
                      + rxt(255)*y(59) + rxt(262)*y(60) + (rxt(264) + rxt(265) &
                      ) * y(63))
         mat(727) = -rxt(141)*y(86)
         mat(426) = -rxt(142)*y(86)
         mat(619) = -rxt(143)*y(86)
         mat(331) = -rxt(146)*y(86)
         mat(116) = -rxt(147)*y(86)
         mat(213) = -rxt(155)*y(86)
         mat(476) = -rxt(169)*y(86)
         mat(251) = -rxt(170)*y(86)
         mat(538) = -rxt(173)*y(86)
         mat(92) = -rxt(176)*y(86)
         mat(382) = -(rxt(186) + rxt(187)) * y(86)
         mat(662) = -rxt(197)*y(86)
         mat(222) = -rxt(201)*y(86)
         mat(240) = -rxt(203)*y(86)
         mat(315) = -rxt(209)*y(86)
         mat(200) = -rxt(217)*y(86)
         mat(104) = -(rxt(238) + rxt(239)) * y(86)
         mat(290) = -rxt(241)*y(86)
         mat(80) = -rxt(245)*y(86)
         mat(75) = -rxt(251)*y(86)
         mat(136) = -rxt(252)*y(86)
         mat(267) = -rxt(255)*y(86)
         mat(229) = -rxt(262)*y(86)
         mat(55) = -(rxt(264) + rxt(265)) * y(86)
         mat(426) = mat(426) + rxt(137)*y(18) + rxt(150)*y(87)
         mat(727) = mat(727) + rxt(148)*y(17) + rxt(242)*y(15) + rxt(149)*y(87) &
                      + rxt(152)*y(19) + rxt(198)*y(26) + rxt(199)*y(27) + rxt(218) &
                      *y(31) + rxt(219)*y(32)
         mat(686) = rxt(133)*y(17) + 2.000_r8*rxt(108)*y(98) + rxt(134)*y(26) &
                      + rxt(135)*y(31)
         mat(331) = mat(331) + rxt(148)*y(2) + rxt(133)*y(97)
         mat(449) = rxt(162)*y(87)
         mat(538) = mat(538) + rxt(174)*y(87)
         mat(290) = mat(290) + rxt(242)*y(2)
         mat(275) = rxt(137)*y(1) + 2.000_r8*rxt(138)*y(87)
         mat(619) = mat(619) + rxt(150)*y(1) + rxt(149)*y(2) + rxt(162)*y(6) &
                      + rxt(174)*y(8) + 2.000_r8*rxt(138)*y(18) + rxt(182)*y(84)
         mat(116) = mat(116) + rxt(152)*y(2)
         mat(563) = 2.000_r8*rxt(108)*y(97) + rxt(221)*y(55)
         mat(354) = rxt(182)*y(87)
         mat(662) = mat(662) + rxt(198)*y(2) + rxt(134)*y(97)
         mat(222) = mat(222) + rxt(199)*y(2)
         mat(200) = mat(200) + rxt(218)*y(2) + rxt(135)*y(97)
         mat(193) = rxt(219)*y(2)
         mat(154) = rxt(221)*y(98)
         mat(539) = -(rxt(167)*y(7) + rxt(171)*y(6) + rxt(172)*y(2) + rxt(173)*y(86) &
                      + rxt(174)*y(87) + rxt(240)*y(15) + rxt(266)*y(63))
         mat(477) = -rxt(167)*y(8)
         mat(450) = -rxt(171)*y(8)
         mat(728) = -rxt(172)*y(8)
         mat(514) = -rxt(173)*y(8)
         mat(620) = -rxt(174)*y(8)
         mat(291) = -rxt(240)*y(8)
         mat(56) = -rxt(266)*y(8)
         mat(427) = rxt(166)*y(7)
         mat(728) = mat(728) + rxt(165)*y(7) + rxt(202)*y(28) + rxt(220)*y(33)
         mat(477) = mat(477) + rxt(166)*y(1) + rxt(165)*y(2)
         mat(514) = mat(514) + rxt(170)*y(9) + rxt(203)*y(28)
         mat(252) = rxt(170)*y(86) + rxt(224)*y(55)
         mat(355) = rxt(204)*y(28)
         mat(241) = rxt(202)*y(2) + rxt(203)*y(86) + rxt(204)*y(84)
         mat(112) = rxt(220)*y(2)
         mat(155) = rxt(224)*y(9)
         mat(247) = -(rxt(170)*y(86) + rxt(224)*y(55))
         mat(501) = -rxt(170)*y(9)
         mat(151) = -rxt(224)*y(9)
         mat(465) = rxt(169)*y(86)
         mat(501) = mat(501) + rxt(169)*y(7)
         mat(526) = rxt(240)*y(15) + rxt(266)*y(63)
         mat(282) = rxt(240)*y(8)
         mat(652) = (rxt(270)+rxt(275)+rxt(281))*y(28)
         mat(236) = (rxt(270)+rxt(275)+rxt(281))*y(26)
         mat(54) = rxt(266)*y(8)
         mat(90) = -(rxt(176)*y(86))
         mat(491) = -rxt(176)*y(10)
         mat(459) = rxt(175)*y(87)
         mat(604) = rxt(175)*y(7)
         mat(458) = rxt(167)*y(8)
         mat(523) = rxt(167)*y(7)
         mat(643) = -(rxt(189)*y(23) + rxt(243)*y(6) + rxt(244)*y(87))
         mat(387) = -rxt(189)*y(13)
         mat(454) = -rxt(243)*y(13)
         mat(624) = -rxt(244)*y(13)
         mat(518) = rxt(245)*y(14)
         mat(82) = rxt(245)*y(86)
         mat(77) = -(rxt(245)*y(86))
         mat(490) = -rxt(245)*y(14)
         mat(628) = rxt(244)*y(87)
         mat(603) = rxt(244)*y(13)
         mat(284) = -(rxt(183)*y(84) + rxt(207)*y(85) + rxt(240)*y(8) + rxt(241)*y(86) &
                      + rxt(242)*y(2))
         mat(347) = -rxt(183)*y(15)
         mat(394) = -rxt(207)*y(15)
         mat(529) = -rxt(240)*y(15)
         mat(504) = -rxt(241)*y(15)
         mat(718) = -rxt(242)*y(15)
         mat(440) = rxt(243)*y(13)
         mat(630) = rxt(243)*y(6) + rxt(189)*y(23)
         mat(373) = rxt(189)*y(13)
         mat(272) = -(rxt(136)*y(3) + rxt(137)*y(1) + (rxt(138) + rxt(139) + rxt(140) &
                      ) * y(87))
         mat(585) = -rxt(136)*y(18)
         mat(417) = -rxt(137)*y(18)
         mat(609) = -(rxt(138) + rxt(139) + rxt(140)) * y(18)
         mat(717) = rxt(148)*y(17) + rxt(141)*y(86)
         mat(677) = rxt(133)*y(17)
         mat(326) = rxt(148)*y(2) + rxt(133)*y(97) + rxt(146)*y(86) + rxt(179)*y(84) &
                      + rxt(222)*y(55)
         mat(103) = rxt(238)*y(86)
         mat(210) = rxt(155)*y(86)
         mat(503) = rxt(141)*y(2) + rxt(146)*y(17) + rxt(238)*y(16) + rxt(155)*y(5) &
                      + rxt(241)*y(15) + rxt(251)*y(57) + rxt(252)*y(58) + rxt(255) &
                      *y(59)
         mat(283) = rxt(241)*y(86)
         mat(346) = rxt(179)*y(17)
         mat(152) = rxt(222)*y(17)
         mat(74) = rxt(251)*y(86)
         mat(134) = rxt(252)*y(86)
         mat(259) = rxt(255)*y(86)
         mat(623) = -((rxt(138) + rxt(139) + rxt(140)) * y(18) + rxt(143)*y(86) &
                      + rxt(149)*y(2) + rxt(150)*y(1) + 4._r8*rxt(151)*y(87) + rxt(162) &
                      *y(6) + rxt(174)*y(8) + rxt(175)*y(7) + (rxt(181) + rxt(182) &
                      ) * y(84) + rxt(188)*y(23) + rxt(206)*y(85) + rxt(210)*y(30) &
                      + rxt(244)*y(13))
         mat(278) = -(rxt(138) + rxt(139) + rxt(140)) * y(87)
         mat(517) = -rxt(143)*y(87)
         mat(731) = -rxt(149)*y(87)
         mat(430) = -rxt(150)*y(87)
         mat(453) = -rxt(162)*y(87)
         mat(542) = -rxt(174)*y(87)
         mat(480) = -rxt(175)*y(87)
         mat(358) = -(rxt(181) + rxt(182)) * y(87)
         mat(386) = -rxt(188)*y(87)
         mat(407) = -rxt(206)*y(87)
         mat(319) = -rxt(210)*y(87)
         mat(642) = -rxt(244)*y(87)
         mat(430) = mat(430) + rxt(142)*y(86)
         mat(731) = mat(731) + rxt(242)*y(15) + rxt(152)*y(19)
         mat(598) = rxt(136)*y(18)
         mat(105) = rxt(239)*y(86)
         mat(453) = mat(453) + rxt(243)*y(13)
         mat(517) = mat(517) + rxt(142)*y(1) + rxt(239)*y(16) + rxt(173)*y(8) &
                      + rxt(147)*y(19) + rxt(186)*y(23) + rxt(209)*y(30) + rxt(262) &
                      *y(60) + .500_r8*rxt(264)*y(63)
         mat(542) = mat(542) + rxt(173)*y(86) + rxt(240)*y(15)
         mat(642) = mat(642) + rxt(243)*y(6) + rxt(189)*y(23)
         mat(294) = rxt(242)*y(2) + rxt(240)*y(8) + rxt(183)*y(84) + rxt(207)*y(85)
         mat(278) = mat(278) + rxt(136)*y(3)
         mat(118) = rxt(152)*y(2) + rxt(147)*y(86) + rxt(180)*y(84)
         mat(358) = mat(358) + rxt(183)*y(15) + rxt(180)*y(19)
         mat(386) = mat(386) + rxt(186)*y(86) + rxt(189)*y(13)
         mat(407) = mat(407) + rxt(207)*y(15)
         mat(319) = mat(319) + rxt(209)*y(86)
         mat(231) = rxt(262)*y(86)
         mat(57) = .500_r8*rxt(264)*y(86)
         mat(114) = -(rxt(147)*y(86) + rxt(152)*y(2) + rxt(180)*y(84))
         mat(493) = -rxt(147)*y(19)
         mat(702) = -rxt(152)*y(19)
         mat(342) = -rxt(180)*y(19)
         mat(493) = mat(493) + 2.000_r8*rxt(145)*y(86)
         mat(605) = 2.000_r8*rxt(151)*y(87)
         mat(565) = -(rxt(108)*y(97) + rxt(221)*y(55) + rxt(263)*y(61))
         mat(688) = -rxt(108)*y(98)
         mat(156) = -rxt(221)*y(98)
         mat(50) = -rxt(263)*y(98)
         mat(333) = rxt(146)*y(86)
         mat(515) = rxt(146)*y(17) + 2.000_r8*rxt(144)*y(86) + rxt(170)*y(9) &
                      + rxt(176)*y(10) + rxt(245)*y(14) + rxt(241)*y(15) + rxt(143) &
                      *y(87) + rxt(147)*y(19) + rxt(197)*y(26) + rxt(201)*y(27) &
                      + rxt(217)*y(31)
         mat(253) = rxt(170)*y(86)
         mat(94) = rxt(176)*y(86)
         mat(81) = rxt(245)*y(86)
         mat(292) = rxt(241)*y(86)
         mat(276) = rxt(140)*y(87)
         mat(621) = rxt(143)*y(86) + rxt(140)*y(18)
         mat(117) = rxt(147)*y(86)
         mat(664) = rxt(197)*y(86) + (rxt(271)+rxt(276)+rxt(282))*y(27) + (rxt(272) &
                       +rxt(283))*y(32)
         mat(223) = rxt(201)*y(86) + (rxt(271)+rxt(276)+rxt(282))*y(26)
         mat(201) = rxt(217)*y(86)
         mat(194) = (rxt(272)+rxt(283))*y(26)
         mat(349) = -(rxt(178)*y(1) + rxt(179)*y(17) + rxt(180)*y(19) + (rxt(181) &
                      + rxt(182)) * y(87) + rxt(183)*y(15) + rxt(200)*y(27) + rxt(204) &
                      *y(28))
         mat(420) = -rxt(178)*y(84)
         mat(328) = -rxt(179)*y(84)
         mat(115) = -rxt(180)*y(84)
         mat(613) = -(rxt(181) + rxt(182)) * y(84)
         mat(286) = -rxt(183)*y(84)
         mat(220) = -rxt(200)*y(84)
         mat(237) = -rxt(204)*y(84)
         mat(721) = rxt(185)*y(23) + rxt(198)*y(26)
         mat(680) = rxt(134)*y(26) + rxt(129)*y(53)
         mat(443) = rxt(190)*y(23)
         mat(507) = rxt(186)*y(23) + rxt(197)*y(26)
         mat(632) = rxt(189)*y(23)
         mat(376) = rxt(185)*y(2) + rxt(190)*y(6) + rxt(186)*y(86) + rxt(189)*y(13) + ( &
                      + 4.000_r8*rxt(192)+2.000_r8*rxt(194))*y(23) + rxt(214)*y(30) &
                      + rxt(259)*y(59)
         mat(656) = rxt(198)*y(2) + rxt(134)*y(97) + rxt(197)*y(86)
         mat(309) = rxt(214)*y(23)
         mat(39) = rxt(129)*y(97)
         mat(261) = rxt(259)*y(23)
      end subroutine nlnmat02
      subroutine nlnmat03( mat, y, rxt )
      use chem_mods, only : gas_pcnst, rxntot, nzcnt
      implicit none
!----------------------------------------------
! ... dummy arguments
!----------------------------------------------
      real(r8), intent(in) :: y(gas_pcnst)
      real(r8), intent(in) :: rxt(rxntot)
      real(r8), intent(inout) :: mat(nzcnt)
!----------------------------------------------
! ... local variables
!----------------------------------------------
!----------------------------------------------
! ... complete matrix entries implicit species
!----------------------------------------------
         mat(340) = rxt(204)*y(28)
         mat(364) = 2.000_r8*rxt(193)*y(23)
         mat(647) = (rxt(271)+rxt(276)+rxt(282))*y(27) + (rxt(270)+rxt(275)+rxt(281)) &
                      *y(28)
         mat(218) = (rxt(271)+rxt(276)+rxt(282))*y(26)
         mat(233) = rxt(204)*y(84) + (rxt(270)+rxt(275)+rxt(281))*y(26)
         mat(377) = -(rxt(185)*y(2) + (rxt(186) + rxt(187)) * y(86) + rxt(188)*y(87) &
                      + rxt(189)*y(13) + rxt(190)*y(6) + rxt(191)*y(7) + (4._r8*rxt(192) &
                      + 4._r8*rxt(193) + 4._r8*rxt(194) + 4._r8*rxt(195)) * y(23) &
                      + (rxt(213) + rxt(214) + rxt(215)) * y(30) + rxt(259)*y(59))
         mat(722) = -rxt(185)*y(23)
         mat(508) = -(rxt(186) + rxt(187)) * y(23)
         mat(614) = -rxt(188)*y(23)
         mat(633) = -rxt(189)*y(23)
         mat(444) = -rxt(190)*y(23)
         mat(471) = -rxt(191)*y(23)
         mat(310) = -(rxt(213) + rxt(214) + rxt(215)) * y(23)
         mat(262) = -rxt(259)*y(23)
         mat(421) = rxt(178)*y(84)
         mat(722) = mat(722) + rxt(199)*y(27) + rxt(202)*y(28)
         mat(508) = mat(508) + rxt(201)*y(27)
         mat(614) = mat(614) + rxt(182)*y(84)
         mat(350) = rxt(178)*y(1) + rxt(182)*y(87) + rxt(200)*y(27)
         mat(67) = rxt(261)*y(59)
         mat(221) = rxt(199)*y(2) + rxt(201)*y(86) + rxt(200)*y(84)
         mat(238) = rxt(202)*y(2)
         mat(262) = mat(262) + rxt(261)*y(24)
         mat(64) = -(rxt(261)*y(59))
         mat(255) = -rxt(261)*y(24)
         mat(366) = 2.000_r8*rxt(194)*y(23) + rxt(213)*y(30)
         mat(300) = rxt(213)*y(23)
         mat(363) = 2.000_r8*rxt(195)*y(23)
         mat(668) = -(rxt(134)*y(97) + rxt(197)*y(86) + rxt(198)*y(2) + (rxt(270) &
                      + rxt(275) + rxt(281)) * y(28) + (rxt(271) + rxt(276) + rxt(282) &
                      ) * y(27) + (rxt(272) + rxt(283)) * y(32))
         mat(692) = -rxt(134)*y(26)
         mat(519) = -rxt(197)*y(26)
         mat(733) = -rxt(198)*y(26)
         mat(243) = -(rxt(270) + rxt(275) + rxt(281)) * y(26)
         mat(224) = -(rxt(271) + rxt(276) + rxt(282)) * y(26)
         mat(195) = -(rxt(272) + rxt(283)) * y(26)
         mat(337) = rxt(179)*y(84)
         mat(519) = mat(519) + rxt(187)*y(23)
         mat(296) = rxt(183)*y(84)
         mat(625) = rxt(181)*y(84)
         mat(119) = rxt(180)*y(84)
         mat(360) = rxt(179)*y(17) + rxt(183)*y(15) + rxt(181)*y(87) + rxt(180)*y(19) &
                      + rxt(200)*y(27)
         mat(388) = rxt(187)*y(86)
         mat(224) = mat(224) + rxt(200)*y(84)
         mat(219) = -(rxt(199)*y(2) + rxt(200)*y(84) + rxt(201)*y(86) + (rxt(271) &
                      + rxt(276) + rxt(282)) * y(26))
         mat(712) = -rxt(199)*y(27)
         mat(343) = -rxt(200)*y(27)
         mat(498) = -rxt(201)*y(27)
         mat(650) = -(rxt(271) + rxt(276) + rxt(282)) * y(27)
         mat(498) = mat(498) + rxt(203)*y(28)
         mat(608) = rxt(188)*y(23)
         mat(367) = rxt(188)*y(87)
         mat(234) = rxt(203)*y(86)
         mat(235) = -(rxt(202)*y(2) + rxt(203)*y(86) + rxt(204)*y(84) + (rxt(270) &
                      + rxt(275) + rxt(281)) * y(26))
         mat(714) = -rxt(202)*y(28)
         mat(500) = -rxt(203)*y(28)
         mat(344) = -rxt(204)*y(28)
         mat(651) = -(rxt(270) + rxt(275) + rxt(281)) * y(28)
         mat(464) = rxt(191)*y(23)
         mat(369) = rxt(191)*y(7)
         mat(365) = rxt(215)*y(30)
         mat(648) = (rxt(272)+rxt(283))*y(32)
         mat(299) = rxt(215)*y(23)
         mat(188) = (rxt(272)+rxt(283))*y(26)
         mat(399) = -(rxt(205)*y(1) + rxt(206)*y(87) + rxt(207)*y(15))
         mat(422) = -rxt(205)*y(85)
         mat(615) = -rxt(206)*y(85)
         mat(287) = -rxt(207)*y(85)
         mat(723) = rxt(208)*y(30) + rxt(218)*y(31)
         mat(682) = rxt(135)*y(31)
         mat(445) = rxt(211)*y(30)
         mat(509) = rxt(209)*y(30) + rxt(217)*y(31)
         mat(378) = (rxt(213)+rxt(214))*y(30)
         mat(311) = rxt(208)*y(2) + rxt(211)*y(6) + rxt(209)*y(86) + (rxt(213) &
                       +rxt(214))*y(23) + 4.000_r8*rxt(216)*y(30) + rxt(260)*y(59)
         mat(199) = rxt(218)*y(2) + rxt(135)*y(97) + rxt(217)*y(86)
         mat(263) = rxt(260)*y(30)
         mat(307) = -(rxt(208)*y(2) + rxt(209)*y(86) + rxt(210)*y(87) + rxt(211)*y(6) &
                      + rxt(212)*y(7) + (rxt(213) + rxt(214) + rxt(215)) * y(23) &
                      + 4._r8*rxt(216)*y(30) + rxt(260)*y(59))
         mat(719) = -rxt(208)*y(30)
         mat(505) = -rxt(209)*y(30)
         mat(611) = -rxt(210)*y(30)
         mat(441) = -rxt(211)*y(30)
         mat(468) = -rxt(212)*y(30)
         mat(374) = -(rxt(213) + rxt(214) + rxt(215)) * y(30)
         mat(260) = -rxt(260)*y(30)
         mat(418) = rxt(205)*y(85)
         mat(719) = mat(719) + rxt(219)*y(32) + rxt(220)*y(33)
         mat(395) = rxt(205)*y(1)
         mat(190) = rxt(219)*y(2)
         mat(109) = rxt(220)*y(2)
         mat(197) = -(rxt(135)*y(97) + rxt(217)*y(86) + rxt(218)*y(2))
         mat(675) = -rxt(135)*y(31)
         mat(496) = -rxt(217)*y(31)
         mat(710) = -rxt(218)*y(31)
         mat(281) = rxt(207)*y(85)
         mat(607) = rxt(206)*y(85)
         mat(392) = rxt(207)*y(15) + rxt(206)*y(87)
         mat(189) = -(rxt(219)*y(2) + (rxt(272) + rxt(283)) * y(26))
         mat(709) = -rxt(219)*y(32)
         mat(649) = -(rxt(272) + rxt(283)) * y(32)
         mat(606) = rxt(210)*y(30)
         mat(302) = rxt(210)*y(87)
         mat(106) = -(rxt(220)*y(2))
         mat(701) = -rxt(220)*y(33)
         mat(460) = rxt(212)*y(30)
         mat(301) = rxt(212)*y(7)
         mat(122) = -((rxt(286) + rxt(287)) * y(2) + rxt(294)*y(3) + rxt(298)*y(93))
         mat(703) = -(rxt(286) + rxt(287)) * y(88)
         mat(576) = -rxt(294)*y(88)
         mat(175) = -rxt(298)*y(88)
         mat(158) = -(rxt(289)*y(5) + rxt(290)*y(6) + rxt(297)*y(93))
         mat(206) = -rxt(289)*y(89)
         mat(435) = -rxt(290)*y(89)
         mat(177) = -rxt(297)*y(89)
         mat(579) = rxt(294)*y(88) + rxt(291)*y(90) + rxt(284)*y(91)
         mat(124) = rxt(294)*y(3)
         mat(85) = rxt(291)*y(3)
         mat(141) = rxt(284)*y(3)
         mat(83) = -((rxt(291) + rxt(292)) * y(3) + rxt(293)*y(2))
         mat(574) = -(rxt(291) + rxt(292)) * y(90)
         mat(698) = -rxt(293)*y(90)
         mat(140) = -(rxt(284)*y(3))
         mat(578) = -rxt(284)*y(91)
         mat(705) = rxt(287)*y(88) + rxt(293)*y(90)
         mat(123) = rxt(287)*y(2)
         mat(84) = rxt(293)*y(2)
         mat(167) = -(rxt(296)*y(93))
         mat(178) = -rxt(296)*y(92)
         mat(707) = rxt(286)*y(88)
         mat(580) = rxt(292)*y(90)
         mat(207) = rxt(289)*y(89)
         mat(436) = rxt(290)*y(89)
         mat(125) = rxt(286)*y(2)
         mat(159) = rxt(289)*y(5) + rxt(290)*y(6)
         mat(86) = rxt(292)*y(3)
         mat(97) = -(rxt(153)*y(3) + rxt(154)*y(2))
         mat(575) = -rxt(153)*y(94)
         mat(699) = -rxt(154)*y(94)
         mat(699) = mat(699) + rxt(286)*y(88)
         mat(121) = rxt(286)*y(2) + .900_r8*rxt(298)*y(93)
         mat(166) = .800_r8*rxt(296)*y(93)
         mat(174) = .900_r8*rxt(298)*y(88) + .800_r8*rxt(296)*y(92)
         mat(179) = -(rxt(296)*y(92) + rxt(297)*y(89) + rxt(298)*y(88))
         mat(168) = -rxt(296)*y(93)
         mat(160) = -rxt(297)*y(93)
         mat(126) = -rxt(298)*y(93)
         mat(24) = -(rxt(128)*y(97))
         mat(671) = -rxt(128)*y(52)
         mat(37) = -(rxt(129)*y(97))
         mat(673) = -rxt(129)*y(53)
         mat(323) = rxt(222)*y(55)
         mat(245) = rxt(224)*y(55)
         mat(548) = rxt(221)*y(55)
         mat(149) = rxt(222)*y(17) + rxt(224)*y(9) + rxt(221)*y(98)
         mat(150) = -(rxt(221)*y(98) + rxt(222)*y(17) + rxt(224)*y(9))
         mat(550) = -rxt(221)*y(55)
         mat(324) = -rxt(222)*y(55)
         mat(246) = -rxt(224)*y(55)
         mat(674) = 2.000_r8*rxt(128)*y(52) + rxt(129)*y(53)
         mat(25) = 2.000_r8*rxt(128)*y(97)
         mat(38) = rxt(129)*y(97)
         mat(69) = -(rxt(250)*y(2) + rxt(251)*y(86))
         mat(697) = -rxt(250)*y(57)
         mat(489) = -rxt(251)*y(57)
         mat(132) = -(rxt(252)*y(86) + rxt(253)*y(3) + rxt(254)*y(1))
         mat(494) = -rxt(252)*y(58)
         mat(577) = -rxt(253)*y(58)
         mat(414) = -rxt(254)*y(58)
         mat(258) = -(rxt(255)*y(86) + rxt(256)*y(3) + rxt(257)*y(1) + rxt(258)*y(7) &
                      + rxt(259)*y(23) + rxt(260)*y(30) + rxt(261)*y(24))
         mat(502) = -rxt(255)*y(59)
         mat(584) = -rxt(256)*y(59)
         mat(416) = -rxt(257)*y(59)
         mat(466) = -rxt(258)*y(59)
         mat(371) = -rxt(259)*y(59)
         mat(305) = -rxt(260)*y(59)
         mat(66) = -rxt(261)*y(59)
         mat(416) = mat(416) + rxt(254)*y(58)
         mat(716) = rxt(250)*y(57)
         mat(584) = mat(584) + rxt(253)*y(58)
         mat(502) = mat(502) + rxt(252)*y(58)
         mat(73) = rxt(250)*y(2)
         mat(133) = rxt(254)*y(1) + rxt(253)*y(3) + rxt(252)*y(86)
         mat(227) = -(rxt(262)*y(86))
         mat(499) = -rxt(262)*y(60)
         mat(415) = rxt(257)*y(59)
         mat(583) = rxt(256)*y(59)
         mat(463) = rxt(258)*y(59)
         mat(499) = mat(499) + rxt(251)*y(57) + rxt(255)*y(59) + (.500_r8*rxt(264) &
                       +rxt(265))*y(63)
         mat(525) = rxt(266)*y(63)
         mat(368) = rxt(259)*y(59)
         mat(65) = rxt(261)*y(59)
         mat(303) = rxt(260)*y(59)
         mat(72) = rxt(251)*y(86)
         mat(257) = rxt(257)*y(1) + rxt(256)*y(3) + rxt(258)*y(7) + rxt(255)*y(86) &
                      + rxt(259)*y(23) + rxt(261)*y(24) + rxt(260)*y(30)
         mat(53) = (.500_r8*rxt(264)+rxt(265))*y(86) + rxt(266)*y(8)
      end subroutine nlnmat03
      subroutine nlnmat04( mat, y, rxt )
      use chem_mods, only : gas_pcnst, rxntot, nzcnt
      implicit none
!----------------------------------------------
! ... dummy arguments
!----------------------------------------------
      real(r8), intent(in) :: y(gas_pcnst)
      real(r8), intent(in) :: rxt(rxntot)
      real(r8), intent(inout) :: mat(nzcnt)
!----------------------------------------------
! ... local variables
!----------------------------------------------
!----------------------------------------------
! ... complete matrix entries implicit species
!----------------------------------------------
         mat(48) = -(rxt(263)*y(98))
         mat(549) = -rxt(263)*y(61)
         mat(487) = rxt(262)*y(60)
         mat(226) = rxt(262)*y(86)
         mat(547) = rxt(263)*y(61)
         mat(47) = rxt(263)*y(98)
         mat(52) = -((rxt(264) + rxt(265)) * y(86) + rxt(266)*y(8))
         mat(488) = -(rxt(264) + rxt(265)) * y(63)
         mat(522) = -rxt(266)*y(63)
      end subroutine nlnmat04
      subroutine nlnmat_finit( mat, lmat, dti )
      use chem_mods, only : gas_pcnst, rxntot, nzcnt
      implicit none
!----------------------------------------------
! ... dummy arguments
!----------------------------------------------
      real(r8), intent(in) :: dti
      real(r8), intent(in) :: lmat(nzcnt)
      real(r8), intent(inout) :: mat(nzcnt)
!----------------------------------------------
! ... local variables
!----------------------------------------------
!----------------------------------------------
! ... complete matrix entries implicit species
!----------------------------------------------
         mat( 1) = lmat( 1)
         mat( 2) = lmat( 2)
         mat( 3) = lmat( 3)
         mat( 4) = lmat( 4)
         mat( 5) = lmat( 5)
         mat( 6) = lmat( 6)
         mat( 7) = lmat( 7)
         mat( 8) = lmat( 8)
         mat( 9) = lmat( 9)
         mat( 10) = lmat( 10)
         mat( 11) = lmat( 11)
         mat( 12) = lmat( 12)
         mat( 13) = lmat( 13)
         mat( 14) = lmat( 14)
         mat( 15) = lmat( 15)
         mat( 16) = lmat( 16)
         mat( 17) = lmat( 17)
         mat( 18) = lmat( 18)
         mat( 19) = lmat( 19)
         mat( 20) = lmat( 20)
         mat( 21) = lmat( 21)
         mat( 22) = lmat( 22)
         mat( 23) = lmat( 23)
         mat( 24) = mat( 24) + lmat( 24)
         mat( 25) = mat( 25) + lmat( 25)
         mat( 27) = lmat( 27)
         mat( 28) = lmat( 28)
         mat( 29) = lmat( 29)
         mat( 30) = mat( 30) + lmat( 30)
         mat( 31) = mat( 31) + lmat( 31)
         mat( 32) = mat( 32) + lmat( 32)
         mat( 33) = mat( 33) + lmat( 33)
         mat( 34) = lmat( 34)
         mat( 35) = lmat( 35)
         mat( 36) = lmat( 36)
         mat( 37) = mat( 37) + lmat( 37)
         mat( 38) = mat( 38) + lmat( 38)
         mat( 39) = mat( 39) + lmat( 39)
         mat( 41) = lmat( 41)
         mat( 42) = lmat( 42)
         mat( 43) = lmat( 43)
         mat( 44) = lmat( 44)
         mat( 45) = lmat( 45)
         mat( 46) = lmat( 46)
         mat( 48) = mat( 48) + lmat( 48)
         mat( 49) = lmat( 49)
         mat( 51) = lmat( 51)
         mat( 52) = mat( 52) + lmat( 52)
         mat( 58) = lmat( 58)
         mat( 59) = lmat( 59)
         mat( 60) = lmat( 60)
         mat( 61) = lmat( 61)
         mat( 62) = lmat( 62)
         mat( 63) = lmat( 63)
         mat( 64) = mat( 64) + lmat( 64)
         mat( 67) = mat( 67) + lmat( 67)
         mat( 68) = lmat( 68)
         mat( 69) = mat( 69) + lmat( 69)
         mat( 70) = mat( 70) + lmat( 70)
         mat( 71) = lmat( 71)
         mat( 77) = mat( 77) + lmat( 77)
         mat( 78) = lmat( 78)
         mat( 79) = lmat( 79)
         mat( 80) = mat( 80) + lmat( 80)
         mat( 83) = mat( 83) + lmat( 83)
         mat( 90) = mat( 90) + lmat( 90)
         mat( 91) = mat( 91) + lmat( 91)
         mat( 92) = mat( 92) + lmat( 92)
         mat( 93) = lmat( 93)
         mat( 96) = lmat( 96)
         mat( 97) = mat( 97) + lmat( 97)
         mat( 102) = mat( 102) + lmat( 102)
         mat( 106) = mat( 106) + lmat( 106)
         mat( 107) = lmat( 107)
         mat( 108) = lmat( 108)
         mat( 109) = mat( 109) + lmat( 109)
         mat( 110) = lmat( 110)
         mat( 111) = lmat( 111)
         mat( 112) = mat( 112) + lmat( 112)
         mat( 114) = mat( 114) + lmat( 114)
         mat( 116) = mat( 116) + lmat( 116)
         mat( 122) = mat( 122) + lmat( 122)
         mat( 132) = mat( 132) + lmat( 132)
         mat( 139) = lmat( 139)
         mat( 140) = mat( 140) + lmat( 140)
         mat( 141) = mat( 141) + lmat( 141)
         mat( 142) = lmat( 142)
         mat( 143) = lmat( 143)
         mat( 149) = mat( 149) + lmat( 149)
         mat( 150) = mat( 150) + lmat( 150)
         mat( 157) = lmat( 157)
         mat( 158) = mat( 158) + lmat( 158)
         mat( 159) = mat( 159) + lmat( 159)
         mat( 162) = mat( 162) + lmat( 162)
         mat( 167) = mat( 167) + lmat( 167)
         mat( 179) = mat( 179) + lmat( 179)
         mat( 189) = mat( 189) + lmat( 189)
         mat( 192) = lmat( 192)
         mat( 193) = mat( 193) + lmat( 193)
         mat( 197) = mat( 197) + lmat( 197)
         mat( 198) = lmat( 198)
         mat( 199) = mat( 199) + lmat( 199)
         mat( 204) = lmat( 204)
         mat( 208) = lmat( 208)
         mat( 209) = mat( 209) + lmat( 209)
         mat( 219) = mat( 219) + lmat( 219)
         mat( 220) = mat( 220) + lmat( 220)
         mat( 222) = mat( 222) + lmat( 222)
         mat( 227) = mat( 227) + lmat( 227)
         mat( 228) = lmat( 228)
         mat( 232) = lmat( 232)
         mat( 234) = mat( 234) + lmat( 234)
         mat( 235) = mat( 235) + lmat( 235)
         mat( 236) = mat( 236) + lmat( 236)
         mat( 237) = mat( 237) + lmat( 237)
         mat( 238) = mat( 238) + lmat( 238)
         mat( 239) = lmat( 239)
         mat( 241) = mat( 241) + lmat( 241)
         mat( 247) = mat( 247) + lmat( 247)
         mat( 250) = lmat( 250)
         mat( 251) = mat( 251) + lmat( 251)
         mat( 256) = lmat( 256)
         mat( 258) = mat( 258) + lmat( 258)
         mat( 271) = mat( 271) + lmat( 271)
         mat( 272) = mat( 272) + lmat( 272)
         mat( 280) = mat( 280) + lmat( 280)
         mat( 283) = mat( 283) + lmat( 283)
         mat( 284) = mat( 284) + lmat( 284)
         mat( 285) = lmat( 285)
         mat( 307) = mat( 307) + lmat( 307)
         mat( 311) = mat( 311) + lmat( 311)
         mat( 322) = mat( 322) + lmat( 322)
         mat( 327) = mat( 327) + lmat( 327)
         mat( 341) = mat( 341) + lmat( 341)
         mat( 349) = mat( 349) + lmat( 349)
         mat( 351) = lmat( 351)
         mat( 358) = mat( 358) + lmat( 358)
         mat( 359) = lmat( 359)
         mat( 360) = mat( 360) + lmat( 360)
         mat( 376) = mat( 376) + lmat( 376)
         mat( 377) = mat( 377) + lmat( 377)
         mat( 390) = mat( 390) + lmat( 390)
         mat( 399) = mat( 399) + lmat( 399)
         mat( 412) = mat( 412) + lmat( 412)
         mat( 423) = mat( 423) + lmat( 423)
         mat( 429) = mat( 429) + lmat( 429)
         mat( 433) = mat( 433) + lmat( 433)
         mat( 434) = mat( 434) + lmat( 434)
         mat( 436) = mat( 436) + lmat( 436)
         mat( 437) = lmat( 437)
         mat( 438) = mat( 438) + lmat( 438)
         mat( 447) = mat( 447) + lmat( 447)
         mat( 457) = mat( 457) + lmat( 457)
         mat( 465) = mat( 465) + lmat( 465)
         mat( 474) = mat( 474) + lmat( 474)
         mat( 475) = mat( 475) + lmat( 475)
         mat( 476) = mat( 476) + lmat( 476)
         mat( 484) = mat( 484) + lmat( 484)
         mat( 485) = lmat( 485)
         mat( 486) = lmat( 486)
         mat( 507) = mat( 507) + lmat( 507)
         mat( 509) = mat( 509) + lmat( 509)
         mat( 513) = mat( 513) + lmat( 513)
         mat( 515) = mat( 515) + lmat( 515)
         mat( 517) = mat( 517) + lmat( 517)
         mat( 518) = mat( 518) + lmat( 518)
         mat( 526) = mat( 526) + lmat( 526)
         mat( 536) = mat( 536) + lmat( 536)
         mat( 537) = mat( 537) + lmat( 537)
         mat( 539) = mat( 539) + lmat( 539)
         mat( 541) = mat( 541) + lmat( 541)
         mat( 546) = mat( 546) + lmat( 546)
         mat( 554) = lmat( 554)
         mat( 556) = lmat( 556)
         mat( 563) = mat( 563) + lmat( 563)
         mat( 565) = mat( 565) + lmat( 565)
         mat( 570) = mat( 570) + lmat( 570)
         mat( 571) = lmat( 571)
         mat( 578) = mat( 578) + lmat( 578)
         mat( 579) = mat( 579) + lmat( 579)
         mat( 581) = lmat( 581)
         mat( 597) = mat( 597) + lmat( 597)
         mat( 601) = mat( 601) + lmat( 601)
         mat( 602) = mat( 602) + lmat( 602)
         mat( 605) = mat( 605) + lmat( 605)
         mat( 623) = mat( 623) + lmat( 623)
         mat( 643) = mat( 643) + lmat( 643)
         mat( 653) = lmat( 653)
         mat( 656) = mat( 656) + lmat( 656)
         mat( 668) = mat( 668) + lmat( 668)
         mat( 671) = mat( 671) + lmat( 671)
         mat( 673) = mat( 673) + lmat( 673)
         mat( 674) = mat( 674) + lmat( 674)
         mat( 677) = mat( 677) + lmat( 677)
         mat( 678) = lmat( 678)
         mat( 679) = mat( 679) + lmat( 679)
         mat( 680) = mat( 680) + lmat( 680)
         mat( 682) = mat( 682) + lmat( 682)
         mat( 684) = lmat( 684)
         mat( 686) = mat( 686) + lmat( 686)
         mat( 689) = mat( 689) + lmat( 689)
         mat( 690) = lmat( 690)
         mat( 691) = lmat( 691)
         mat( 693) = mat( 693) + lmat( 693)
         mat( 694) = mat( 694) + lmat( 694)
         mat( 705) = mat( 705) + lmat( 705)
         mat( 708) = lmat( 708)
         mat( 735) = mat( 735) + lmat( 735)
         mat( 128) = 0._r8
         mat( 130) = 0._r8
         mat( 144) = 0._r8
         mat( 145) = 0._r8
         mat( 147) = 0._r8
         mat( 170) = 0._r8
         mat( 171) = 0._r8
         mat( 172) = 0._r8
         mat( 176) = 0._r8
         mat( 181) = 0._r8
         mat( 182) = 0._r8
         mat( 183) = 0._r8
         mat( 184) = 0._r8
         mat( 185) = 0._r8
         mat( 191) = 0._r8
         mat( 205) = 0._r8
         mat( 215) = 0._r8
         mat( 216) = 0._r8
         mat( 230) = 0._r8
         mat( 242) = 0._r8
         mat( 248) = 0._r8
         mat( 249) = 0._r8
         mat( 254) = 0._r8
         mat( 268) = 0._r8
         mat( 270) = 0._r8
         mat( 288) = 0._r8
         mat( 289) = 0._r8
         mat( 293) = 0._r8
         mat( 295) = 0._r8
         mat( 297) = 0._r8
         mat( 304) = 0._r8
         mat( 306) = 0._r8
         mat( 308) = 0._r8
         mat( 312) = 0._r8
         mat( 316) = 0._r8
         mat( 317) = 0._r8
         mat( 320) = 0._r8
         mat( 321) = 0._r8
         mat( 325) = 0._r8
         mat( 329) = 0._r8
         mat( 330) = 0._r8
         mat( 332) = 0._r8
         mat( 334) = 0._r8
         mat( 335) = 0._r8
         mat( 336) = 0._r8
         mat( 345) = 0._r8
         mat( 353) = 0._r8
         mat( 356) = 0._r8
         mat( 361) = 0._r8
         mat( 362) = 0._r8
         mat( 370) = 0._r8
         mat( 372) = 0._r8
         mat( 375) = 0._r8
         mat( 379) = 0._r8
         mat( 383) = 0._r8
         mat( 384) = 0._r8
         mat( 389) = 0._r8
         mat( 393) = 0._r8
         mat( 396) = 0._r8
         mat( 397) = 0._r8
         mat( 398) = 0._r8
         mat( 401) = 0._r8
         mat( 402) = 0._r8
         mat( 403) = 0._r8
         mat( 404) = 0._r8
         mat( 405) = 0._r8
         mat( 408) = 0._r8
         mat( 409) = 0._r8
         mat( 410) = 0._r8
         mat( 411) = 0._r8
         mat( 419) = 0._r8
         mat( 428) = 0._r8
         mat( 431) = 0._r8
         mat( 432) = 0._r8
         mat( 439) = 0._r8
         mat( 442) = 0._r8
         mat( 451) = 0._r8
         mat( 455) = 0._r8
         mat( 456) = 0._r8
         mat( 461) = 0._r8
         mat( 467) = 0._r8
         mat( 469) = 0._r8
         mat( 470) = 0._r8
         mat( 472) = 0._r8
         mat( 478) = 0._r8
         mat( 481) = 0._r8
         mat( 482) = 0._r8
         mat( 483) = 0._r8
         mat( 495) = 0._r8
         mat( 520) = 0._r8
         mat( 527) = 0._r8
         mat( 528) = 0._r8
         mat( 530) = 0._r8
         mat( 531) = 0._r8
         mat( 532) = 0._r8
         mat( 533) = 0._r8
         mat( 534) = 0._r8
         mat( 535) = 0._r8
         mat( 540) = 0._r8
         mat( 543) = 0._r8
         mat( 544) = 0._r8
         mat( 545) = 0._r8
         mat( 551) = 0._r8
         mat( 552) = 0._r8
         mat( 553) = 0._r8
         mat( 555) = 0._r8
         mat( 557) = 0._r8
         mat( 558) = 0._r8
         mat( 559) = 0._r8
         mat( 560) = 0._r8
         mat( 561) = 0._r8
         mat( 562) = 0._r8
         mat( 564) = 0._r8
         mat( 566) = 0._r8
         mat( 567) = 0._r8
         mat( 568) = 0._r8
         mat( 569) = 0._r8
         mat( 586) = 0._r8
         mat( 587) = 0._r8
         mat( 588) = 0._r8
         mat( 589) = 0._r8
         mat( 590) = 0._r8
         mat( 593) = 0._r8
         mat( 594) = 0._r8
         mat( 595) = 0._r8
         mat( 596) = 0._r8
         mat( 599) = 0._r8
         mat( 600) = 0._r8
         mat( 610) = 0._r8
         mat( 626) = 0._r8
         mat( 629) = 0._r8
         mat( 631) = 0._r8
         mat( 634) = 0._r8
         mat( 635) = 0._r8
         mat( 638) = 0._r8
         mat( 639) = 0._r8
         mat( 640) = 0._r8
         mat( 644) = 0._r8
         mat( 645) = 0._r8
         mat( 646) = 0._r8
         mat( 654) = 0._r8
         mat( 655) = 0._r8
         mat( 657) = 0._r8
         mat( 658) = 0._r8
         mat( 659) = 0._r8
         mat( 660) = 0._r8
         mat( 661) = 0._r8
         mat( 663) = 0._r8
         mat( 665) = 0._r8
         mat( 666) = 0._r8
         mat( 667) = 0._r8
         mat( 676) = 0._r8
         mat( 681) = 0._r8
         mat( 685) = 0._r8
         mat( 687) = 0._r8
         mat( 704) = 0._r8
         mat( 706) = 0._r8
         mat( 713) = 0._r8
         mat( 715) = 0._r8
         mat( 729) = 0._r8
         mat( 732) = 0._r8
         mat( 734) = 0._r8
         mat( 1) = mat( 1) - dti
         mat( 2) = mat( 2) - dti
         mat( 3) = mat( 3) - dti
         mat( 4) = mat( 4) - dti
         mat( 5) = mat( 5) - dti
         mat( 6) = mat( 6) - dti
         mat( 7) = mat( 7) - dti
         mat( 8) = mat( 8) - dti
         mat( 9) = mat( 9) - dti
         mat( 10) = mat( 10) - dti
         mat( 11) = mat( 11) - dti
         mat( 12) = mat( 12) - dti
         mat( 13) = mat( 13) - dti
         mat( 14) = mat( 14) - dti
         mat( 15) = mat( 15) - dti
         mat( 16) = mat( 16) - dti
         mat( 17) = mat( 17) - dti
         mat( 18) = mat( 18) - dti
         mat( 19) = mat( 19) - dti
         mat( 20) = mat( 20) - dti
         mat( 21) = mat( 21) - dti
         mat( 24) = mat( 24) - dti
         mat( 27) = mat( 27) - dti
         mat( 30) = mat( 30) - dti
         mat( 33) = mat( 33) - dti
         mat( 35) = mat( 35) - dti
         mat( 37) = mat( 37) - dti
         mat( 41) = mat( 41) - dti
         mat( 44) = mat( 44) - dti
         mat( 48) = mat( 48) - dti
         mat( 52) = mat( 52) - dti
         mat( 58) = mat( 58) - dti
         mat( 64) = mat( 64) - dti
         mat( 69) = mat( 69) - dti
         mat( 77) = mat( 77) - dti
         mat( 83) = mat( 83) - dti
         mat( 90) = mat( 90) - dti
         mat( 97) = mat( 97) - dti
         mat( 102) = mat( 102) - dti
         mat( 106) = mat( 106) - dti
         mat( 114) = mat( 114) - dti
         mat( 122) = mat( 122) - dti
         mat( 132) = mat( 132) - dti
         mat( 140) = mat( 140) - dti
         mat( 150) = mat( 150) - dti
         mat( 158) = mat( 158) - dti
         mat( 167) = mat( 167) - dti
         mat( 179) = mat( 179) - dti
         mat( 189) = mat( 189) - dti
         mat( 197) = mat( 197) - dti
         mat( 209) = mat( 209) - dti
         mat( 219) = mat( 219) - dti
         mat( 227) = mat( 227) - dti
         mat( 235) = mat( 235) - dti
         mat( 247) = mat( 247) - dti
         mat( 258) = mat( 258) - dti
         mat( 272) = mat( 272) - dti
         mat( 284) = mat( 284) - dti
         mat( 307) = mat( 307) - dti
         mat( 327) = mat( 327) - dti
         mat( 349) = mat( 349) - dti
         mat( 377) = mat( 377) - dti
         mat( 399) = mat( 399) - dti
         mat( 423) = mat( 423) - dti
         mat( 447) = mat( 447) - dti
         mat( 475) = mat( 475) - dti
         mat( 513) = mat( 513) - dti
         mat( 539) = mat( 539) - dti
         mat( 565) = mat( 565) - dti
         mat( 597) = mat( 597) - dti
         mat( 623) = mat( 623) - dti
         mat( 643) = mat( 643) - dti
         mat( 668) = mat( 668) - dti
         mat( 693) = mat( 693) - dti
         mat( 735) = mat( 735) - dti
      end subroutine nlnmat_finit
      subroutine nlnmat( mat, y, rxt, lmat, dti )
      use chem_mods, only : gas_pcnst, rxntot, nzcnt
      implicit none
!----------------------------------------------
! ... dummy arguments
!----------------------------------------------
      real(r8), intent(in) :: dti
      real(r8), intent(in) :: lmat(nzcnt)
      real(r8), intent(in) :: y(gas_pcnst)
      real(r8), intent(in) :: rxt(rxntot)
      real(r8), intent(inout) :: mat(nzcnt)
      call nlnmat01( mat, y, rxt )
      call nlnmat02( mat, y, rxt )
      call nlnmat03( mat, y, rxt )
      call nlnmat04( mat, y, rxt )
      call nlnmat_finit( mat, lmat, dti )
      end subroutine nlnmat
      end module mo_nln_matrix
