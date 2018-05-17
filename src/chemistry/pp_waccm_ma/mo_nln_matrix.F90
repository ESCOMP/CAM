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
         mat(610) = -(rxt(89)*y(2) + rxt(107)*y(3) + rxt(133)*y(21) + rxt(138)*y(22) &
                      + rxt(146)*y(23) + rxt(159)*y(9) + rxt(162)*y(10) + rxt(174) &
                      *y(28) + rxt(201)*y(37))
         mat(504) = -rxt(89)*y(1)
         mat(570) = -rxt(107)*y(1)
         mat(339) = -rxt(133)*y(1)
         mat(465) = -rxt(138)*y(1)
         mat(410) = -rxt(146)*y(1)
         mat(592) = -rxt(159)*y(1)
         mat(385) = -rxt(162)*y(1)
         mat(434) = -rxt(174)*y(1)
         mat(360) = -rxt(201)*y(1)
         mat(504) = mat(504) + rxt(88)*y(4)
         mat(233) = rxt(88)*y(2)
         mat(499) = -(rxt(88)*y(4) + rxt(89)*y(1) + 4._r8*rxt(90)*y(2) + rxt(137) &
                      *y(22) + rxt(144)*y(20) + rxt(145)*y(23) + rxt(148)*y(24) &
                      + rxt(157)*y(9) + (rxt(160) + rxt(161)) * y(10) + rxt(168)*y(11) &
                      + rxt(181)*y(30) + rxt(194)*y(33) + rxt(195)*y(34) + rxt(198) &
                      *y(35) + rxt(204)*y(38) + rxt(214)*y(39) + rxt(215)*y(40) &
                      + rxt(216)*y(41) + rxt(238)*y(18) + (rxt(265) + rxt(266) &
                      ) * y(65) + rxt(272)*y(67))
         mat(230) = -rxt(88)*y(2)
         mat(605) = -rxt(89)*y(2)
         mat(460) = -rxt(137)*y(2)
         mat(272) = -rxt(144)*y(2)
         mat(405) = -rxt(145)*y(2)
         mat(81) = -rxt(148)*y(2)
         mat(587) = -rxt(157)*y(2)
         mat(380) = -(rxt(160) + rxt(161)) * y(2)
         mat(519) = -rxt(168)*y(2)
         mat(319) = -rxt(181)*y(2)
         mat(542) = -rxt(194)*y(2)
         mat(189) = -rxt(195)*y(2)
         mat(200) = -rxt(198)*y(2)
         mat(293) = -rxt(204)*y(2)
         mat(167) = -rxt(214)*y(2)
         mat(118) = -rxt(215)*y(2)
         mat(74) = -rxt(216)*y(2)
         mat(258) = -rxt(238)*y(2)
         mat(98) = -(rxt(265) + rxt(266)) * y(2)
         mat(89) = -rxt(272)*y(2)
         mat(565) = (rxt(102)+rxt(103))*y(4)
         mat(230) = mat(230) + (rxt(102)+rxt(103))*y(3) + rxt(152)*y(8) + rxt(271) &
                      *y(67) + rxt(263)*y(68) + rxt(284)*y(73)
         mat(180) = rxt(152)*y(4) + rxt(153)*y(9) + rxt(154)*y(10) + rxt(268)*y(66)
         mat(587) = mat(587) + rxt(153)*y(8)
         mat(380) = mat(380) + rxt(154)*y(8)
         mat(460) = mat(460) + 2.000_r8*rxt(140)*y(22)
         mat(334) = rxt(136)*y(23)
         mat(405) = mat(405) + rxt(136)*y(21)
         mat(125) = rxt(268)*y(8) + 1.150_r8*rxt(276)*y(70)
         mat(89) = mat(89) + rxt(271)*y(4)
         mat(110) = rxt(263)*y(4)
         mat(133) = rxt(275)*y(70)
         mat(150) = 1.150_r8*rxt(276)*y(66) + rxt(275)*y(69)
         mat(58) = rxt(284)*y(4)
         mat(568) = -((rxt(102) + rxt(103)) * y(4) + rxt(104)*y(74) + rxt(107)*y(1) &
                      + rxt(124)*y(60) + rxt(125)*y(61) + rxt(129)*y(20) + rxt(130) &
                      *y(33) + rxt(131)*y(39))
         mat(231) = -(rxt(102) + rxt(103)) * y(3)
         mat(244) = -rxt(104)*y(3)
         mat(608) = -rxt(107)*y(3)
         mat(14) = -rxt(124)*y(3)
         mat(20) = -rxt(125)*y(3)
         mat(275) = -rxt(129)*y(3)
         mat(545) = -rxt(130)*y(3)
         mat(168) = -rxt(131)*y(3)
         mat(231) = mat(231) + rxt(149)*y(71)
         mat(126) = .850_r8*rxt(276)*y(70)
         mat(62) = rxt(149)*y(4)
         mat(151) = .850_r8*rxt(276)*y(66)
         mat(225) = -(rxt(88)*y(2) + rxt(98)*y(6) + rxt(102)*y(3) + rxt(132)*y(21) &
                      + rxt(149)*y(71) + rxt(152)*y(8) + rxt(263)*y(68) + (rxt(270) &
                      + rxt(271)) * y(67) + rxt(273)*y(65) + rxt(284)*y(73))
         mat(487) = -rxt(88)*y(4)
         mat(5) = -rxt(98)*y(4)
         mat(555) = -rxt(102)*y(4)
         mat(326) = -rxt(132)*y(4)
         mat(61) = -rxt(149)*y(4)
         mat(175) = -rxt(152)*y(4)
         mat(106) = -rxt(263)*y(4)
         mat(88) = -(rxt(270) + rxt(271)) * y(4)
         mat(97) = -rxt(273)*y(4)
         mat(57) = -rxt(284)*y(4)
         mat(596) = 2.000_r8*rxt(89)*y(2) + 2.000_r8*rxt(107)*y(3) + rxt(159)*y(9) &
                      + rxt(162)*y(10) + rxt(138)*y(22) + rxt(133)*y(21) &
                      + 2.000_r8*rxt(146)*y(23) + rxt(174)*y(28) + rxt(201)*y(37)
         mat(487) = mat(487) + 2.000_r8*rxt(89)*y(1) + 2.000_r8*rxt(90)*y(2) + rxt(97) &
                      *y(6) + rxt(160)*y(10) + rxt(137)*y(22) + rxt(168)*y(11) &
                      + rxt(145)*y(23) + rxt(181)*y(30) + rxt(204)*y(38)
         mat(555) = mat(555) + 2.000_r8*rxt(107)*y(1)
         mat(225) = mat(225) + 2.000_r8*rxt(98)*y(6)
         mat(5) = mat(5) + rxt(97)*y(2) + 2.000_r8*rxt(98)*y(4)
         mat(175) = mat(175) + rxt(156)*y(10)
         mat(576) = rxt(159)*y(1) + rxt(269)*y(66)
         mat(369) = rxt(162)*y(1) + rxt(160)*y(2) + rxt(156)*y(8)
         mat(448) = rxt(138)*y(1) + rxt(137)*y(2) + rxt(172)*y(13) + rxt(139)*y(23) &
                      + rxt(183)*y(30)
         mat(509) = rxt(168)*y(2) + rxt(170)*y(23)
         mat(40) = rxt(172)*y(22)
         mat(613) = rxt(240)*y(23)
         mat(326) = mat(326) + rxt(133)*y(1) + rxt(135)*y(23)
         mat(393) = 2.000_r8*rxt(146)*y(1) + rxt(145)*y(2) + rxt(139)*y(22) + rxt(170) &
                      *y(11) + rxt(240)*y(16) + rxt(135)*y(21) + 2.000_r8*rxt(147) &
                      *y(23) + rxt(177)*y(28) + rxt(184)*y(30) + rxt(202)*y(37) &
                      + rxt(206)*y(38)
         mat(418) = rxt(174)*y(1) + rxt(177)*y(23)
         mat(307) = rxt(181)*y(2) + rxt(183)*y(22) + rxt(184)*y(23) + ( &
                      + 2.000_r8*rxt(188)+2.000_r8*rxt(189))*y(30) + (rxt(210) &
                       +rxt(211))*y(38)
         mat(343) = rxt(201)*y(1) + rxt(202)*y(23)
         mat(282) = rxt(204)*y(2) + rxt(206)*y(23) + (rxt(210)+rxt(211))*y(30) &
                      + 2.000_r8*rxt(212)*y(38)
         mat(124) = rxt(269)*y(9)
         mat(7) = -(rxt(91)*y(2) + rxt(92)*y(4) + rxt(94)*y(1))
         mat(468) = -rxt(91)*y(5)
         mat(214) = -rxt(92)*y(5)
         mat(595) = -rxt(94)*y(5)
         mat(549) = rxt(102)*y(4)
         mat(214) = mat(214) + rxt(102)*y(3)
         mat(4) = -(rxt(97)*y(2) + rxt(98)*y(4))
         mat(467) = -rxt(97)*y(6)
         mat(213) = -rxt(98)*y(6)
         mat(594) = rxt(94)*y(5)
         mat(467) = mat(467) + rxt(91)*y(5)
         mat(213) = mat(213) + rxt(92)*y(5)
         mat(6) = rxt(94)*y(1) + rxt(91)*y(2) + rxt(92)*y(4)
         mat(267) = -(rxt(129)*y(3) + rxt(142)*y(22) + rxt(144)*y(2) + rxt(175)*y(28) &
                      + rxt(218)*y(63))
         mat(558) = -rxt(129)*y(20)
         mat(451) = -rxt(142)*y(20)
         mat(490) = -rxt(144)*y(20)
         mat(421) = -rxt(175)*y(20)
         mat(157) = -rxt(218)*y(20)
         mat(328) = rxt(135)*y(23)
         mat(396) = rxt(135)*y(21)
         mat(64) = -((rxt(234) + rxt(235)) * y(22))
         mat(440) = -(rxt(234) + rxt(235)) * y(19)
         mat(472) = rxt(238)*y(18)
         mat(440) = mat(440) + rxt(237)*y(18)
         mat(507) = rxt(236)*y(18)
         mat(246) = rxt(238)*y(2) + rxt(237)*y(22) + rxt(236)*y(11) + rxt(179)*y(28) &
                      + rxt(203)*y(37)
         mat(413) = rxt(179)*y(18)
         mat(341) = rxt(203)*y(18)
         mat(174) = -(rxt(151)*y(22) + rxt(152)*y(4) + rxt(153)*y(9) + (rxt(154) &
                      + rxt(155) + rxt(156)) * y(10) + rxt(268)*y(66))
         mat(444) = -rxt(151)*y(8)
         mat(224) = -rxt(152)*y(8)
         mat(575) = -rxt(153)*y(8)
         mat(366) = -(rxt(154) + rxt(155) + rxt(156)) * y(8)
         mat(123) = -rxt(268)*y(8)
         mat(483) = rxt(272)*y(67) + rxt(150)*y(71)
         mat(224) = mat(224) + rxt(270)*y(67)
         mat(96) = 1.100_r8*rxt(277)*y(70)
         mat(87) = rxt(272)*y(2) + rxt(270)*y(4)
         mat(131) = .200_r8*rxt(275)*y(70)
         mat(60) = rxt(150)*y(2)
         mat(145) = 1.100_r8*rxt(277)*y(65) + .200_r8*rxt(275)*y(69)
         mat(591) = -(rxt(153)*y(8) + rxt(157)*y(2) + rxt(158)*y(23) + rxt(159)*y(1) &
                      + rxt(167)*y(11) + rxt(186)*y(30) + rxt(207)*y(38) + rxt(239) &
                      *y(16) + rxt(269)*y(66))
         mat(182) = -rxt(153)*y(9)
         mat(503) = -rxt(157)*y(9)
         mat(409) = -rxt(158)*y(9)
         mat(609) = -rxt(159)*y(9)
         mat(523) = -rxt(167)*y(9)
         mat(323) = -rxt(186)*y(9)
         mat(297) = -rxt(207)*y(9)
         mat(628) = -rxt(239)*y(9)
         mat(127) = -rxt(269)*y(9)
         mat(503) = mat(503) + rxt(160)*y(10)
         mat(232) = rxt(152)*y(8) + rxt(149)*y(71)
         mat(182) = mat(182) + rxt(152)*y(4) + 2.000_r8*rxt(155)*y(10) + rxt(151) &
                      *y(22)
         mat(384) = rxt(160)*y(2) + 2.000_r8*rxt(155)*y(8)
         mat(464) = rxt(151)*y(8)
         mat(63) = rxt(149)*y(4)
         mat(376) = -((rxt(154) + rxt(155) + rxt(156)) * y(8) + (rxt(160) + rxt(161) &
                      ) * y(2) + rxt(162)*y(1) + rxt(163)*y(11) + rxt(165)*y(22) &
                      + rxt(171)*y(23) + rxt(187)*y(30) + rxt(208)*y(38))
         mat(177) = -(rxt(154) + rxt(155) + rxt(156)) * y(10)
         mat(495) = -(rxt(160) + rxt(161)) * y(10)
         mat(601) = -rxt(162)*y(10)
         mat(515) = -rxt(163)*y(10)
         mat(456) = -rxt(165)*y(10)
         mat(401) = -rxt(171)*y(10)
         mat(315) = -rxt(187)*y(10)
         mat(289) = -rxt(208)*y(10)
         mat(601) = mat(601) + rxt(159)*y(9)
         mat(495) = mat(495) + rxt(157)*y(9) + rxt(168)*y(11)
         mat(583) = rxt(159)*y(1) + rxt(157)*y(2) + 2.000_r8*rxt(167)*y(11) + rxt(239) &
                      *y(16) + rxt(158)*y(23) + rxt(186)*y(30) + rxt(207)*y(38)
         mat(456) = mat(456) + rxt(169)*y(11) + rxt(172)*y(13)
         mat(515) = mat(515) + rxt(168)*y(2) + 2.000_r8*rxt(167)*y(9) + rxt(169)*y(22) &
                      + rxt(170)*y(23)
         mat(42) = rxt(172)*y(22)
         mat(620) = rxt(239)*y(9)
         mat(401) = mat(401) + rxt(158)*y(9) + rxt(170)*y(11)
         mat(315) = mat(315) + rxt(186)*y(9)
         mat(289) = mat(289) + rxt(207)*y(9)
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
         mat(459) = -(rxt(137)*y(2) + rxt(138)*y(1) + rxt(139)*y(23) + (4._r8*rxt(140) &
                      + 4._r8*rxt(141)) * y(22) + rxt(142)*y(20) + rxt(143)*y(24) &
                      + rxt(151)*y(8) + rxt(165)*y(10) + rxt(166)*y(12) + rxt(169) &
                      *y(11) + rxt(172)*y(13) + (rxt(182) + rxt(183)) * y(30) + rxt(193) &
                      *y(33) + rxt(197)*y(34) + rxt(199)*y(35) + rxt(205)*y(38) &
                      + rxt(213)*y(39) + (rxt(234) + rxt(235)) * y(19) + rxt(237) &
                      *y(18) + rxt(241)*y(17))
         mat(498) = -rxt(137)*y(22)
         mat(604) = -rxt(138)*y(22)
         mat(404) = -rxt(139)*y(22)
         mat(271) = -rxt(142)*y(22)
         mat(80) = -rxt(143)*y(22)
         mat(179) = -rxt(151)*y(22)
         mat(379) = -rxt(165)*y(22)
         mat(210) = -rxt(166)*y(22)
         mat(518) = -rxt(169)*y(22)
         mat(44) = -rxt(172)*y(22)
         mat(318) = -(rxt(182) + rxt(183)) * y(22)
         mat(541) = -rxt(193)*y(22)
         mat(188) = -rxt(197)*y(22)
         mat(199) = -rxt(199)*y(22)
         mat(292) = -rxt(205)*y(22)
         mat(166) = -rxt(213)*y(22)
         mat(67) = -(rxt(234) + rxt(235)) * y(22)
         mat(257) = -rxt(237)*y(22)
         mat(37) = -rxt(241)*y(22)
         mat(604) = mat(604) + rxt(133)*y(21) + rxt(146)*y(23)
         mat(498) = mat(498) + rxt(144)*y(20) + rxt(238)*y(18) + rxt(145)*y(23) &
                      + rxt(148)*y(24) + rxt(194)*y(33) + rxt(195)*y(34) + rxt(214) &
                      *y(39) + rxt(215)*y(40)
         mat(564) = rxt(129)*y(20) + 2.000_r8*rxt(104)*y(74) + rxt(130)*y(33) &
                      + rxt(131)*y(39)
         mat(271) = mat(271) + rxt(144)*y(2) + rxt(129)*y(3)
         mat(586) = rxt(158)*y(23)
         mat(518) = mat(518) + rxt(170)*y(23)
         mat(257) = mat(257) + rxt(238)*y(2)
         mat(333) = rxt(133)*y(1) + 2.000_r8*rxt(134)*y(23)
         mat(404) = mat(404) + rxt(146)*y(1) + rxt(145)*y(2) + rxt(158)*y(9) &
                      + rxt(170)*y(11) + 2.000_r8*rxt(134)*y(21) + rxt(178)*y(28)
         mat(80) = mat(80) + rxt(148)*y(2)
         mat(241) = 2.000_r8*rxt(104)*y(3) + rxt(217)*y(63)
         mat(428) = rxt(178)*y(23)
         mat(541) = mat(541) + rxt(194)*y(2) + rxt(130)*y(3)
         mat(188) = mat(188) + rxt(195)*y(2)
         mat(166) = mat(166) + rxt(214)*y(2) + rxt(131)*y(3)
         mat(117) = rxt(215)*y(2)
         mat(159) = rxt(217)*y(74)
         mat(520) = -(rxt(163)*y(10) + rxt(167)*y(9) + rxt(168)*y(2) + rxt(169)*y(22) &
                      + rxt(170)*y(23) + rxt(236)*y(18))
         mat(381) = -rxt(163)*y(11)
         mat(588) = -rxt(167)*y(11)
         mat(500) = -rxt(168)*y(11)
         mat(461) = -rxt(169)*y(11)
         mat(406) = -rxt(170)*y(11)
         mat(259) = -rxt(236)*y(11)
         mat(606) = rxt(162)*y(10)
         mat(500) = mat(500) + rxt(161)*y(10) + rxt(198)*y(35) + rxt(216)*y(41)
         mat(381) = mat(381) + rxt(162)*y(1) + rxt(161)*y(2)
         mat(461) = mat(461) + rxt(166)*y(12) + rxt(199)*y(35)
         mat(211) = rxt(166)*y(22) + rxt(220)*y(63)
         mat(430) = rxt(200)*y(35)
         mat(201) = rxt(198)*y(2) + rxt(199)*y(22) + rxt(200)*y(28)
         mat(75) = rxt(216)*y(2)
         mat(160) = rxt(220)*y(12)
         mat(205) = -(rxt(166)*y(22) + rxt(220)*y(63))
         mat(447) = -rxt(166)*y(12)
         mat(155) = -rxt(220)*y(12)
         mat(368) = rxt(165)*y(22)
         mat(447) = mat(447) + rxt(165)*y(10)
         mat(508) = rxt(236)*y(18)
         mat(248) = rxt(236)*y(11)
         mat(531) = (rxt(249)+rxt(254)+rxt(260))*y(35)
         mat(194) = (rxt(249)+rxt(254)+rxt(260))*y(33)
         mat(39) = -(rxt(172)*y(22))
         mat(439) = -rxt(172)*y(13)
         mat(363) = rxt(171)*y(23)
         mat(388) = rxt(171)*y(10)
         mat(362) = rxt(163)*y(11)
         mat(506) = rxt(163)*y(10)
         mat(630) = -(rxt(185)*y(30) + rxt(239)*y(9) + rxt(240)*y(23))
         mat(325) = -rxt(185)*y(16)
         mat(593) = -rxt(239)*y(16)
         mat(411) = -rxt(240)*y(16)
         mat(466) = rxt(241)*y(17)
         mat(38) = rxt(241)*y(22)
         mat(33) = -(rxt(241)*y(22))
         mat(438) = -rxt(241)*y(17)
         mat(612) = rxt(240)*y(23)
         mat(387) = rxt(240)*y(16)
         mat(250) = -(rxt(179)*y(28) + rxt(203)*y(37) + rxt(236)*y(11) + rxt(237) &
                      *y(22) + rxt(238)*y(2))
         mat(420) = -rxt(179)*y(18)
         mat(345) = -rxt(203)*y(18)
         mat(511) = -rxt(236)*y(18)
         mat(450) = -rxt(237)*y(18)
         mat(489) = -rxt(238)*y(18)
         mat(577) = rxt(239)*y(16)
         mat(615) = rxt(239)*y(9) + rxt(185)*y(30)
         mat(309) = rxt(185)*y(16)
         mat(329) = -(rxt(132)*y(4) + rxt(133)*y(1) + (rxt(134) + rxt(135) + rxt(136) &
                      ) * y(23))
         mat(226) = -rxt(132)*y(21)
         mat(599) = -rxt(133)*y(21)
         mat(399) = -(rxt(134) + rxt(135) + rxt(136)) * y(21)
         mat(493) = rxt(144)*y(20) + rxt(137)*y(22)
         mat(559) = rxt(129)*y(20)
         mat(268) = rxt(144)*y(2) + rxt(129)*y(3) + rxt(142)*y(22) + rxt(175)*y(28) &
                      + rxt(218)*y(63)
         mat(65) = rxt(234)*y(22)
         mat(176) = rxt(151)*y(22)
         mat(454) = rxt(137)*y(2) + rxt(142)*y(20) + rxt(234)*y(19) + rxt(151)*y(8) &
                      + rxt(237)*y(18)
         mat(252) = rxt(237)*y(22)
         mat(423) = rxt(175)*y(20)
         mat(158) = rxt(218)*y(20)
         mat(402) = -((rxt(134) + rxt(135) + rxt(136)) * y(21) + rxt(139)*y(22) &
                      + rxt(145)*y(2) + rxt(146)*y(1) + 4._r8*rxt(147)*y(23) + rxt(158) &
                      *y(9) + rxt(170)*y(11) + rxt(171)*y(10) + (rxt(177) + rxt(178) &
                      ) * y(28) + rxt(184)*y(30) + rxt(202)*y(37) + rxt(206)*y(38) &
                      + rxt(240)*y(16))
         mat(331) = -(rxt(134) + rxt(135) + rxt(136)) * y(23)
         mat(457) = -rxt(139)*y(23)
         mat(496) = -rxt(145)*y(23)
         mat(602) = -rxt(146)*y(23)
         mat(584) = -rxt(158)*y(23)
         mat(516) = -rxt(170)*y(23)
         mat(377) = -rxt(171)*y(23)
         mat(426) = -(rxt(177) + rxt(178)) * y(23)
         mat(316) = -rxt(184)*y(23)
         mat(352) = -rxt(202)*y(23)
         mat(290) = -rxt(206)*y(23)
         mat(621) = -rxt(240)*y(23)
         mat(602) = mat(602) + rxt(138)*y(22)
         mat(496) = mat(496) + rxt(238)*y(18) + rxt(148)*y(24)
         mat(228) = rxt(132)*y(21)
         mat(66) = rxt(235)*y(22)
         mat(584) = mat(584) + rxt(239)*y(16)
         mat(457) = mat(457) + rxt(138)*y(1) + rxt(235)*y(19) + rxt(169)*y(11) &
                      + rxt(143)*y(24) + rxt(182)*y(30) + rxt(205)*y(38)
         mat(516) = mat(516) + rxt(169)*y(22) + rxt(236)*y(18)
         mat(621) = mat(621) + rxt(239)*y(9) + rxt(185)*y(30)
         mat(255) = rxt(238)*y(2) + rxt(236)*y(11) + rxt(179)*y(28) + rxt(203)*y(37)
         mat(331) = mat(331) + rxt(132)*y(4)
         mat(78) = rxt(148)*y(2) + rxt(143)*y(22) + rxt(176)*y(28)
         mat(426) = mat(426) + rxt(179)*y(18) + rxt(176)*y(24)
         mat(316) = mat(316) + rxt(182)*y(22) + rxt(185)*y(16)
         mat(352) = mat(352) + rxt(203)*y(18)
         mat(290) = mat(290) + rxt(205)*y(22)
         mat(76) = -(rxt(143)*y(22) + rxt(148)*y(2) + rxt(176)*y(28))
         mat(441) = -rxt(143)*y(24)
         mat(474) = -rxt(148)*y(24)
         mat(414) = -rxt(176)*y(24)
         mat(441) = mat(441) + 2.000_r8*rxt(141)*y(22)
         mat(389) = 2.000_r8*rxt(147)*y(23)
         mat(237) = -(rxt(104)*y(3) + rxt(217)*y(63))
         mat(556) = -rxt(104)*y(74)
         mat(156) = -rxt(217)*y(74)
         mat(266) = rxt(142)*y(22)
         mat(449) = rxt(142)*y(20) + 2.000_r8*rxt(140)*y(22) + rxt(166)*y(12) &
                      + rxt(172)*y(13) + rxt(241)*y(17) + rxt(237)*y(18) + rxt(139) &
                      *y(23) + rxt(143)*y(24) + rxt(193)*y(33) + rxt(197)*y(34) &
                      + rxt(213)*y(39)
         mat(206) = rxt(166)*y(22)
         mat(41) = rxt(172)*y(22)
         mat(34) = rxt(241)*y(22)
         mat(249) = rxt(237)*y(22)
         mat(327) = rxt(136)*y(23)
         mat(394) = rxt(139)*y(22) + rxt(136)*y(21)
         mat(77) = rxt(143)*y(22)
         mat(532) = rxt(193)*y(22) + (rxt(250)+rxt(255)+rxt(261))*y(34) + (rxt(251) &
                       +rxt(262))*y(40)
         mat(185) = rxt(197)*y(22) + (rxt(250)+rxt(255)+rxt(261))*y(33)
         mat(163) = rxt(213)*y(22)
         mat(113) = (rxt(251)+rxt(262))*y(33)
         mat(427) = -(rxt(174)*y(1) + rxt(175)*y(20) + rxt(176)*y(24) + (rxt(177) &
                      + rxt(178)) * y(23) + rxt(179)*y(18) + rxt(196)*y(34) + rxt(200) &
                      *y(35))
         mat(603) = -rxt(174)*y(28)
         mat(270) = -rxt(175)*y(28)
         mat(79) = -rxt(176)*y(28)
         mat(403) = -(rxt(177) + rxt(178)) * y(28)
         mat(256) = -rxt(179)*y(28)
         mat(187) = -rxt(196)*y(28)
         mat(198) = -rxt(200)*y(28)
         mat(497) = rxt(181)*y(30) + rxt(194)*y(33)
         mat(563) = rxt(130)*y(33) + rxt(125)*y(61)
         mat(585) = rxt(186)*y(30)
         mat(458) = rxt(182)*y(30) + rxt(193)*y(33)
         mat(622) = rxt(185)*y(30)
         mat(317) = rxt(181)*y(2) + rxt(186)*y(9) + rxt(182)*y(22) + rxt(185)*y(16) + ( &
                      + 4.000_r8*rxt(188)+2.000_r8*rxt(190))*y(30) + rxt(210)*y(38)
         mat(540) = rxt(194)*y(2) + rxt(130)*y(3) + rxt(193)*y(22)
         mat(291) = rxt(210)*y(30)
         mat(19) = rxt(125)*y(3)
         mat(412) = rxt(200)*y(35)
         mat(302) = 2.000_r8*rxt(189)*y(30)
         mat(526) = (rxt(250)+rxt(255)+rxt(261))*y(34) + (rxt(249)+rxt(254)+rxt(260)) &
                      *y(35)
         mat(183) = (rxt(250)+rxt(255)+rxt(261))*y(33)
         mat(191) = rxt(200)*y(28) + (rxt(249)+rxt(254)+rxt(260))*y(33)
         mat(312) = -(rxt(181)*y(2) + (rxt(182) + rxt(183)) * y(22) + rxt(184)*y(23) &
                      + rxt(185)*y(16) + rxt(186)*y(9) + rxt(187)*y(10) + (4._r8*rxt(188) &
                      + 4._r8*rxt(189) + 4._r8*rxt(190) + 4._r8*rxt(191)) * y(30) &
                      + (rxt(209) + rxt(210) + rxt(211)) * y(38))
         mat(492) = -rxt(181)*y(30)
         mat(453) = -(rxt(182) + rxt(183)) * y(30)
         mat(398) = -rxt(184)*y(30)
         mat(617) = -rxt(185)*y(30)
         mat(580) = -rxt(186)*y(30)
         mat(373) = -rxt(187)*y(30)
         mat(286) = -(rxt(209) + rxt(210) + rxt(211)) * y(30)
         mat(598) = rxt(174)*y(28)
         mat(492) = mat(492) + rxt(195)*y(34) + rxt(198)*y(35)
         mat(453) = mat(453) + rxt(197)*y(34)
         mat(398) = mat(398) + rxt(178)*y(28)
         mat(422) = rxt(174)*y(1) + rxt(178)*y(23) + rxt(196)*y(34)
         mat(186) = rxt(195)*y(2) + rxt(197)*y(22) + rxt(196)*y(28)
         mat(196) = rxt(198)*y(2)
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
         mat(301) = 2.000_r8*rxt(190)*y(30) + rxt(209)*y(38)
         mat(277) = rxt(209)*y(30)
         mat(300) = 2.000_r8*rxt(191)*y(30)
         mat(544) = -(rxt(130)*y(3) + rxt(193)*y(22) + rxt(194)*y(2) + (rxt(249) &
                      + rxt(254) + rxt(260)) * y(35) + (rxt(250) + rxt(255) + rxt(261) &
                      ) * y(34) + (rxt(251) + rxt(262)) * y(40))
         mat(567) = -rxt(130)*y(33)
         mat(462) = -rxt(193)*y(33)
         mat(501) = -rxt(194)*y(33)
         mat(202) = -(rxt(249) + rxt(254) + rxt(260)) * y(33)
         mat(190) = -(rxt(250) + rxt(255) + rxt(261)) * y(33)
         mat(119) = -(rxt(251) + rxt(262)) * y(33)
         mat(274) = rxt(175)*y(28)
         mat(462) = mat(462) + rxt(183)*y(30)
         mat(260) = rxt(179)*y(28)
         mat(407) = rxt(177)*y(28)
         mat(82) = rxt(176)*y(28)
         mat(431) = rxt(175)*y(20) + rxt(179)*y(18) + rxt(177)*y(23) + rxt(176)*y(24) &
                      + rxt(196)*y(34)
         mat(321) = rxt(183)*y(22)
         mat(190) = mat(190) + rxt(196)*y(28)
         mat(184) = -(rxt(195)*y(2) + rxt(196)*y(28) + rxt(197)*y(22) + (rxt(250) &
                      + rxt(255) + rxt(261)) * y(33))
         mat(484) = -rxt(195)*y(34)
         mat(415) = -rxt(196)*y(34)
         mat(445) = -rxt(197)*y(34)
         mat(529) = -(rxt(250) + rxt(255) + rxt(261)) * y(34)
         mat(445) = mat(445) + rxt(199)*y(35)
         mat(392) = rxt(184)*y(30)
         mat(304) = rxt(184)*y(23)
         mat(192) = rxt(199)*y(22)
         mat(193) = -(rxt(198)*y(2) + rxt(199)*y(22) + rxt(200)*y(28) + (rxt(249) &
                      + rxt(254) + rxt(260)) * y(33))
         mat(485) = -rxt(198)*y(35)
         mat(446) = -rxt(199)*y(35)
         mat(416) = -rxt(200)*y(35)
         mat(530) = -(rxt(249) + rxt(254) + rxt(260)) * y(35)
         mat(367) = rxt(187)*y(30)
         mat(305) = rxt(187)*y(10)
         mat(303) = rxt(211)*y(38)
         mat(527) = (rxt(251)+rxt(262))*y(40)
         mat(278) = rxt(211)*y(30)
         mat(111) = (rxt(251)+rxt(262))*y(33)
         mat(350) = -(rxt(201)*y(1) + rxt(202)*y(23) + rxt(203)*y(18))
         mat(600) = -rxt(201)*y(37)
         mat(400) = -rxt(202)*y(37)
         mat(253) = -rxt(203)*y(37)
         mat(494) = rxt(204)*y(38) + rxt(214)*y(39)
         mat(560) = rxt(131)*y(39)
         mat(582) = rxt(207)*y(38)
         mat(455) = rxt(205)*y(38) + rxt(213)*y(39)
         mat(314) = (rxt(209)+rxt(210))*y(38)
         mat(288) = rxt(204)*y(2) + rxt(207)*y(9) + rxt(205)*y(22) + (rxt(209) &
                       +rxt(210))*y(30) + 4.000_r8*rxt(212)*y(38)
         mat(165) = rxt(214)*y(2) + rxt(131)*y(3) + rxt(213)*y(22)
         mat(285) = -(rxt(204)*y(2) + rxt(205)*y(22) + rxt(206)*y(23) + rxt(207)*y(9) &
                      + rxt(208)*y(10) + (rxt(209) + rxt(210) + rxt(211)) * y(30) &
                      + 4._r8*rxt(212)*y(38))
         mat(491) = -rxt(204)*y(38)
         mat(452) = -rxt(205)*y(38)
         mat(397) = -rxt(206)*y(38)
         mat(579) = -rxt(207)*y(38)
         mat(372) = -rxt(208)*y(38)
         mat(311) = -(rxt(209) + rxt(210) + rxt(211)) * y(38)
         mat(597) = rxt(201)*y(37)
         mat(491) = mat(491) + rxt(215)*y(40) + rxt(216)*y(41)
         mat(347) = rxt(201)*y(1)
         mat(114) = rxt(215)*y(2)
         mat(71) = rxt(216)*y(2)
         mat(162) = -(rxt(131)*y(3) + rxt(213)*y(22) + rxt(214)*y(2))
         mat(553) = -rxt(131)*y(39)
         mat(443) = -rxt(213)*y(39)
         mat(482) = -rxt(214)*y(39)
         mat(247) = rxt(203)*y(37)
         mat(391) = rxt(202)*y(37)
         mat(342) = rxt(203)*y(18) + rxt(202)*y(23)
         mat(112) = -(rxt(215)*y(2) + (rxt(251) + rxt(262)) * y(33))
         mat(478) = -rxt(215)*y(40)
         mat(528) = -(rxt(251) + rxt(262)) * y(40)
         mat(390) = rxt(206)*y(38)
         mat(280) = rxt(206)*y(23)
         mat(68) = -(rxt(216)*y(2))
         mat(473) = -rxt(216)*y(41)
         mat(364) = rxt(208)*y(38)
         mat(279) = rxt(208)*y(10)
         mat(91) = -((rxt(265) + rxt(266)) * y(2) + rxt(273)*y(4) + rxt(277)*y(70))
         mat(476) = -(rxt(265) + rxt(266)) * y(65)
         mat(219) = -rxt(273)*y(65)
         mat(140) = -rxt(277)*y(65)
         mat(120) = -(rxt(268)*y(8) + rxt(269)*y(9) + rxt(276)*y(70))
         mat(171) = -rxt(268)*y(66)
         mat(572) = -rxt(269)*y(66)
         mat(142) = -rxt(276)*y(66)
         mat(221) = rxt(273)*y(65) + rxt(270)*y(67) + rxt(263)*y(68) + rxt(284)*y(73)
         mat(93) = rxt(273)*y(4)
         mat(85) = rxt(270)*y(4)
         mat(103) = rxt(263)*y(4)
         mat(56) = rxt(284)*y(4)
         mat(83) = -((rxt(270) + rxt(271)) * y(4) + rxt(272)*y(2))
         mat(218) = -(rxt(270) + rxt(271)) * y(67)
         mat(475) = -rxt(272)*y(67)
         mat(102) = -(rxt(263)*y(4))
         mat(220) = -rxt(263)*y(68)
         mat(477) = rxt(266)*y(65) + rxt(272)*y(67) + rxt(280)*y(72) + rxt(283)*y(73)
         mat(92) = rxt(266)*y(2)
         mat(84) = rxt(272)*y(2)
         mat(141) = rxt(282)*y(72) + rxt(286)*y(73)
         mat(50) = rxt(280)*y(2) + rxt(282)*y(70)
         mat(55) = rxt(283)*y(2) + rxt(286)*y(70)
         mat(129) = -(rxt(275)*y(70))
         mat(143) = -rxt(275)*y(69)
         mat(480) = rxt(265)*y(65)
         mat(222) = rxt(271)*y(67)
         mat(172) = rxt(268)*y(66)
         mat(573) = rxt(269)*y(66)
         mat(94) = rxt(265)*y(2)
         mat(121) = rxt(268)*y(8) + rxt(269)*y(9)
         mat(86) = rxt(271)*y(4)
         mat(59) = -(rxt(149)*y(4) + rxt(150)*y(2))
         mat(217) = -rxt(149)*y(71)
         mat(471) = -rxt(150)*y(71)
         mat(471) = mat(471) + rxt(265)*y(65)
         mat(90) = rxt(265)*y(2) + .900_r8*rxt(277)*y(70)
         mat(128) = .800_r8*rxt(275)*y(70)
         mat(138) = .900_r8*rxt(277)*y(65) + .800_r8*rxt(275)*y(69)
         mat(144) = -(rxt(275)*y(69) + rxt(276)*y(66) + rxt(277)*y(65))
         mat(130) = -rxt(275)*y(70)
         mat(122) = -rxt(276)*y(70)
         mat(95) = -rxt(277)*y(70)
         mat(46) = -(rxt(280)*y(2) + (rxt(281) + rxt(282)) * y(70))
         mat(469) = -rxt(280)*y(72)
         mat(136) = -(rxt(281) + rxt(282)) * y(72)
         mat(53) = -(rxt(283)*y(2) + rxt(284)*y(4) + rxt(286)*y(70))
         mat(470) = -rxt(283)*y(73)
         mat(216) = -rxt(284)*y(73)
         mat(137) = -rxt(286)*y(73)
         mat(137) = mat(137) + rxt(281)*y(72)
         mat(47) = rxt(281)*y(70)
         mat(12) = -(rxt(124)*y(3))
         mat(550) = -rxt(124)*y(60)
         mat(17) = -(rxt(125)*y(3))
         mat(551) = -rxt(125)*y(61)
         mat(263) = rxt(218)*y(63)
         mat(203) = rxt(220)*y(63)
         mat(234) = rxt(217)*y(63)
         mat(153) = rxt(218)*y(20) + rxt(220)*y(12) + rxt(217)*y(74)
         mat(154) = -(rxt(217)*y(74) + rxt(218)*y(20) + rxt(220)*y(12))
         mat(235) = -rxt(217)*y(63)
         mat(264) = -rxt(218)*y(63)
         mat(204) = -rxt(220)*y(63)
         mat(552) = 2.000_r8*rxt(124)*y(60) + rxt(125)*y(61)
         mat(13) = 2.000_r8*rxt(124)*y(3)
         mat(18) = rxt(125)*y(3)
      end subroutine nlnmat03
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
         mat( 4) = mat( 4) + lmat( 4)
         mat( 5) = mat( 5) + lmat( 5)
         mat( 6) = mat( 6) + lmat( 6)
         mat( 7) = mat( 7) + lmat( 7)
         mat( 8) = lmat( 8)
         mat( 9) = lmat( 9)
         mat( 10) = lmat( 10)
         mat( 11) = lmat( 11)
         mat( 12) = mat( 12) + lmat( 12)
         mat( 13) = mat( 13) + lmat( 13)
         mat( 15) = lmat( 15)
         mat( 16) = lmat( 16)
         mat( 17) = mat( 17) + lmat( 17)
         mat( 18) = mat( 18) + lmat( 18)
         mat( 19) = mat( 19) + lmat( 19)
         mat( 21) = lmat( 21)
         mat( 22) = lmat( 22)
         mat( 23) = lmat( 23)
         mat( 24) = lmat( 24)
         mat( 25) = lmat( 25)
         mat( 26) = lmat( 26)
         mat( 27) = lmat( 27)
         mat( 28) = lmat( 28)
         mat( 29) = lmat( 29)
         mat( 30) = lmat( 30)
         mat( 31) = lmat( 31)
         mat( 32) = lmat( 32)
         mat( 33) = mat( 33) + lmat( 33)
         mat( 35) = lmat( 35)
         mat( 36) = lmat( 36)
         mat( 37) = mat( 37) + lmat( 37)
         mat( 39) = mat( 39) + lmat( 39)
         mat( 42) = mat( 42) + lmat( 42)
         mat( 43) = lmat( 43)
         mat( 44) = mat( 44) + lmat( 44)
         mat( 45) = lmat( 45)
         mat( 46) = mat( 46) + lmat( 46)
         mat( 47) = mat( 47) + lmat( 47)
         mat( 48) = lmat( 48)
         mat( 49) = lmat( 49)
         mat( 50) = mat( 50) + lmat( 50)
         mat( 51) = lmat( 51)
         mat( 52) = lmat( 52)
         mat( 53) = mat( 53) + lmat( 53)
         mat( 54) = lmat( 54)
         mat( 55) = mat( 55) + lmat( 55)
         mat( 58) = mat( 58) + lmat( 58)
         mat( 59) = mat( 59) + lmat( 59)
         mat( 64) = mat( 64) + lmat( 64)
         mat( 68) = mat( 68) + lmat( 68)
         mat( 69) = lmat( 69)
         mat( 70) = lmat( 70)
         mat( 71) = mat( 71) + lmat( 71)
         mat( 72) = lmat( 72)
         mat( 73) = lmat( 73)
         mat( 75) = mat( 75) + lmat( 75)
         mat( 76) = mat( 76) + lmat( 76)
         mat( 80) = mat( 80) + lmat( 80)
         mat( 83) = mat( 83) + lmat( 83)
         mat( 91) = mat( 91) + lmat( 91)
         mat( 101) = lmat( 101)
         mat( 102) = mat( 102) + lmat( 102)
         mat( 103) = mat( 103) + lmat( 103)
         mat( 104) = lmat( 104)
         mat( 105) = lmat( 105)
         mat( 112) = mat( 112) + lmat( 112)
         mat( 115) = lmat( 115)
         mat( 117) = mat( 117) + lmat( 117)
         mat( 120) = mat( 120) + lmat( 120)
         mat( 121) = mat( 121) + lmat( 121)
         mat( 127) = mat( 127) + lmat( 127)
         mat( 129) = mat( 129) + lmat( 129)
         mat( 144) = mat( 144) + lmat( 144)
         mat( 153) = mat( 153) + lmat( 153)
         mat( 154) = mat( 154) + lmat( 154)
         mat( 161) = lmat( 161)
         mat( 162) = mat( 162) + lmat( 162)
         mat( 164) = lmat( 164)
         mat( 165) = mat( 165) + lmat( 165)
         mat( 169) = lmat( 169)
         mat( 173) = lmat( 173)
         mat( 174) = mat( 174) + lmat( 174)
         mat( 184) = mat( 184) + lmat( 184)
         mat( 187) = mat( 187) + lmat( 187)
         mat( 188) = mat( 188) + lmat( 188)
         mat( 192) = mat( 192) + lmat( 192)
         mat( 193) = mat( 193) + lmat( 193)
         mat( 194) = mat( 194) + lmat( 194)
         mat( 196) = mat( 196) + lmat( 196)
         mat( 197) = lmat( 197)
         mat( 198) = mat( 198) + lmat( 198)
         mat( 201) = mat( 201) + lmat( 201)
         mat( 205) = mat( 205) + lmat( 205)
         mat( 209) = lmat( 209)
         mat( 210) = mat( 210) + lmat( 210)
         mat( 215) = lmat( 215)
         mat( 216) = mat( 216) + lmat( 216)
         mat( 220) = mat( 220) + lmat( 220)
         mat( 221) = mat( 221) + lmat( 221)
         mat( 223) = lmat( 223)
         mat( 225) = mat( 225) + lmat( 225)
         mat( 230) = mat( 230) + lmat( 230)
         mat( 231) = mat( 231) + lmat( 231)
         mat( 237) = mat( 237) + lmat( 237)
         mat( 238) = lmat( 238)
         mat( 239) = lmat( 239)
         mat( 241) = mat( 241) + lmat( 241)
         mat( 242) = lmat( 242)
         mat( 244) = mat( 244) + lmat( 244)
         mat( 246) = mat( 246) + lmat( 246)
         mat( 250) = mat( 250) + lmat( 250)
         mat( 251) = lmat( 251)
         mat( 252) = mat( 252) + lmat( 252)
         mat( 267) = mat( 267) + lmat( 267)
         mat( 285) = mat( 285) + lmat( 285)
         mat( 288) = mat( 288) + lmat( 288)
         mat( 293) = mat( 293) + lmat( 293)
         mat( 312) = mat( 312) + lmat( 312)
         mat( 317) = mat( 317) + lmat( 317)
         mat( 319) = mat( 319) + lmat( 319)
         mat( 329) = mat( 329) + lmat( 329)
         mat( 350) = mat( 350) + lmat( 350)
         mat( 368) = mat( 368) + lmat( 368)
         mat( 376) = mat( 376) + lmat( 376)
         mat( 379) = mat( 379) + lmat( 379)
         mat( 380) = mat( 380) + lmat( 380)
         mat( 384) = mat( 384) + lmat( 384)
         mat( 389) = mat( 389) + lmat( 389)
         mat( 402) = mat( 402) + lmat( 402)
         mat( 413) = mat( 413) + lmat( 413)
         mat( 424) = lmat( 424)
         mat( 426) = mat( 426) + lmat( 426)
         mat( 427) = mat( 427) + lmat( 427)
         mat( 431) = mat( 431) + lmat( 431)
         mat( 435) = lmat( 435)
         mat( 436) = lmat( 436)
         mat( 437) = lmat( 437)
         mat( 449) = mat( 449) + lmat( 449)
         mat( 455) = mat( 455) + lmat( 455)
         mat( 457) = mat( 457) + lmat( 457)
         mat( 458) = mat( 458) + lmat( 458)
         mat( 459) = mat( 459) + lmat( 459)
         mat( 466) = mat( 466) + lmat( 466)
         mat( 469) = mat( 469) + lmat( 469)
         mat( 470) = mat( 470) + lmat( 470)
         mat( 477) = mat( 477) + lmat( 477)
         mat( 481) = lmat( 481)
         mat( 499) = mat( 499) + lmat( 499)
         mat( 508) = mat( 508) + lmat( 508)
         mat( 509) = mat( 509) + lmat( 509)
         mat( 515) = mat( 515) + lmat( 515)
         mat( 519) = mat( 519) + lmat( 519)
         mat( 520) = mat( 520) + lmat( 520)
         mat( 523) = mat( 523) + lmat( 523)
         mat( 536) = lmat( 536)
         mat( 540) = mat( 540) + lmat( 540)
         mat( 544) = mat( 544) + lmat( 544)
         mat( 550) = mat( 550) + lmat( 550)
         mat( 551) = mat( 551) + lmat( 551)
         mat( 552) = mat( 552) + lmat( 552)
         mat( 555) = mat( 555) + lmat( 555)
         mat( 557) = lmat( 557)
         mat( 558) = mat( 558) + lmat( 558)
         mat( 559) = mat( 559) + lmat( 559)
         mat( 560) = mat( 560) + lmat( 560)
         mat( 562) = lmat( 562)
         mat( 563) = mat( 563) + lmat( 563)
         mat( 564) = mat( 564) + lmat( 564)
         mat( 565) = mat( 565) + lmat( 565)
         mat( 568) = mat( 568) + lmat( 568)
         mat( 569) = lmat( 569)
         mat( 571) = lmat( 571)
         mat( 573) = mat( 573) + lmat( 573)
         mat( 574) = lmat( 574)
         mat( 575) = mat( 575) + lmat( 575)
         mat( 587) = mat( 587) + lmat( 587)
         mat( 591) = mat( 591) + lmat( 591)
         mat( 594) = mat( 594) + lmat( 594)
         mat( 596) = mat( 596) + lmat( 596)
         mat( 605) = mat( 605) + lmat( 605)
         mat( 608) = mat( 608) + lmat( 608)
         mat( 610) = mat( 610) + lmat( 610)
         mat( 630) = mat( 630) + lmat( 630)
         mat( 99) = 0._r8
         mat( 100) = 0._r8
         mat( 107) = 0._r8
         mat( 108) = 0._r8
         mat( 109) = 0._r8
         mat( 116) = 0._r8
         mat( 132) = 0._r8
         mat( 134) = 0._r8
         mat( 135) = 0._r8
         mat( 139) = 0._r8
         mat( 146) = 0._r8
         mat( 147) = 0._r8
         mat( 148) = 0._r8
         mat( 149) = 0._r8
         mat( 152) = 0._r8
         mat( 170) = 0._r8
         mat( 178) = 0._r8
         mat( 181) = 0._r8
         mat( 195) = 0._r8
         mat( 207) = 0._r8
         mat( 208) = 0._r8
         mat( 212) = 0._r8
         mat( 227) = 0._r8
         mat( 229) = 0._r8
         mat( 236) = 0._r8
         mat( 240) = 0._r8
         mat( 243) = 0._r8
         mat( 245) = 0._r8
         mat( 254) = 0._r8
         mat( 261) = 0._r8
         mat( 262) = 0._r8
         mat( 265) = 0._r8
         mat( 269) = 0._r8
         mat( 273) = 0._r8
         mat( 276) = 0._r8
         mat( 281) = 0._r8
         mat( 283) = 0._r8
         mat( 284) = 0._r8
         mat( 287) = 0._r8
         mat( 294) = 0._r8
         mat( 295) = 0._r8
         mat( 296) = 0._r8
         mat( 298) = 0._r8
         mat( 299) = 0._r8
         mat( 306) = 0._r8
         mat( 308) = 0._r8
         mat( 310) = 0._r8
         mat( 313) = 0._r8
         mat( 320) = 0._r8
         mat( 322) = 0._r8
         mat( 324) = 0._r8
         mat( 330) = 0._r8
         mat( 332) = 0._r8
         mat( 335) = 0._r8
         mat( 336) = 0._r8
         mat( 337) = 0._r8
         mat( 338) = 0._r8
         mat( 340) = 0._r8
         mat( 344) = 0._r8
         mat( 346) = 0._r8
         mat( 348) = 0._r8
         mat( 349) = 0._r8
         mat( 351) = 0._r8
         mat( 353) = 0._r8
         mat( 354) = 0._r8
         mat( 355) = 0._r8
         mat( 356) = 0._r8
         mat( 357) = 0._r8
         mat( 358) = 0._r8
         mat( 359) = 0._r8
         mat( 361) = 0._r8
         mat( 365) = 0._r8
         mat( 370) = 0._r8
         mat( 371) = 0._r8
         mat( 374) = 0._r8
         mat( 375) = 0._r8
         mat( 378) = 0._r8
         mat( 382) = 0._r8
         mat( 383) = 0._r8
         mat( 386) = 0._r8
         mat( 395) = 0._r8
         mat( 408) = 0._r8
         mat( 417) = 0._r8
         mat( 419) = 0._r8
         mat( 425) = 0._r8
         mat( 429) = 0._r8
         mat( 432) = 0._r8
         mat( 433) = 0._r8
         mat( 442) = 0._r8
         mat( 463) = 0._r8
         mat( 479) = 0._r8
         mat( 486) = 0._r8
         mat( 488) = 0._r8
         mat( 502) = 0._r8
         mat( 505) = 0._r8
         mat( 510) = 0._r8
         mat( 512) = 0._r8
         mat( 513) = 0._r8
         mat( 514) = 0._r8
         mat( 517) = 0._r8
         mat( 521) = 0._r8
         mat( 522) = 0._r8
         mat( 524) = 0._r8
         mat( 525) = 0._r8
         mat( 533) = 0._r8
         mat( 534) = 0._r8
         mat( 535) = 0._r8
         mat( 537) = 0._r8
         mat( 538) = 0._r8
         mat( 539) = 0._r8
         mat( 543) = 0._r8
         mat( 546) = 0._r8
         mat( 547) = 0._r8
         mat( 548) = 0._r8
         mat( 554) = 0._r8
         mat( 561) = 0._r8
         mat( 566) = 0._r8
         mat( 578) = 0._r8
         mat( 581) = 0._r8
         mat( 589) = 0._r8
         mat( 590) = 0._r8
         mat( 607) = 0._r8
         mat( 611) = 0._r8
         mat( 614) = 0._r8
         mat( 616) = 0._r8
         mat( 618) = 0._r8
         mat( 619) = 0._r8
         mat( 623) = 0._r8
         mat( 624) = 0._r8
         mat( 625) = 0._r8
         mat( 626) = 0._r8
         mat( 627) = 0._r8
         mat( 629) = 0._r8
         mat( 1) = mat( 1) - dti
         mat( 4) = mat( 4) - dti
         mat( 7) = mat( 7) - dti
         mat( 9) = mat( 9) - dti
         mat( 12) = mat( 12) - dti
         mat( 15) = mat( 15) - dti
         mat( 17) = mat( 17) - dti
         mat( 21) = mat( 21) - dti
         mat( 24) = mat( 24) - dti
         mat( 27) = mat( 27) - dti
         mat( 33) = mat( 33) - dti
         mat( 39) = mat( 39) - dti
         mat( 46) = mat( 46) - dti
         mat( 53) = mat( 53) - dti
         mat( 59) = mat( 59) - dti
         mat( 64) = mat( 64) - dti
         mat( 68) = mat( 68) - dti
         mat( 76) = mat( 76) - dti
         mat( 83) = mat( 83) - dti
         mat( 91) = mat( 91) - dti
         mat( 102) = mat( 102) - dti
         mat( 112) = mat( 112) - dti
         mat( 120) = mat( 120) - dti
         mat( 129) = mat( 129) - dti
         mat( 144) = mat( 144) - dti
         mat( 154) = mat( 154) - dti
         mat( 162) = mat( 162) - dti
         mat( 174) = mat( 174) - dti
         mat( 184) = mat( 184) - dti
         mat( 193) = mat( 193) - dti
         mat( 205) = mat( 205) - dti
         mat( 225) = mat( 225) - dti
         mat( 237) = mat( 237) - dti
         mat( 250) = mat( 250) - dti
         mat( 267) = mat( 267) - dti
         mat( 285) = mat( 285) - dti
         mat( 312) = mat( 312) - dti
         mat( 329) = mat( 329) - dti
         mat( 350) = mat( 350) - dti
         mat( 376) = mat( 376) - dti
         mat( 402) = mat( 402) - dti
         mat( 427) = mat( 427) - dti
         mat( 459) = mat( 459) - dti
         mat( 499) = mat( 499) - dti
         mat( 520) = mat( 520) - dti
         mat( 544) = mat( 544) - dti
         mat( 568) = mat( 568) - dti
         mat( 591) = mat( 591) - dti
         mat( 610) = mat( 610) - dti
         mat( 630) = mat( 630) - dti
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
      call nlnmat_finit( mat, lmat, dti )
      end subroutine nlnmat
      end module mo_nln_matrix
