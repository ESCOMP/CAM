
      module mo_sim_dat

      private
      public :: set_sim_dat

      contains

      subroutine set_sim_dat

      use chem_mods,     only : clscnt, cls_rxt_cnt, clsmap, permute, adv_mass, fix_mass, crb_mass
      use chem_mods,     only : diag_map
      use chem_mods,     only : phtcnt, rxt_tag_cnt, rxt_tag_lst, rxt_tag_map
      use chem_mods,     only : pht_alias_lst, pht_alias_mult
      use chem_mods,     only : extfrc_lst, inv_lst, slvd_lst
      use chem_mods,     only : enthalpy_cnt, cph_enthalpy, cph_rid, num_rnts, rxntot
      use cam_abortutils,only : endrun
      use mo_tracname,   only : solsym
      use chem_mods,     only : frc_from_dataset
      use chem_mods,     only : is_scalar, is_vector
      use shr_kind_mod,  only : r8 => shr_kind_r8
      use cam_logfile,   only : iulog

      implicit none

!--------------------------------------------------------------
!      ... local variables
!--------------------------------------------------------------
      integer :: ios

      is_scalar = .true.
      is_vector = .false.

      clscnt(:) = (/     22,     0,     0,   113,     0 /)

      cls_rxt_cnt(:,1) = (/      1,    57,     0,    22 /)
      cls_rxt_cnt(:,4) = (/     29,   151,   413,   113 /)

      solsym(:135) = (/ 'O3              ','O               ','O2              ','N2O             ','N               ', &
                        'NO              ','NO2             ','NO3             ','HNO3            ','HO2NO2          ', &
                        'N2O5            ','CH4             ','CH3O2           ','CH3OOH          ','CH2O            ', &
                        'CO              ','H2              ','H               ','H2O2            ','HONO            ', &
                        'CLY             ','BRY             ','CL2             ','CLO             ','OCLO            ', &
                        'CL2O2           ','HCL             ','HOCL            ','CLONO2          ','BRCL            ', &
                        'BRO             ','HBR             ','HOBR            ','BRONO2          ','CH3CL           ', &
                        'CH3BR           ','CFC11           ','CFC12           ','CFC113          ','HCFC22          ', &
                        'CCL4            ','CH3CCL3         ','CF3BR           ','CF2CLBR         ','HCFC141B        ', &
                        'HCFC142B        ','CFC114          ','CFC115          ','H1202           ','H2402           ', &
                        'CHBR3           ','CH2BR2          ','COF2            ','COFCL           ','HF              ', &
                        'F               ','CO2             ','OCS             ','S               ','SO              ', &
                        'SO2             ','SO3             ','H2SO4           ','DMS             ','SOAG            ', &
                        'so4_a1          ','pom_a1          ','soa_a1          ','bc_a1           ','dst_a1          ', &
                        'ncl_a1          ','num_a1          ','so4_a2          ','soa_a2          ','ncl_a2          ', &
                        'num_a2          ','dst_a2          ','dst_a3          ','ncl_a3          ','so4_a3          ', &
                        'num_a3          ','pom_a4          ','bc_a4           ','num_a4          ','CL              ', &
                        'BR              ','OH              ','HO2             ','N2p             ','O2p             ', &
                        'O4p             ','O2p_H2O         ','Hp_H2O          ','Hp_2H2O         ','Hp_3H2O         ', &
                        'Hp_4H2O         ','Hp_5H2O         ','H3Op_OH         ','Hp_3N1          ','Hp_4N1          ', &
                        'NOp_H2O         ','NOp_2H2O        ','NOp_3H2O        ','NOp_CO2         ','NOp_N2          ', &
                        'Om              ','O2m             ','O3m             ','O4m             ','OHm             ', &
                        'CO3m            ','CO4m            ','NO2m            ','NO3m            ','CO3m_H2O        ', &
                        'CO3m2H2O        ','NO2m_H2O        ','NO3m_H2O        ','NO3m2H2O        ','NO3mHNO3        ', &
                        'NO3m_HCL        ','HCO3m           ','CLm             ','CLOm            ','CLm_H2O         ', &
                        'CLm_HCL         ','Np              ','Op              ','NOp             ','e               ', &
                        'N2D             ','O2_1S           ','O2_1D           ','O1D             ','H2O             ' /)

      adv_mass(:135) = (/    47.998200_r8,    15.999400_r8,    31.998800_r8,    44.012880_r8,    14.006740_r8, &
                             30.006140_r8,    46.005540_r8,    62.004940_r8,    63.012340_r8,    79.011740_r8, &
                            108.010480_r8,    16.040600_r8,    47.032000_r8,    48.039400_r8,    30.025200_r8, &
                             28.010400_r8,     2.014800_r8,     1.007400_r8,    34.013600_r8,    47.012940_r8, &
                            100.916850_r8,    99.716850_r8,    70.905400_r8,    51.452100_r8,    67.451500_r8, &
                            102.904200_r8,    36.460100_r8,    52.459500_r8,    97.457640_r8,   115.356700_r8, &
                             95.903400_r8,    80.911400_r8,    96.910800_r8,   141.908940_r8,    50.485900_r8, &
                             94.937200_r8,   137.367503_r8,   120.913206_r8,   187.375310_r8,    86.467906_r8, &
                            153.821800_r8,   133.402300_r8,   148.910210_r8,   165.364506_r8,   116.948003_r8, &
                            100.493706_r8,   170.921013_r8,   154.466716_r8,   209.815806_r8,   259.823613_r8, &
                            252.730400_r8,   173.833800_r8,    66.007206_r8,    82.461503_r8,    20.005803_r8, &
                             18.998403_r8,    44.009800_r8,    60.076400_r8,    32.066000_r8,    48.065400_r8, &
                             64.064800_r8,    80.064200_r8,    98.078400_r8,    62.132400_r8,    12.011000_r8, &
                            115.107340_r8,    12.011000_r8,    12.011000_r8,    12.011000_r8,   135.064039_r8, &
                             58.442468_r8,     1.007400_r8,   115.107340_r8,    12.011000_r8,    58.442468_r8, &
                              1.007400_r8,   135.064039_r8,   135.064039_r8,    58.442468_r8,   115.107340_r8, &
                              1.007400_r8,    12.011000_r8,    12.011000_r8,     1.007400_r8,    35.452700_r8, &
                             79.904000_r8,    17.006800_r8,    33.006200_r8,    28.013480_r8,    31.998800_r8, &
                             63.997600_r8,    50.013000_r8,    19.021600_r8,    37.035800_r8,    55.050000_r8, &
                             73.064200_r8,    91.078400_r8,    36.028400_r8,   118.062340_r8,   136.076540_r8, &
                             48.020340_r8,    66.034540_r8,    68.049340_r8,    74.015940_r8,    58.019620_r8, &
                             15.999400_r8,    31.998800_r8,    47.998200_r8,    63.997600_r8,    17.006800_r8, &
                             60.009200_r8,    76.008600_r8,    46.005540_r8,    62.004940_r8,    78.023400_r8, &
                             96.037600_r8,    64.019740_r8,    80.019140_r8,    98.033340_r8,   125.017280_r8, &
                             98.465040_r8,    61.016600_r8,    35.452700_r8,    51.452100_r8,    53.466900_r8, &
                             71.912800_r8,    14.006740_r8,    15.999400_r8,    30.006140_r8, 0.548567E-03_r8, &
                             14.006740_r8,    31.998800_r8,    31.998800_r8,    15.999400_r8,    18.014200_r8 /)

      crb_mass(:135) = (/     0.000000_r8,     0.000000_r8,     0.000000_r8,     0.000000_r8,     0.000000_r8, &
                              0.000000_r8,     0.000000_r8,     0.000000_r8,     0.000000_r8,     0.000000_r8, &
                              0.000000_r8,    12.011000_r8,    12.011000_r8,    12.011000_r8,    12.011000_r8, &
                             12.011000_r8,     0.000000_r8,     0.000000_r8,     0.000000_r8,     0.000000_r8, &
                             12.011000_r8,     0.000000_r8,     0.000000_r8,     0.000000_r8,     0.000000_r8, &
                              0.000000_r8,     0.000000_r8,     0.000000_r8,     0.000000_r8,     0.000000_r8, &
                              0.000000_r8,     0.000000_r8,     0.000000_r8,     0.000000_r8,    12.011000_r8, &
                             12.011000_r8,    12.011000_r8,    12.011000_r8,    24.022000_r8,    12.011000_r8, &
                             12.011000_r8,    24.022000_r8,    12.011000_r8,    12.011000_r8,    24.022000_r8, &
                             24.022000_r8,    24.022000_r8,    24.022000_r8,    12.011000_r8,    24.022000_r8, &
                             12.011000_r8,    12.011000_r8,    12.011000_r8,    12.011000_r8,     0.000000_r8, &
                              0.000000_r8,    12.011000_r8,    12.011000_r8,     0.000000_r8,     0.000000_r8, &
                              0.000000_r8,     0.000000_r8,     0.000000_r8,    24.022000_r8,    12.011000_r8, &
                              0.000000_r8,    12.011000_r8,    12.011000_r8,    12.011000_r8,     0.000000_r8, &
                              0.000000_r8,     0.000000_r8,     0.000000_r8,    12.011000_r8,     0.000000_r8, &
                              0.000000_r8,     0.000000_r8,     0.000000_r8,     0.000000_r8,     0.000000_r8, &
                              0.000000_r8,    12.011000_r8,    12.011000_r8,     0.000000_r8,     0.000000_r8, &
                              0.000000_r8,     0.000000_r8,     0.000000_r8,     0.000000_r8,     0.000000_r8, &
                              0.000000_r8,     0.000000_r8,     0.000000_r8,     0.000000_r8,     0.000000_r8, &
                              0.000000_r8,     0.000000_r8,     0.000000_r8,     0.000000_r8,     0.000000_r8, &
                              0.000000_r8,     0.000000_r8,     0.000000_r8,    12.011000_r8,     0.000000_r8, &
                              0.000000_r8,     0.000000_r8,     0.000000_r8,     0.000000_r8,     0.000000_r8, &
                             12.011000_r8,    12.011000_r8,     0.000000_r8,     0.000000_r8,    12.011000_r8, &
                             12.011000_r8,     0.000000_r8,     0.000000_r8,     0.000000_r8,     0.000000_r8, &
                              0.000000_r8,    12.011000_r8,     0.000000_r8,     0.000000_r8,     0.000000_r8, &
                              0.000000_r8,     0.000000_r8,     0.000000_r8,     0.000000_r8,     0.000000_r8, &
                              0.000000_r8,     0.000000_r8,     0.000000_r8,     0.000000_r8,     0.000000_r8 /)

      fix_mass(:  2) = (/ 0.00000000_r8, 28.0134800_r8 /)

      clsmap(: 22,1) = (/   12,   4,  35,  36,  37,  38,  39,  47,  48,  40, &
                            45,  46,  41,  42,  43,  44,  49,  50,  51,  52, &
                            21,  22 /)
      clsmap(:113,4) = (/    1,   2, 134,   3, 132, 133,  17,  16,  57,   5, &
                             6,   7,  87,   8,  20,   9,  10,  11,  13,  14, &
                            15,  18,  88,  19, 135,  85,  23,  24,  25,  26, &
                            27,  28,  29,  30,  86,  31,  32,  33,  34,  89, &
                            90, 127, 128, 129, 131,  91,  92,  93,  94,  95, &
                            96,  97,  98,  99, 100, 101, 102, 103, 104, 105, &
                           106, 107, 108, 109, 111, 112, 113, 114, 110, 122, &
                           115, 116, 117, 118, 119, 120, 121, 123, 124, 125, &
                           126, 130,  53,  54,  55,  56,  58,  59,  60,  61, &
                            62,  63,  64,  65,  66,  67,  68,  69,  70,  71, &
                            72,  73,  74,  75,  76,  77,  78,  79,  80,  81, &
                            82,  83,  84 /)

      permute(:113,4) = (/   94,  95,  82,  97,  32,  61,  80,  50, 102,  64, &
                            104, 105,  88, 107,  51, 109,  39,  62,  63,  37, &
                             73,  91,  86,  43, 100, 112,  24, 106,  33,  21, &
                             93,  56,  59,  26,  81,  75,  55,  52,  41,  45, &
                            110,  38,  57, 113,  40,  65,  49,  54,  44, 101, &
                            108, 103,  35,  28,  29,  98,  90,  36,  47,  46, &
                             99, 111,  74,  42,  89,  84,  87,  85,  92,  66, &
                             78,  68,  71,  79,  77,  76,  72,  83,  60,  70, &
                             69,  96,  22,  25,  27,  53,  34,  48,  67,  58, &
                             30,  23,  31,   1,   2,   3,   4,   5,   6,   7, &
                              8,   9,  10,  11,  12,  13,  14,  15,  16,  17, &
                             18,  19,  20 /)

      diag_map(:113) = (/    1,   2,   3,   4,   5,   6,   7,   8,   9,  10, &
                            11,  12,  13,  14,  15,  16,  17,  18,  19,  20, &
                            21,  24,  27,  30,  32,  36,  39,  42,  46,  51, &
                            55,  61,  64,  69,  77,  83,  89,  95, 102, 109, &
                           114, 122, 129, 136, 143, 153, 161, 168, 177, 186, &
                           191, 197, 206, 215, 220, 228, 236, 247, 255, 265, &
                           273, 283, 296, 309, 324, 336, 353, 367, 382, 398, &
                           413, 427, 445, 459, 480, 496, 512, 532, 554, 575, &
                           600, 631, 655, 677, 710, 744, 770, 819, 849, 893, &
                           931, 965,1010,1051,1110,1151,1198,1241,1278,1338, &
                          1382,1424,1468,1512,1559,1606,1648,1692,1733,1780, &
                          1817,1860,1907 /)

      extfrc_lst(: 23) = (/ 'NO              ','NO2             ','CO              ','SO2             ','DMS             ', &
                            'so4_a1          ','so4_a2          ','pom_a1          ','pom_a4          ','bc_a1           ', &
                            'bc_a4           ','num_a1          ','num_a2          ','num_a4          ','Op              ', &
                            'O2p             ','Np              ','N2p             ','N2D             ','O               ', &
                            'N               ','e               ','OH              ' /)

      frc_from_dataset(: 23) = (/ .true., .true., .true., .true., .true., &
                                  .true., .true., .true., .true., .true., &
                                  .true., .true., .true., .true., .false., &
                                  .false., .false., .false., .false., .false., &
                                  .false., .false., .false. /)

      inv_lst(:  2) = (/ 'M               ', 'N2              ' /)

      slvd_lst(: 50) = (/ 'CL              ', 'BR              ', 'OH              ', 'HO2             ', 'Op              ', &
                          'O2p             ', 'NOp             ', 'Np              ', 'N2p             ', 'e               ', &
                          'O2_1S           ', 'O2_1D           ', 'N2D             ', 'O1D             ', 'O4p             ', &
                          'O2p_H2O         ', 'Hp_H2O          ', 'Hp_2H2O         ', 'Hp_3H2O         ', 'Hp_4H2O         ', &
                          'Hp_5H2O         ', 'H3Op_OH         ', 'Hp_3N1          ', 'Hp_4N1          ', 'NOp_H2O         ', &
                          'NOp_2H2O        ', 'NOp_3H2O        ', 'NOp_CO2         ', 'NOp_N2          ', 'Om              ', &
                          'O2m             ', 'O3m             ', 'O4m             ', 'CO3m            ', 'CO4m            ', &
                          'NO2m            ', 'NO3m            ', 'OHm             ', 'HCO3m           ', 'CO3m_H2O        ', &
                          'CO3m2H2O        ', 'NO2m_H2O        ', 'NO3m_H2O        ', 'NO3m2H2O        ', 'NO3mHNO3        ', &
                          'NO3m_HCL        ', 'CLm             ', 'CLOm            ', 'CLm_H2O         ', 'CLm_HCL         ' /)

      if( allocated( rxt_tag_lst ) ) then
         deallocate( rxt_tag_lst )
      end if
      allocate( rxt_tag_lst(rxt_tag_cnt),stat=ios )
      if( ios /= 0 ) then
         write(iulog,*) 'set_sim_dat: failed to allocate rxt_tag_lst; error = ',ios
         call endrun
      end if
      if( allocated( rxt_tag_map ) ) then
         deallocate( rxt_tag_map )
      end if
      allocate( rxt_tag_map(rxt_tag_cnt),stat=ios )
      if( ios /= 0 ) then
         write(iulog,*) 'set_sim_dat: failed to allocate rxt_tag_map; error = ',ios
         call endrun
      end if
      rxt_tag_lst(     1:   200) = (/ 'jo2_a                           ', 'jo2_b                           ', &
                                      'jo3_a                           ', 'jo3_b                           ', &
                                      'jn2o                            ', 'jno                             ', &
                                      'jno_i                           ', 'jno2                            ', &
                                      'jn2o5_a                         ', 'jn2o5_b                         ', &
                                      'jhno3                           ', 'jno3_a                          ', &
                                      'jno3_b                          ', 'jho2no2_a                       ', &
                                      'jho2no2_b                       ', 'jch3ooh                         ', &
                                      'jch2o_a                         ', 'jch2o_b                         ', &
                                      'jh2o_a                          ', 'jh2o_b                          ', &
                                      'jh2o_c                          ', 'jh2o2                           ', &
                                      'jcl2                            ', 'jclo                            ', &
                                      'joclo                           ', 'jcl2o2                          ', &
                                      'jhocl                           ', 'jhcl                            ', &
                                      'jclono2_a                       ', 'jclono2_b                       ', &
                                      'jbrcl                           ', 'jbro                            ', &
                                      'jhobr                           ', 'jhbr                            ', &
                                      'jbrono2_a                       ', 'jbrono2_b                       ', &
                                      'jch3cl                          ', 'jccl4                           ', &
                                      'jch3ccl3                        ', 'jcfcl3                          ', &
                                      'jcf2cl2                         ', 'jcfc113                         ', &
                                      'jcfc114                         ', 'jcfc115                         ', &
                                      'jhcfc22                         ', 'jhcfc141b                       ', &
                                      'jhcfc142b                       ', 'jch3br                          ', &
                                      'jcf3br                          ', 'jcf2clbr                        ', &
                                      'jchbr3                          ', 'jch2br2                         ', &
                                      'jh1202                          ', 'jh2402                          ', &
                                      'jcof2                           ', 'jcofcl                          ', &
                                      'jhf                             ', 'jco2                            ', &
                                      'jch4_a                          ', 'jch4_b                          ', &
                                      'jh2so4                          ', 'jso2                            ', &
                                      'jso3                            ', 'jocs                            ', &
                                      'jso                             ', 'jhono                           ', &
                                      'jeuv_1                          ', 'jeuv_2                          ', &
                                      'jeuv_3                          ', 'jeuv_4                          ', &
                                      'jeuv_5                          ', 'jeuv_6                          ', &
                                      'jeuv_7                          ', 'jeuv_8                          ', &
                                      'jeuv_9                          ', 'jeuv_10                         ', &
                                      'jeuv_11                         ', 'jeuv_12                         ', &
                                      'jeuv_13                         ', 'jeuv_14                         ', &
                                      'jeuv_15                         ', 'jeuv_16                         ', &
                                      'jeuv_17                         ', 'jeuv_18                         ', &
                                      'jeuv_19                         ', 'jeuv_20                         ', &
                                      'jeuv_21                         ', 'jeuv_22                         ', &
                                      'jeuv_23                         ', 'jeuv_24                         ', &
                                      'jeuv_25                         ', 'jeuv_26                         ', &
                                      'jppi                            ', 'jepn1                           ', &
                                      'jepn2                           ', 'jepn3                           ', &
                                      'jepn4                           ', 'jepn6                           ', &
                                      'jepn7                           ', 'jpni1                           ', &
                                      'jpni2                           ', 'jpni3                           ', &
                                      'jpni4                           ', 'jpni5                           ', &
                                      'usr_O_O2                        ', 'O_O3                            ', &
                                      'usr_O_O                         ', 'O2_1S_O                         ', &
                                      'O2_1S_O2                        ', 'O2_1S_N2                        ', &
                                      'O2_1S_O3                        ', 'O2_1S_CO2                       ', &
                                      'ag2                             ', 'O2_1D_O                         ', &
                                      'O2_1D_O2                        ', 'O2_1D_N2                        ', &
                                      'ag1                             ', 'O1D_N2                          ', &
                                      'O1D_O2                          ', 'O1D_O2b                         ', &
                                      'O1D_H2O                         ', 'O1D_N2Oa                        ', &
                                      'O1D_N2Ob                        ', 'O1D_O3                          ', &
                                      'O1D_CFC11                       ', 'O1D_CFC12                       ', &
                                      'O1D_CFC113                      ', 'O1D_CFC114                      ', &
                                      'O1D_CFC115                      ', 'O1D_HCFC22                      ', &
                                      'O1D_HCFC141B                    ', 'O1D_HCFC142B                    ', &
                                      'O1D_CCL4                        ', 'O1D_CH3BR                       ', &
                                      'O1D_CF2CLBR                     ', 'O1D_CF3BR                       ', &
                                      'O1D_H1202                       ', 'O1D_H2402                       ', &
                                      'O1D_CHBR3                       ', 'O1D_CH2BR2                      ', &
                                      'O1D_COF2                        ', 'O1D_COFCL                       ', &
                                      'O1D_CH4a                        ', 'O1D_CH4b                        ', &
                                      'O1D_CH4c                        ', 'O1D_H2                          ', &
                                      'O1D_HCL                         ', 'O1D_HBR                         ', &
                                      'H_O2                            ', 'H_O3                            ', &
                                      'H_HO2a                          ', 'H_HO2                           ', &
                                      'H_HO2b                          ', 'OH_O                            ', &
                                      'OH_O3                           ', 'OH_HO2                          ', &
                                      'OH_OH                           ', 'OH_OH_M                         ', &
                                      'OH_H2                           ', 'OH_H2O2                         ', &
                                      'H2_O                            ', 'HO2_O                           ', &
                                      'HO2_O3                          ', 'usr_HO2_HO2                     ', &
                                      'H2O2_O                          ', 'HONO1                           ', &
                                      'HONO2                           ', 'N2D_O2                          ', &
                                      'N2D_O                           ', 'N_OH                            ', &
                                      'N_O2                            ', 'N_NO                            ', &
                                      'N_NO2a                          ', 'N_NO2b                          ', &
                                      'N_NO2c                          ', 'NO_O                            ', &
                                      'NO_HO2                          ', 'NO_O3                           ', &
                                      'NO2_O                           ', 'NO2_O_M                         ', &
                                      'NO2_O3                          ', 'tag_NO2_NO3                     ', &
                                      'usr_N2O5_M                      ', 'tag_NO2_OH                      ', &
                                      'usr_HNO3_OH                     ', 'NO3_NO                          ', &
                                      'NO3_O                           ', 'NO3_OH                          ', &
                                      'NO3_HO2                         ', 'tag_NO2_HO2                     ', &
                                      'HO2NO2_OH                       ', 'usr_HO2NO2_M                    ', &
                                      'CL_O3                           ', 'CL_H2                           ', &
                                      'CL_H2O2                         ', 'CL_HO2a                         ', &
                                      'CL_HO2b                         ', 'CL_CH2O                         ', &
                                      'CL_CH4                          ', 'CLO_O                           ' /)
      rxt_tag_lst(   201:   400) = (/ 'CLO_OHa                         ', 'CLO_OHb                         ', &
                                      'CLO_HO2                         ', 'CLO_CH3O2                       ', &
                                      'CLO_NO                          ', 'CLO_NO2_M                       ', &
                                      'CLO_CLOa                        ', 'CLO_CLOb                        ', &
                                      'CLO_CLOc                        ', 'tag_CLO_CLO_M                   ', &
                                      'usr_CL2O2_M                     ', 'HCL_OH                          ', &
                                      'HCL_O                           ', 'HOCL_O                          ', &
                                      'HOCL_CL                         ', 'HOCL_OH                         ', &
                                      'CLONO2_O                        ', 'CLONO2_OH                       ', &
                                      'CLONO2_CL                       ', 'BR_O3                           ', &
                                      'BR_HO2                          ', 'BR_CH2O                         ', &
                                      'BRO_O                           ', 'BRO_OH                          ', &
                                      'BRO_HO2                         ', 'BRO_NO                          ', &
                                      'BRO_NO2_M                       ', 'BRO_CLOa                        ', &
                                      'BRO_CLOb                        ', 'BRO_CLOc                        ', &
                                      'BRO_BRO                         ', 'HBR_OH                          ', &
                                      'HBR_O                           ', 'HOBR_O                          ', &
                                      'BRONO2_O                        ', 'F_H2O                           ', &
                                      'F_H2                            ', 'F_CH4                           ', &
                                      'F_HNO3                          ', 'CH3CL_CL                        ', &
                                      'CH3CL_OH                        ', 'CH3CCL3_OH                      ', &
                                      'HCFC22_OH                       ', 'CH3BR_OH                        ', &
                                      'CH3BR_CL                        ', 'HCFC141B_OH                     ', &
                                      'HCFC142B_OH                     ', 'CH2BR2_OH                       ', &
                                      'CHBR3_OH                        ', 'CH2BR2_CL                       ', &
                                      'CHBR3_CL                        ', 'CH4_OH                          ', &
                                      'usr_CO_OH_b                     ', 'CO_OH_M                         ', &
                                      'CH2O_NO3                        ', 'CH2O_OH                         ', &
                                      'CH2O_O                          ', 'CH3O2_NO                        ', &
                                      'CH3O2_HO2                       ', 'CH3OOH_OH                       ', &
                                      'usr_N2O5_aer                    ', 'usr_NO3_aer                     ', &
                                      'usr_NO2_aer                     ', 'usr_HO2_aer                     ', &
                                      'OCS_O                           ', 'OCS_OH                          ', &
                                      'S_OH                            ', 'S_O2                            ', &
                                      'S_O3                            ', 'SO_OH                           ', &
                                      'SO_O2                           ', 'SO_O3                           ', &
                                      'SO_NO2                          ', 'SO_CLO                          ', &
                                      'SO_BRO                          ', 'SO_OCLO                         ', &
                                      'usr_SO2_OH                      ', 'usr_SO3_H2O                     ', &
                                      'usr_DMS_OH                      ', 'DMS_OHb                         ', &
                                      'DMS_NO3                         ', 'het1                            ', &
                                      'het2                            ', 'het3                            ', &
                                      'het4                            ', 'het5                            ', &
                                      'het6                            ', 'het7                            ', &
                                      'het8                            ', 'het9                            ', &
                                      'het10                           ', 'het11                           ', &
                                      'het12                           ', 'het13                           ', &
                                      'het14                           ', 'het15                           ', &
                                      'het16                           ', 'het17                           ', &
                                      'ion_Op_O2                       ', 'ion_Op_N2                       ', &
                                      'ion_N2p_Oa                      ', 'ion_N2p_Ob                      ', &
                                      'ion_Op_CO2                      ', 'ion_O2p_N                       ', &
                                      'ion_O2p_NO                      ', 'ion_Np_O2a                      ', &
                                      'ion_Np_O2b                      ', 'ion_Np_O                        ', &
                                      'ion_N2p_O2                      ', 'ion_O2p_N2                      ', &
                                      'elec1                           ', 'elec2                           ', &
                                      'elec3                           ', 'ean1                            ', &
                                      'ean2                            ', 'ean3                            ', &
                                      'rpe1                            ', 'rpe2                            ', &
                                      'rpe3                            ', 'rpe4                            ', &
                                      'rpe5                            ', 'pir1                            ', &
                                      'pir2                            ', 'pir3                            ', &
                                      'pir4                            ', 'pir5                            ', &
                                      'pir6                            ', 'pir7                            ', &
                                      'pir8                            ', 'pir9                            ', &
                                      'pir10                           ', 'pir11                           ', &
                                      'pir12                           ', 'pir13                           ', &
                                      'pir14                           ', 'pir15                           ', &
                                      'pir16                           ', 'edn1                            ', &
                                      'edn2                            ', 'nir1                            ', &
                                      'nir2                            ', 'nir3                            ', &
                                      'nir4                            ', 'nir5                            ', &
                                      'nir6                            ', 'usr_CLm_H2O_M                   ', &
                                      'usr_CLm_HCL_M                   ', 'nir7                            ', &
                                      'nir8                            ', 'nir9                            ', &
                                      'nir10                           ', 'nir11                           ', &
                                      'nir12                           ', 'nir13                           ', &
                                      'iira1                           ', 'iira2                           ', &
                                      'iira3                           ', 'iira4                           ', &
                                      'iira5                           ', 'iira6                           ', &
                                      'iira7                           ', 'iira8                           ', &
                                      'iira9                           ', 'iira10                          ', &
                                      'iira11                          ', 'iira12                          ', &
                                      'iira13                          ', 'iira14                          ', &
                                      'iira15                          ', 'iira16                          ', &
                                      'iira17                          ', 'iira18                          ', &
                                      'iira19                          ', 'iira20                          ', &
                                      'iira21                          ', 'iira22                          ', &
                                      'iira23                          ', 'iira24                          ', &
                                      'iira25                          ', 'iira26                          ', &
                                      'iira27                          ', 'iira28                          ', &
                                      'iira29                          ', 'iira30                          ', &
                                      'iira31                          ', 'iira32                          ', &
                                      'iira33                          ', 'iira34                          ', &
                                      'iira35                          ', 'iira36                          ', &
                                      'iira37                          ', 'iira38                          ', &
                                      'iira39                          ', 'iira40                          ', &
                                      'iira41                          ', 'iira42                          ', &
                                      'iira43                          ', 'iira44                          ', &
                                      'iira45                          ', 'iira46                          ' /)
      rxt_tag_lst(   401:   480) = (/ 'iira47                          ', 'iira48                          ', &
                                      'iira49                          ', 'iira50                          ', &
                                      'iira51                          ', 'iira52                          ', &
                                      'iira53                          ', 'iira54                          ', &
                                      'iira55                          ', 'iira56                          ', &
                                      'iira57                          ', 'iira58                          ', &
                                      'iira59                          ', 'iira60                          ', &
                                      'iira61                          ', 'iira62                          ', &
                                      'iira63                          ', 'iira64                          ', &
                                      'iira65                          ', 'iira66                          ', &
                                      'iira67                          ', 'iira68                          ', &
                                      'iira69                          ', 'iira70                          ', &
                                      'iira71                          ', 'iira72                          ', &
                                      'iira73                          ', 'iira74                          ', &
                                      'iira75                          ', 'iira76                          ', &
                                      'iira77                          ', 'iira78                          ', &
                                      'iira79                          ', 'iira80                          ', &
                                      'iira81                          ', 'iira82                          ', &
                                      'iira83                          ', 'iira84                          ', &
                                      'iira85                          ', 'iira86                          ', &
                                      'iira87                          ', 'iira88                          ', &
                                      'iira89                          ', 'iira90                          ', &
                                      'iira91                          ', 'iira92                          ', &
                                      'iira93                          ', 'iira94                          ', &
                                      'iira95                          ', 'iira96                          ', &
                                      'iira97                          ', 'iira98                          ', &
                                      'iira99                          ', 'iira100                         ', &
                                      'iira101                         ', 'iira102                         ', &
                                      'iira103                         ', 'iira104                         ', &
                                      'iira105                         ', 'iira106                         ', &
                                      'iira107                         ', 'iira108                         ', &
                                      'iira109                         ', 'iira110                         ', &
                                      'iira111                         ', 'iira112                         ', &
                                      'iirb1                           ', 'iirb2                           ', &
                                      'iirb3                           ', 'iirb4                           ', &
                                      'iirb5                           ', 'iirb6                           ', &
                                      'iirb7                           ', 'iirb8                           ', &
                                      'iirb9                           ', 'iirb10                          ', &
                                      'iirb11                          ', 'iirb12                          ', &
                                      'iirb13                          ', 'iirb14                          ' /)
      rxt_tag_map(:rxt_tag_cnt) = (/    1,   2,   3,   4,   5,   6,   7,   8,   9,  10, &
                                       11,  12,  13,  14,  15,  16,  17,  18,  19,  20, &
                                       21,  22,  23,  24,  25,  26,  27,  28,  29,  30, &
                                       31,  32,  33,  34,  35,  36,  37,  38,  39,  40, &
                                       41,  42,  43,  44,  45,  46,  47,  48,  49,  50, &
                                       51,  52,  53,  54,  55,  56,  57,  58,  59,  60, &
                                       61,  62,  63,  64,  65,  66,  67,  68,  69,  70, &
                                       71,  72,  73,  74,  75,  76,  77,  78,  79,  80, &
                                       81,  82,  83,  84,  85,  86,  87,  88,  89,  90, &
                                       91,  92,  93,  94,  95,  96,  97,  98,  99, 100, &
                                      101, 102, 103, 104, 105, 106, 107, 108, 109, 110, &
                                      111, 112, 113, 114, 115, 116, 117, 118, 119, 120, &
                                      121, 122, 123, 124, 125, 126, 127, 128, 129, 130, &
                                      131, 132, 133, 134, 135, 136, 137, 138, 139, 140, &
                                      141, 142, 143, 144, 145, 146, 147, 148, 149, 150, &
                                      151, 152, 153, 154, 155, 156, 157, 158, 159, 160, &
                                      161, 162, 163, 164, 165, 166, 167, 168, 169, 170, &
                                      171, 172, 173, 174, 175, 176, 177, 178, 179, 180, &
                                      181, 182, 183, 184, 185, 186, 187, 188, 189, 190, &
                                      191, 192, 193, 194, 195, 196, 197, 198, 199, 200, &
                                      201, 202, 203, 204, 205, 206, 207, 208, 209, 210, &
                                      211, 212, 213, 214, 215, 216, 217, 218, 219, 220, &
                                      221, 222, 223, 224, 225, 226, 227, 228, 229, 230, &
                                      231, 232, 233, 234, 235, 236, 237, 238, 239, 240, &
                                      241, 242, 243, 244, 245, 246, 247, 248, 249, 250, &
                                      251, 252, 253, 254, 255, 256, 257, 258, 259, 260, &
                                      261, 262, 263, 264, 265, 266, 267, 268, 269, 270, &
                                      271, 272, 273, 274, 275, 276, 277, 278, 279, 280, &
                                      281, 282, 283, 284, 285, 286, 287, 288, 289, 290, &
                                      291, 292, 293, 294, 295, 296, 297, 298, 299, 300, &
                                      301, 302, 303, 304, 305, 306, 307, 308, 309, 310, &
                                      311, 312, 313, 314, 315, 316, 317, 324, 325, 326, &
                                      329, 332, 340, 341, 342, 343, 344, 345, 346, 347, &
                                      352, 353, 357, 360, 361, 362, 363, 365, 371, 399, &
                                      400, 401, 402, 414, 415, 445, 446, 449, 452, 455, &
                                      457, 460, 462, 465, 468, 469, 470, 471, 472, 473, &
                                      474, 475, 476, 477, 478, 479, 480, 481, 482, 483, &
                                      484, 485, 486, 487, 488, 489, 490, 491, 492, 493, &
                                      494, 495, 496, 497, 498, 499, 500, 501, 502, 503, &
                                      504, 505, 506, 507, 508, 509, 510, 511, 512, 513, &
                                      514, 515, 516, 517, 518, 519, 520, 521, 522, 523, &
                                      524, 525, 526, 527, 528, 529, 530, 531, 532, 533, &
                                      534, 535, 536, 537, 538, 539, 540, 541, 542, 543, &
                                      544, 545, 546, 547, 548, 549, 550, 551, 552, 553, &
                                      554, 555, 556, 557, 558, 559, 560, 561, 562, 563, &
                                      564, 565, 566, 567, 568, 569, 570, 571, 572, 573, &
                                      574, 575, 576, 577, 578, 579, 580, 581, 582, 583, &
                                      584, 585, 586, 587, 588, 589, 590, 591, 592, 593 /)
      if( allocated( pht_alias_lst ) ) then
         deallocate( pht_alias_lst )
      end if
      allocate( pht_alias_lst(phtcnt,2),stat=ios )
      if( ios /= 0 ) then
         write(iulog,*) 'set_sim_dat: failed to allocate pht_alias_lst; error = ',ios
         call endrun
      end if
      if( allocated( pht_alias_mult ) ) then
         deallocate( pht_alias_mult )
      end if
      allocate( pht_alias_mult(phtcnt,2),stat=ios )
      if( ios /= 0 ) then
         write(iulog,*) 'set_sim_dat: failed to allocate pht_alias_mult; error = ',ios
         call endrun
      end if
      pht_alias_lst(:,1) = (/ 'userdefined     ', 'userdefined     ', '                ', '                ', &
                              '                ', 'userdefined     ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', 'userdefined     ', 'userdefined     ', &
                              'userdefined     ', 'userdefined     ', 'userdefined     ', 'userdefined     ', &
                              'userdefined     ', 'userdefined     ', 'userdefined     ', 'userdefined     ', &
                              'userdefined     ', 'userdefined     ', 'userdefined     ', 'userdefined     ', &
                              'userdefined     ', 'userdefined     ', 'userdefined     ', 'userdefined     ', &
                              'userdefined     ', 'userdefined     ', 'userdefined     ', 'userdefined     ', &
                              'userdefined     ', 'userdefined     ', 'userdefined     ', 'userdefined     ', &
                              'userdefined     ', 'userdefined     ', 'userdefined     ', 'userdefined     ', &
                              'userdefined     ', 'userdefined     ', 'userdefined     ', 'userdefined     ', &
                              'userdefined     ', 'userdefined     ', 'userdefined     ', 'userdefined     ' /)
      pht_alias_lst(:,2) = (/ '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', 'userdefined     ', 'userdefined     ', &
                              'userdefined     ', 'userdefined     ', 'userdefined     ', 'userdefined     ', &
                              'userdefined     ', 'userdefined     ', 'userdefined     ', 'userdefined     ', &
                              'userdefined     ', 'userdefined     ', 'userdefined     ', 'userdefined     ', &
                              'userdefined     ', 'userdefined     ', 'userdefined     ', 'userdefined     ', &
                              'userdefined     ', 'userdefined     ', 'userdefined     ', 'userdefined     ', &
                              'userdefined     ', 'userdefined     ', 'userdefined     ', 'userdefined     ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ' /)
      pht_alias_mult(:,1) = (/ 1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8 /)
      pht_alias_mult(:,2) = (/ 1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8 /)
      allocate( cph_enthalpy(enthalpy_cnt),stat=ios )
      if( ios /= 0 ) then
         write(iulog,*) 'set_sim_dat: failed to allocate cph_enthalpy; error = ',ios
         call endrun
      end if
      allocate( cph_rid(enthalpy_cnt),stat=ios )
      if( ios /= 0 ) then
         write(iulog,*) 'set_sim_dat: failed to allocate cph_rid; error = ',ios
         call endrun
      end if
      cph_rid(:)      = (/             105,            106,            107,            108,            109, &
                                       110,            111,            114,            115,            116, &
                                       118,            119,            120,            149,            150, &
                                       152,            154,            155,            156,            162, &
                                       163,            164,            168,            169,            171, &
                                       172,            177,            178,            179,            299, &
                                       300,            301,            304,            305,            306, &
                                       307,            308,            309,            311,            312, &
                                       313 /)
      cph_enthalpy(:) = (/   101.390000_r8,  392.190000_r8,  493.580000_r8,   62.600000_r8,   62.600000_r8, &
                              62.600000_r8,   62.600000_r8,   94.300000_r8,   94.300000_r8,   94.300000_r8, &
                             189.910000_r8,   32.910000_r8,  189.810000_r8,  203.400000_r8,  194.710000_r8, &
                             232.590000_r8,   67.670000_r8,  165.300000_r8,  293.620000_r8,  226.580000_r8, &
                             120.100000_r8,  165.510000_r8,  177.510000_r8,  229.610000_r8,  133.750000_r8, &
                             313.750000_r8,   34.470000_r8,  199.170000_r8,  193.020000_r8,  150.110000_r8, &
                             105.040000_r8,   67.530000_r8,  406.160000_r8,  271.380000_r8,  239.840000_r8, &
                             646.280000_r8,   95.550000_r8,  339.590000_r8,   82.389000_r8,  508.950000_r8, &
                             354.830000_r8 /)
      allocate( num_rnts(rxntot-phtcnt),stat=ios )
      if( ios /= 0 ) then
         write(iulog,*) 'set_sim_dat: failed to allocate num_rnts; error = ',ios
         call endrun
      end if
      num_rnts(:) = (/      3,     2,     3,     2,     2,     2,     2,     2,     1,     2, &
                            2,     2,     1,     2,     2,     2,     2,     2,     2,     2, &
                            2,     2,     2,     2,     2,     2,     2,     2,     2,     2, &
                            2,     2,     2,     2,     2,     2,     2,     2,     2,     2, &
                            2,     2,     2,     2,     3,     2,     2,     2,     2,     2, &
                            2,     2,     2,     3,     2,     2,     2,     2,     2,     2, &
                            2,     3,     2,     2,     2,     2,     2,     2,     2,     2, &
                            2,     3,     2,     2,     2,     3,     2,     3,     2,     3, &
                            2,     2,     2,     2,     2,     3,     2,     2,     2,     2, &
                            2,     2,     2,     2,     2,     2,     2,     2,     2,     2, &
                            2,     3,     2,     2,     2,     3,     2,     2,     2,     2, &
                            2,     2,     2,     2,     2,     2,     2,     2,     2,     2, &
                            2,     2,     3,     2,     2,     2,     2,     2,     2,     2, &
                            2,     2,     2,     2,     2,     2,     2,     2,     2,     2, &
                            2,     2,     2,     2,     2,     2,     2,     2,     2,     3, &
                            2,     2,     2,     2,     2,     2,     1,     1,     1,     1, &
                            2,     2,     2,     2,     2,     2,     2,     2,     2,     2, &
                            2,     2,     2,     2,     2,     2,     2,     1,     1,     1, &
                            2,     2,     2,     1,     1,     2,     2,     1,     1,     1, &
                            1,     2,     2,     2,     2,     2,     2,     2,     2,     2, &
                            2,     2,     2,     2,     2,     2,     2,     2,     2,     3, &
                            2,     3,     2,     2,     2,     2,     2,     2,     2,     2, &
                            2,     2,     2,     2,     2,     2,     2,     3,     3,     2, &
                            2,     2,     2,     2,     2,     3,     2,     3,     2,     3, &
                            2,     3,     2,     2,     2,     2,     2,     3,     3,     2, &
                            2,     2,     3,     2,     2,     3,     2,     3,     2,     2, &
                            2,     2,     2,     2,     2,     2,     2,     2,     2,     2, &
                            2,     2,     2,     3,     2,     2,     3,     2,     2,     2, &
                            2,     2,     2,     2,     3,     2,     3,     2,     2,     2, &
                            2,     2,     2,     2,     2,     2,     2,     2,     2,     2, &
                            2,     2,     3,     2,     2,     2,     2,     2,     2,     2, &
                            2,     2,     2,     2,     2,     2,     2,     2,     2,     2, &
                            2,     2,     2,     2,     2,     2,     2,     2,     2,     2, &
                            2,     2,     2,     2,     3,     3,     2,     2,     2,     2, &
                            2,     2,     3,     3,     2,     2,     2,     2,     2,     3, &
                            2,     2,     2,     3,     3,     2,     2,     2,     3,     2, &
                            3,     2,     2,     2,     2,     2,     2,     2,     2,     2, &
                            2,     2,     2,     2,     2,     2,     2,     2,     2,     2, &
                            2,     2,     2,     2,     2,     2,     2,     2,     2,     2, &
                            2,     2,     2,     2,     2,     2,     2,     2,     2,     2, &
                            2,     2,     2,     2,     2,     2,     2,     2,     2,     2, &
                            2,     2,     2,     2,     2,     2,     2,     2,     2,     2, &
                            2,     2,     2,     2,     2,     2,     2,     2,     2,     2, &
                            2,     2,     2,     2,     2,     2,     2,     2,     2,     2, &
                            2,     2,     2,     2,     2,     2,     2,     2,     2,     2, &
                            2,     2,     2,     2,     2,     2,     2,     2,     2,     2, &
                            2,     2,     2,     2,     2,     2,     2,     2,     2,     2, &
                            2,     2,     2,     2,     2,     3,     3,     3,     3,     3, &
                            3,     3,     3,     3,     3,     3,     3,     3,     3 /)

      end subroutine set_sim_dat

      end module mo_sim_dat
