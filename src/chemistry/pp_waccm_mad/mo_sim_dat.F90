
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

      is_scalar = .false.
      is_vector = .true.

      clscnt(:) = (/     22,     0,     0,    88,     0 /)

      cls_rxt_cnt(:,1) = (/      1,    56,     0,    22 /)
      cls_rxt_cnt(:,4) = (/     28,   151,   405,    88 /)

      solsym(:110) = (/ 'BRCL            ','BRO             ','BRONO2          ','BRY             ','CCL4            ', &
                        'CF2CLBR         ','CF3BR           ','CFC11           ','CFC113          ','CFC114          ', &
                        'CFC115          ','CFC12           ','CH2BR2          ','CH2O            ','CH3BR           ', &
                        'CH3CCL3         ','CH3CL           ','CH3O2           ','CH3OOH          ','CH4             ', &
                        'CHBR3           ','CL2             ','CL2O2           ','CLO             ','CLONO2          ', &
                        'CLY             ','CO              ','CO2             ','COF2            ','COFCL           ', &
                        'F               ','H               ','H2              ','H2402           ','H2O2            ', &
                        'HBR             ','HCFC141B        ','HCFC142B        ','HCFC22          ','HCL             ', &
                        'HF              ','HNO3            ','HO2             ','HO2NO2          ','HOBR            ', &
                        'HOCL            ','HONO            ','N               ','N2O             ','N2O5            ', &
                        'NO              ','NO2             ','NO3             ','O               ','O2              ', &
                        'O3              ','OCLO            ','SF6             ','BR              ','CL              ', &
                        'CLm             ','CLm_H2O         ','CLm_HCL         ','CLOm            ','CO3m            ', &
                        'CO3m2H2O        ','CO3m_H2O        ','CO4m            ','e               ','H3Op_OH         ', &
                        'HCO3m           ','Hp_2H2O         ','Hp_3H2O         ','Hp_3N1          ','Hp_4H2O         ', &
                        'Hp_4N1          ','Hp_5H2O         ','Hp_H2O          ','N2D             ','N2p             ', &
                        'NO2m            ','NO2m_H2O        ','NO3m            ','NO3m2H2O        ','NO3m_H2O        ', &
                        'NO3m_HCL        ','NO3mHNO3        ','NOp             ','NOp_2H2O        ','NOp_3H2O        ', &
                        'NOp_CO2         ','NOp_H2O         ','NOp_N2          ','Np              ','O1D             ', &
                        'O2_1D           ','O2_1S           ','O2m             ','O2p             ','O2p_H2O         ', &
                        'O3m             ','O4m             ','O4p             ','OH              ','OHm             ', &
                        'Om              ','Op              ','Op2D            ','Op2P            ','H2O             ' /)

      adv_mass(:110) = (/   115.356700_r8,    95.903400_r8,   141.908940_r8,    99.716850_r8,   153.821800_r8, &
                            165.364506_r8,   148.910210_r8,   137.367503_r8,   187.375310_r8,   170.921013_r8, &
                            154.466716_r8,   120.913206_r8,   173.833800_r8,    30.025200_r8,    94.937200_r8, &
                            133.402300_r8,    50.485900_r8,    47.032000_r8,    48.039400_r8,    16.040600_r8, &
                            252.730400_r8,    70.905400_r8,   102.904200_r8,    51.452100_r8,    97.457640_r8, &
                            100.916850_r8,    28.010400_r8,    44.009800_r8,    66.007206_r8,    82.461503_r8, &
                             18.998403_r8,     1.007400_r8,     2.014800_r8,   259.823613_r8,    34.013600_r8, &
                             80.911400_r8,   116.948003_r8,   100.493706_r8,    86.467906_r8,    36.460100_r8, &
                             20.005803_r8,    63.012340_r8,    33.006200_r8,    79.011740_r8,    96.910800_r8, &
                             52.459500_r8,    47.012940_r8,    14.006740_r8,    44.012880_r8,   108.010480_r8, &
                             30.006140_r8,    46.005540_r8,    62.004940_r8,    15.999400_r8,    31.998800_r8, &
                             47.998200_r8,    67.451500_r8,   146.056419_r8,    79.904000_r8,    35.452700_r8, &
                             35.452700_r8,    53.466900_r8,    71.912800_r8,    51.452100_r8,    60.009200_r8, &
                             96.037600_r8,    78.023400_r8,    76.008600_r8, 0.548567E-03_r8,    36.028400_r8, &
                             61.016600_r8,    37.035800_r8,    55.050000_r8,   118.062340_r8,    73.064200_r8, &
                            136.076540_r8,    91.078400_r8,    19.021600_r8,    14.006740_r8,    28.013480_r8, &
                             46.005540_r8,    64.019740_r8,    62.004940_r8,    98.033340_r8,    80.019140_r8, &
                             98.465040_r8,   125.017280_r8,    30.006140_r8,    66.034540_r8,    68.049340_r8, &
                             74.015940_r8,    48.020340_r8,    58.019620_r8,    14.006740_r8,    15.999400_r8, &
                             31.998800_r8,    31.998800_r8,    31.998800_r8,    31.998800_r8,    50.013000_r8, &
                             47.998200_r8,    63.997600_r8,    63.997600_r8,    17.006800_r8,    17.006800_r8, &
                             15.999400_r8,    15.999400_r8,    15.999400_r8,    15.999400_r8,    18.014200_r8 /)

      crb_mass(:110) = (/     0.000000_r8,     0.000000_r8,     0.000000_r8,     0.000000_r8,    12.011000_r8, &
                             12.011000_r8,    12.011000_r8,    12.011000_r8,    24.022000_r8,    24.022000_r8, &
                             24.022000_r8,    12.011000_r8,    12.011000_r8,    12.011000_r8,    12.011000_r8, &
                             24.022000_r8,    12.011000_r8,    12.011000_r8,    12.011000_r8,    12.011000_r8, &
                             12.011000_r8,     0.000000_r8,     0.000000_r8,     0.000000_r8,     0.000000_r8, &
                             12.011000_r8,    12.011000_r8,    12.011000_r8,    12.011000_r8,    12.011000_r8, &
                              0.000000_r8,     0.000000_r8,     0.000000_r8,    24.022000_r8,     0.000000_r8, &
                              0.000000_r8,    24.022000_r8,    24.022000_r8,    12.011000_r8,     0.000000_r8, &
                              0.000000_r8,     0.000000_r8,     0.000000_r8,     0.000000_r8,     0.000000_r8, &
                              0.000000_r8,     0.000000_r8,     0.000000_r8,     0.000000_r8,     0.000000_r8, &
                              0.000000_r8,     0.000000_r8,     0.000000_r8,     0.000000_r8,     0.000000_r8, &
                              0.000000_r8,     0.000000_r8,     0.000000_r8,     0.000000_r8,     0.000000_r8, &
                              0.000000_r8,     0.000000_r8,     0.000000_r8,     0.000000_r8,    12.011000_r8, &
                             12.011000_r8,    12.011000_r8,    12.011000_r8,     0.000000_r8,     0.000000_r8, &
                             12.011000_r8,     0.000000_r8,     0.000000_r8,     0.000000_r8,     0.000000_r8, &
                              0.000000_r8,     0.000000_r8,     0.000000_r8,     0.000000_r8,     0.000000_r8, &
                              0.000000_r8,     0.000000_r8,     0.000000_r8,     0.000000_r8,     0.000000_r8, &
                              0.000000_r8,     0.000000_r8,     0.000000_r8,     0.000000_r8,     0.000000_r8, &
                             12.011000_r8,     0.000000_r8,     0.000000_r8,     0.000000_r8,     0.000000_r8, &
                              0.000000_r8,     0.000000_r8,     0.000000_r8,     0.000000_r8,     0.000000_r8, &
                              0.000000_r8,     0.000000_r8,     0.000000_r8,     0.000000_r8,     0.000000_r8, &
                              0.000000_r8,     0.000000_r8,     0.000000_r8,     0.000000_r8,     0.000000_r8 /)

      fix_mass(:  2) = (/ 0.00000000_r8, 28.0134800_r8 /)

      clsmap(: 22,1) = (/    4,   5,   6,   7,   8,   9,  10,  11,  12,  13, &
                            15,  16,  17,  20,  21,  26,  34,  37,  38,  39, &
                            49,  58 /)
      clsmap(: 88,4) = (/    1,   2,   3,  14,  18,  19,  22,  23,  24,  25, &
                            27,  29,  30,  28,  31,  32,  33,  35,  36,  40, &
                            41,  42,  43,  44,  45,  46,  47,  48,  50,  51, &
                            52,  53,  54,  55,  56,  57,  59,  60,  61,  62, &
                            63,  64,  65,  66,  67,  68,  69,  70,  71,  72, &
                            73,  74,  75,  76,  77,  78,  79,  80,  81,  82, &
                            83,  84,  85,  86,  87,  88,  89,  90,  91,  92, &
                            93,  94,  95,  96,  97,  98,  99, 100, 101, 102, &
                           103, 104, 105, 106, 107, 108, 109, 110 /)

      permute(: 88,4) = (/    6,  52,  17,  56,  39,  12,   4,   1,  81,  34, &
                             25,   2,   5,  78,  30,  86,  50,  18,  33,  80, &
                              7,  71,  61,  14,  28,  32,  24,  40,  37,  87, &
                             72,  74,  68,  75,  79,   3,  44,  85,  58,  46, &
                             45,  35,  64,  43,  54,  59,  69,  11,  42,  21, &
                             65,   8,  73,   9,  82,  26,  29,  19,  62,  47, &
                             60,  53,  55,  48,  51,  66,  67,  13,  23,  70, &
                             22,  31,  57,  36,  10,  76,  77,  27,  49,  20, &
                             41,  63,  83,  84,  38,  16,  15,  88 /)

      diag_map(: 88) = (/    1,   4,   7,  10,  12,  16,  19,  22,  26,  30, &
                            33,  39,  45,  51,  58,  65,  71,  79,  86,  95, &
                           102, 108, 116, 123, 128, 134, 142, 151, 159, 168, &
                           176, 184, 191, 201, 211, 219, 229, 244, 257, 270, &
                           285, 297, 311, 328, 340, 356, 371, 385, 401, 419, &
                           438, 457, 477, 497, 519, 541, 573, 597, 619, 652, &
                           686, 712, 755, 785, 829, 878, 922, 979,1023,1066, &
                          1107,1153,1197,1233,1277,1314,1360,1403,1442,1488, &
                          1532,1576,1610,1647,1691,1728,1771,1825 /)

      extfrc_lst(: 12) = (/ 'NO2             ','NO              ','CO              ','N2p             ','Np              ', &
                            'O               ','O2p             ','OH              ','Op              ','e               ', &
                            'N               ','N2D             ' /)

      frc_from_dataset(: 12) = (/ .true., .true., .true., .false., .false., &
                                  .false., .false., .false., .false., .false., &
                                  .false., .false. /)

      inv_lst(:  2) = (/ 'M               ', 'N2              ' /)

      slvd_lst(: 51) = (/ 'BR              ', 'CL              ', 'CLm             ', 'CLm_H2O         ', 'CLm_HCL         ', &
                          'CLOm            ', 'CO3m            ', 'CO3m2H2O        ', 'CO3m_H2O        ', 'CO4m            ', &
                          'e               ', 'H3Op_OH         ', 'HCO3m           ', 'Hp_2H2O         ', 'Hp_3H2O         ', &
                          'Hp_3N1          ', 'Hp_4H2O         ', 'Hp_4N1          ', 'Hp_5H2O         ', 'Hp_H2O          ', &
                          'N2D             ', 'N2p             ', 'NO2m            ', 'NO2m_H2O        ', 'NO3m            ', &
                          'NO3m2H2O        ', 'NO3m_H2O        ', 'NO3m_HCL        ', 'NO3mHNO3        ', 'NOp             ', &
                          'NOp_2H2O        ', 'NOp_3H2O        ', 'NOp_CO2         ', 'NOp_H2O         ', 'NOp_N2          ', &
                          'Np              ', 'O1D             ', 'O2_1D           ', 'O2_1S           ', 'O2m             ', &
                          'O2p             ', 'O2p_H2O         ', 'O3m             ', 'O4m             ', 'O4p             ', &
                          'OH              ', 'OHm             ', 'Om              ', 'Op              ', 'Op2D            ', &
                          'Op2P            ' /)

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
      rxt_tag_lst(     1:   200) = (/ 'jh2o_b                          ', 'jh2o_c                          ', &
                                      'jh2o_a                          ', 'jh2o2                           ', &
                                      'jo2_b                           ', 'jo2_a                           ', &
                                      'jo3_a                           ', 'jo3_b                           ', &
                                      'jhno3                           ', 'jho2no2_a                       ', &
                                      'jho2no2_b                       ', 'jhono                           ', &
                                      'jn2o                            ', 'jn2o5_a                         ', &
                                      'jn2o5_b                         ', 'jno_i                           ', &
                                      'jno                             ', 'jno2                            ', &
                                      'jno3_a                          ', 'jno3_b                          ', &
                                      'jch2o_a                         ', 'jch2o_b                         ', &
                                      'jch3ooh                         ', 'jch4_a                          ', &
                                      'jch4_b                          ', 'jco2                            ', &
                                      'jbrcl                           ', 'jbro                            ', &
                                      'jbrono2_b                       ', 'jbrono2_a                       ', &
                                      'jccl4                           ', 'jcf2clbr                        ', &
                                      'jcf3br                          ', 'jcfcl3                          ', &
                                      'jcfc113                         ', 'jcfc114                         ', &
                                      'jcfc115                         ', 'jcf2cl2                         ', &
                                      'jch2br2                         ', 'jch3br                          ', &
                                      'jch3ccl3                        ', 'jch3cl                          ', &
                                      'jchbr3                          ', 'jcl2                            ', &
                                      'jcl2o2                          ', 'jclo                            ', &
                                      'jclono2_a                       ', 'jclono2_b                       ', &
                                      'jcof2                           ', 'jcofcl                          ', &
                                      'jh2402                          ', 'jhbr                            ', &
                                      'jhcfc141b                       ', 'jhcfc142b                       ', &
                                      'jhcfc22                         ', 'jhcl                            ', &
                                      'jhf                             ', 'jhobr                           ', &
                                      'jhocl                           ', 'joclo                           ', &
                                      'jsf6                            ', 'jeuv_26                         ', &
                                      'jpni3                           ', 'jpni5                           ', &
                                      'jpni4                           ', 'jeuv_4                          ', &
                                      'jeuv_25                         ', 'jeuv_18                         ', &
                                      'jeuv_13                         ', 'jeuv_23                         ', &
                                      'jeuv_11                         ', 'jeuv_6                          ', &
                                      'jeuv_10                         ', 'jeuv_22                         ', &
                                      'jepn6                           ', 'jepn7                           ', &
                                      'jeuv_16                         ', 'jeuv_1                          ', &
                                      'jeuv_14                         ', 'jeuv_15                         ', &
                                      'jeuv_2                          ', 'jeuv_3                          ', &
                                      'jeuv_21                         ', 'jeuv_24                         ', &
                                      'jeuv_19                         ', 'jeuv_5                          ', &
                                      'jeuv_7                          ', 'jeuv_17                         ', &
                                      'jeuv_12                         ', 'jeuv_20                         ', &
                                      'jeuv_8                          ', 'jeuv_9                          ', &
                                      'jepn2                           ', 'jppi                            ', &
                                      'jpni1                           ', 'jepn3                           ', &
                                      'jpni2                           ', 'jepn4                           ', &
                                      'jepn1                           ', 'CLm_H                           ', &
                                      'CLmH2O_HCL                      ', 'CLm_H2O_Ma                      ', &
                                      'CLmHCL_M                        ', 'CLm_HNO3                        ', &
                                      'CLm_NO2                         ', 'CLOm_NOa                        ', &
                                      'CLOm_NOb                        ', 'CLOm_O                          ', &
                                      'CO3m_CLa                        ', 'CO3m_CLb                        ', &
                                      'CO3m_CLO                        ', 'CO3m_H                          ', &
                                      'CO3mH2O_H2O_M                   ', 'CO3m_H2O_M                      ', &
                                      'CO3mH2O_NO2a                    ', 'CO3mH2O_NO2b                    ', &
                                      'CO3mH2O_NOa                     ', 'CO3mH2O_NOb                     ', &
                                      'CO3m_HNO3                       ', 'CO3m_O                          ', &
                                      'CO3m_O2                         ', 'CO4m_CL                         ', &
                                      'CO4m_CLO                        ', 'CO4m_H                          ', &
                                      'CO4m_HCL                        ', 'CO4m_O                          ', &
                                      'CO4m_O3                         ', 'ean1                            ', &
                                      'ean2                            ', 'ean3                            ', &
                                      'edn1                            ', 'edn2                            ', &
                                      'H3OpOH_e                        ', 'H3OpOH_H2O                      ', &
                                      'Hp3N1_H2O                       ', 'Hp4H2O_e                        ', &
                                      'Hp4H2O_N2O5                     ', 'Hp4N1_H2O                       ', &
                                      'Hp5H2O_e                        ', 'Hp5H2O_N2O5                     ', &
                                      'iira1                           ', 'iira10                          ', &
                                      'iira100                         ', 'iira101                         ', &
                                      'iira102                         ', 'iira103                         ', &
                                      'iira104                         ', 'iira105                         ', &
                                      'iira106                         ', 'iira107                         ', &
                                      'iira108                         ', 'iira109                         ', &
                                      'iira11                          ', 'iira110                         ', &
                                      'iira111                         ', 'iira112                         ', &
                                      'iira12                          ', 'iira13                          ', &
                                      'iira14                          ', 'iira15                          ', &
                                      'iira16                          ', 'iira17                          ', &
                                      'iira18                          ', 'iira19                          ', &
                                      'iira2                           ', 'iira20                          ', &
                                      'iira21                          ', 'iira22                          ', &
                                      'iira23                          ', 'iira24                          ', &
                                      'iira25                          ', 'iira26                          ', &
                                      'iira27                          ', 'iira28                          ', &
                                      'iira29                          ', 'iira3                           ', &
                                      'iira30                          ', 'iira31                          ', &
                                      'iira32                          ', 'iira33                          ', &
                                      'iira34                          ', 'iira35                          ', &
                                      'iira36                          ', 'iira37                          ', &
                                      'iira38                          ', 'iira39                          ', &
                                      'iira4                           ', 'iira40                          ', &
                                      'iira41                          ', 'iira42                          ', &
                                      'iira43                          ', 'iira44                          ', &
                                      'iira45                          ', 'iira46                          ', &
                                      'iira47                          ', 'iira48                          ', &
                                      'iira49                          ', 'iira5                           ', &
                                      'iira50                          ', 'iira51                          ' /)
      rxt_tag_lst(   201:   400) = (/ 'iira52                          ', 'iira53                          ', &
                                      'iira54                          ', 'iira55                          ', &
                                      'iira56                          ', 'iira57                          ', &
                                      'iira58                          ', 'iira59                          ', &
                                      'iira6                           ', 'iira60                          ', &
                                      'iira61                          ', 'iira62                          ', &
                                      'iira63                          ', 'iira64                          ', &
                                      'iira65                          ', 'iira66                          ', &
                                      'iira67                          ', 'iira68                          ', &
                                      'iira69                          ', 'iira7                           ', &
                                      'iira70                          ', 'iira71                          ', &
                                      'iira72                          ', 'iira73                          ', &
                                      'iira74                          ', 'iira75                          ', &
                                      'iira76                          ', 'iira77                          ', &
                                      'iira78                          ', 'iira79                          ', &
                                      'iira8                           ', 'iira80                          ', &
                                      'iira81                          ', 'iira82                          ', &
                                      'iira83                          ', 'iira84                          ', &
                                      'iira85                          ', 'iira86                          ', &
                                      'iira87                          ', 'iira88                          ', &
                                      'iira89                          ', 'iira9                           ', &
                                      'iira90                          ', 'iira91                          ', &
                                      'iira92                          ', 'iira93                          ', &
                                      'iira94                          ', 'iira95                          ', &
                                      'iira96                          ', 'iira97                          ', &
                                      'iira98                          ', 'iira99                          ', &
                                      'iirb1                           ', 'iirb10                          ', &
                                      'iirb11                          ', 'iirb12                          ', &
                                      'iirb13                          ', 'iirb14                          ', &
                                      'iirb2                           ', 'iirb3                           ', &
                                      'iirb4                           ', 'iirb5                           ', &
                                      'iirb6                           ', 'iirb7                           ', &
                                      'iirb8                           ', 'iirb9                           ', &
                                      'nir1                            ', 'nir10                           ', &
                                      'nir11                           ', 'nir12                           ', &
                                      'nir13                           ', 'nir2                            ', &
                                      'nir3                            ', 'nir4                            ', &
                                      'nir5                            ', 'nir6                            ', &
                                      'nir7                            ', 'nir8                            ', &
                                      'nir9                            ', 'NO2m_CL                         ', &
                                      'NO2m_CLO                        ', 'NO2m_H                          ', &
                                      'NO2m_H2O_M                      ', 'NO2m_HCL                        ', &
                                      'NO2m_HNO3                       ', 'NO2m_NO2                        ', &
                                      'NO2m_O3                         ', 'NO3m2H2O_N2O5                   ', &
                                      'NO3mH2O_H2O_M                   ', 'NO3mH2O_HNO3                    ', &
                                      'NO3m_H2O_M                      ', 'NO3mH2O_N2O5                    ', &
                                      'NO3m_HCLa                       ', 'NO3mHCL_HNO3                    ', &
                                      'NO3m_HNO3_M                     ', 'NO3m_O                          ', &
                                      'NO3m_O3                         ', 'NOp2H2O_e                       ', &
                                      'NOp3H2O_e                       ', 'NOp3H2O_H2O                     ', &
                                      'NOpCO2_e                        ', 'NOpCO2_H2O                      ', &
                                      'NOpH2O_e                        ', 'NOpH2O_H                        ', &
                                      'NOpH2O_HO2                      ', 'NOpH2O_OH                       ', &
                                      'NOpN2_CO2                       ', 'NOpN2_H2O                       ', &
                                      'O2m_CL                          ', 'O2m_CLO                         ', &
                                      'O2m_CO2_M                       ', 'O2m_H                           ', &
                                      'O2m_HCL                         ', 'O2m_HNO3                        ', &
                                      'O2m_NO2                         ', 'O2m_O21D                        ', &
                                      'O2m_O2_M                        ', 'O2m_O3                          ', &
                                      'O2m_O_a                         ', 'O2m_O_b                         ', &
                                      'O2pH2O_e                        ', 'O2pH2O_H2Oa                     ', &
                                      'O2pH2O_H2Ob                     ', 'O2p_H2O_M                       ', &
                                      'O3m_CO2                         ', 'O3m_H                           ', &
                                      'O3m_O3                          ', 'O3m_O_a                         ', &
                                      'O3m_O_b                         ', 'O4m_CO2                         ', &
                                      'O4m_O                           ', 'O4p_H2O                         ', &
                                      'O4p_O                           ', 'O4p_O21D                        ', &
                                      'OH_HONO                         ', 'OHm_CL                          ', &
                                      'OHm_CLO                         ', 'OHm_CO2                         ', &
                                      'OHm_H                           ', 'OHm_HCL                         ', &
                                      'OHm_NO2                         ', 'OHm_O                           ', &
                                      'OHm_O3                          ', 'OH_NO_M                         ', &
                                      'Om_CL                           ', 'Om_CLO                          ', &
                                      'Om_CO2_M                        ', 'Om_H2_a                         ', &
                                      'Om_H2_b                         ', 'Om_H2O                          ', &
                                      'Om_HCL                          ', 'Om_HNO3                         ', &
                                      'Om_M                            ', 'Om_NO2                          ', &
                                      'Om_O                            ', 'Om_O21D                         ', &
                                      'Om_O2_M                         ', 'Om_O3                           ', &
                                      'pir1                            ', 'pir10                           ', &
                                      'pir11                           ', 'pir12                           ', &
                                      'pir13                           ', 'pir14                           ', &
                                      'pir15                           ', 'pir16                           ', &
                                      'pir2                            ', 'pir3                            ', &
                                      'pir4                            ', 'pir5                            ', &
                                      'pir6                            ', 'pir7                            ', &
                                      'pir8                            ', 'pir9                            ', &
                                      'rpe1                            ', 'rpe2                            ', &
                                      'rpe3                            ', 'rpe4                            ', &
                                      'rpe5                            ', 'usr_CLm_H2O_M                   ', &
                                      'usr_CLm_HCL_M                   ', 'ag1                             ', &
                                      'ag2                             ', 'O1D_H2                          ', &
                                      'O1D_H2O                         ', 'O1D_N2                          ', &
                                      'O1D_O2                          ', 'O1D_O2b                         ', &
                                      'O1D_O3                          ', 'O2_1D_N2                        ', &
                                      'O2_1D_O                         ', 'O2_1D_O2                        ', &
                                      'O2_1S_CO2                       ', 'O2_1S_N2                        ', &
                                      'O2_1S_O                         ', 'O2_1S_O2                        ', &
                                      'O2_1S_O3                        ', 'O_O3                            ', &
                                      'usr_O_O                         ', 'usr_O_O2                        ' /)
      rxt_tag_lst(   401:   585) = (/ 'H2_O                            ', 'H2O2_O                          ', &
                                      'H_HO2                           ', 'H_HO2a                          ', &
                                      'H_HO2b                          ', 'H_O2                            ', &
                                      'HO2_O                           ', 'HO2_O3                          ', &
                                      'H_O3                            ', 'OH_H2                           ', &
                                      'OH_H2O2                         ', 'OH_HO2                          ', &
                                      'OH_O                            ', 'OH_O3                           ', &
                                      'OH_OH                           ', 'OH_OH_M                         ', &
                                      'usr_HO2_HO2                     ', 'HO2NO2_OH                       ', &
                                      'N2D_O                           ', 'N2D_O2                          ', &
                                      'N_NO                            ', 'N_NO2a                          ', &
                                      'N_NO2b                          ', 'N_NO2c                          ', &
                                      'N_O2                            ', 'NO2_O                           ', &
                                      'NO2_O3                          ', 'NO2_O_M                         ', &
                                      'NO3_HO2                         ', 'NO3_NO                          ', &
                                      'NO3_O                           ', 'NO3_OH                          ', &
                                      'N_OH                            ', 'NO_HO2                          ', &
                                      'NO_O3                           ', 'NO_O_M                          ', &
                                      'O1D_N2Oa                        ', 'O1D_N2Ob                        ', &
                                      'tag_NO2_HO2                     ', 'tag_NO2_NO3                     ', &
                                      'tag_NO2_OH                      ', 'usr_HNO3_OH                     ', &
                                      'usr_HO2NO2_M                    ', 'usr_N2O5_M                      ', &
                                      'CL_CH2O                         ', 'CL_CH4                          ', &
                                      'CL_H2                           ', 'CL_H2O2                         ', &
                                      'CL_HO2a                         ', 'CL_HO2b                         ', &
                                      'CL_O3                           ', 'CLO_CH3O2                       ', &
                                      'CLO_CLOa                        ', 'CLO_CLOb                        ', &
                                      'CLO_CLOc                        ', 'CLO_HO2                         ', &
                                      'CLO_NO                          ', 'CLONO2_CL                       ', &
                                      'CLO_NO2_M                       ', 'CLONO2_O                        ', &
                                      'CLONO2_OH                       ', 'CLO_O                           ', &
                                      'CLO_OHa                         ', 'CLO_OHb                         ', &
                                      'HCL_O                           ', 'HCL_OH                          ', &
                                      'HOCL_CL                         ', 'HOCL_O                          ', &
                                      'HOCL_OH                         ', 'O1D_CCL4                        ', &
                                      'O1D_CF2CLBR                     ', 'O1D_CFC11                       ', &
                                      'O1D_CFC113                      ', 'O1D_CFC114                      ', &
                                      'O1D_CFC115                      ', 'O1D_CFC12                       ', &
                                      'O1D_HCLa                        ', 'O1D_HCLb                        ', &
                                      'tag_CLO_CLO_M                   ', 'usr_CL2O2_M                     ', &
                                      'BR_CH2O                         ', 'BR_HO2                          ', &
                                      'BR_O3                           ', 'BRO_BRO                         ', &
                                      'BRO_CLOa                        ', 'BRO_CLOb                        ', &
                                      'BRO_CLOc                        ', 'BRO_HO2                         ', &
                                      'BRO_NO                          ', 'BRO_NO2_M                       ', &
                                      'BRONO2_O                        ', 'BRO_O                           ', &
                                      'BRO_OH                          ', 'HBR_O                           ', &
                                      'HBR_OH                          ', 'HOBR_O                          ', &
                                      'O1D_CF3BR                       ', 'O1D_CHBR3                       ', &
                                      'O1D_H2402                       ', 'O1D_HBRa                        ', &
                                      'O1D_HBRb                        ', 'F_CH4                           ', &
                                      'F_H2                            ', 'F_H2O                           ', &
                                      'F_HNO3                          ', 'O1D_COF2                        ', &
                                      'O1D_COFCL                       ', 'CH2BR2_CL                       ', &
                                      'CH2BR2_OH                       ', 'CH3BR_CL                        ', &
                                      'CH3BR_OH                        ', 'CH3CCL3_OH                      ', &
                                      'CH3CL_CL                        ', 'CH3CL_OH                        ', &
                                      'CHBR3_CL                        ', 'CHBR3_OH                        ', &
                                      'HCFC141B_OH                     ', 'HCFC142B_OH                     ', &
                                      'HCFC22_OH                       ', 'O1D_CH2BR2                      ', &
                                      'O1D_CH3BR                       ', 'O1D_HCFC141B                    ', &
                                      'O1D_HCFC142B                    ', 'O1D_HCFC22                      ', &
                                      'CH2O_NO3                        ', 'CH2O_O                          ', &
                                      'CH2O_OH                         ', 'CH3O2_HO2                       ', &
                                      'CH3O2_NO                        ', 'CH3OOH_OH                       ', &
                                      'CH4_OH                          ', 'CO_OH_M                         ', &
                                      'O1D_CH4a                        ', 'O1D_CH4b                        ', &
                                      'O1D_CH4c                        ', 'usr_CO_OH_b                     ', &
                                      'usr_HO2_aer                     ', 'usr_N2O5_aer                    ', &
                                      'usr_NO2_aer                     ', 'usr_NO3_aer                     ', &
                                      'het1                            ', 'het10                           ', &
                                      'het11                           ', 'het12                           ', &
                                      'het13                           ', 'het14                           ', &
                                      'het15                           ', 'het16                           ', &
                                      'het17                           ', 'het2                            ', &
                                      'het3                            ', 'het4                            ', &
                                      'het5                            ', 'het6                            ', &
                                      'het7                            ', 'het8                            ', &
                                      'het9                            ', 'ag247nm                         ', &
                                      'ag373nm                         ', 'ag732nm                         ', &
                                      'elec1                           ', 'elec2                           ', &
                                      'elec3                           ', 'ion_N2p_O2                      ', &
                                      'ion_N2p_Oa                      ', 'ion_N2p_Ob                      ', &
                                      'ion_Np_O                        ', 'ion_Np_O2a                      ', &
                                      'ion_Np_O2b                      ', 'ion_O2p_N                       ', &
                                      'ion_O2p_N2                      ', 'ion_O2p_NO                      ', &
                                      'ion_Op_CO2                      ', 'ion_Op_N2                       ', &
                                      'ion_Op_N2D                      ', 'ion_Op_O2                       ', &
                                      'Op2D_e                          ', 'Op2D_N2                         ', &
                                      'Op2D_O                          ', 'Op2D_O2                         ', &
                                      'Op2P_ea                         ', 'Op2P_eb                         ', &
                                      'Op2P_N2a                        ', 'Op2P_N2b                        ', &
                                      'Op2P_O                          ' /)
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
                                      311, 312, 313, 314, 315, 316, 317, 318, 319, 320, &
                                      321, 322, 323, 324, 325, 326, 327, 328, 329, 330, &
                                      331, 332, 333, 334, 335, 336, 337, 338, 339, 340, &
                                      341, 342, 343, 344, 345, 346, 347, 348, 349, 350, &
                                      351, 352, 353, 354, 355, 356, 357, 358, 359, 360, &
                                      361, 362, 363, 364, 365, 366, 367, 368, 369, 370, &
                                      371, 372, 373, 374, 375, 376, 377, 378, 379, 380, &
                                      381, 382, 383, 384, 385, 386, 387, 388, 389, 390, &
                                      391, 392, 393, 394, 395, 396, 397, 398, 399, 400, &
                                      401, 402, 403, 404, 405, 406, 407, 408, 409, 410, &
                                      411, 412, 413, 414, 415, 416, 417, 418, 419, 420, &
                                      421, 422, 423, 424, 425, 426, 427, 428, 429, 430, &
                                      431, 432, 433, 434, 435, 436, 437, 438, 439, 440, &
                                      441, 442, 443, 444, 445, 446, 447, 448, 449, 450, &
                                      451, 452, 453, 454, 455, 456, 457, 458, 459, 460, &
                                      461, 462, 463, 464, 465, 466, 467, 468, 469, 470, &
                                      471, 472, 473, 474, 475, 476, 477, 478, 479, 480, &
                                      481, 482, 483, 484, 485, 486, 487, 488, 489, 490, &
                                      491, 492, 493, 494, 495, 496, 497, 498, 499, 500, &
                                      501, 502, 503, 504, 505, 506, 507, 508, 509, 510, &
                                      511, 512, 513, 514, 515, 516, 517, 518, 519, 520, &
                                      521, 522, 523, 524, 525, 526, 527, 528, 529, 530, &
                                      531, 532, 533, 534, 535, 536, 537, 538, 539, 540, &
                                      541, 542, 543, 544, 545, 546, 547, 548, 549, 550, &
                                      551, 552, 553, 554, 555, 556, 557, 558, 559, 560, &
                                      561, 562, 563, 564, 565, 566, 567, 568, 569, 570, &
                                      571, 572, 573, 574, 575, 576, 577, 578, 579, 580, &
                                      581, 582, 583, 584, 585 /)
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
      pht_alias_lst(:,1) = (/ '                ', '                ', '                ', '                ', &
                              'userdefined     ', 'userdefined     ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              'userdefined     ', '                ', '                ', '                ', &
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
                              '                ', 'userdefined     ', 'userdefined     ', 'userdefined     ', &
                              'userdefined     ', 'userdefined     ', 'userdefined     ', 'userdefined     ', &
                              'userdefined     ', 'userdefined     ', 'userdefined     ', 'userdefined     ', &
                              'userdefined     ', 'userdefined     ', 'userdefined     ', 'userdefined     ', &
                              'userdefined     ', 'userdefined     ', 'userdefined     ', 'userdefined     ', &
                              'userdefined     ', 'userdefined     ', 'userdefined     ', 'userdefined     ', &
                              'userdefined     ', 'userdefined     ', 'userdefined     ', 'userdefined     ', &
                              'userdefined     ', 'userdefined     ', 'userdefined     ', 'userdefined     ', &
                              'userdefined     ', 'userdefined     ', 'userdefined     ', 'userdefined     ', &
                              'userdefined     ', 'userdefined     ', 'userdefined     ' /)
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
                              '                ', 'userdefined     ', '                ', '                ', &
                              '                ', 'userdefined     ', 'userdefined     ', 'userdefined     ', &
                              'userdefined     ', 'userdefined     ', 'userdefined     ', 'userdefined     ', &
                              'userdefined     ', 'userdefined     ', '                ', '                ', &
                              'userdefined     ', 'userdefined     ', 'userdefined     ', 'userdefined     ', &
                              'userdefined     ', 'userdefined     ', 'userdefined     ', 'userdefined     ', &
                              'userdefined     ', 'userdefined     ', 'userdefined     ', 'userdefined     ', &
                              'userdefined     ', 'userdefined     ', 'userdefined     ', 'userdefined     ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ' /)
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
      cph_rid(:)      = (/             386,            387,            388,            390,            391, &
                                       392,            394,            395,            396,            397, &
                                       398,            399,            400,            403,            406, &
                                       407,            408,            409,            412,            413, &
                                       414,            417,            419,            420,            421, &
                                       425,            426,            434,            435,            558, &
                                       559,            560,            561,            562,            563, &
                                       564,            565,            567,            568,            569, &
                                       570,            572,            574,            575,            576, &
                                       577,            578,            579,            580,            581, &
                                       582,            583,            584,            585 /)
      cph_enthalpy(:) = (/   189.810000_r8,   32.910000_r8,  189.810000_r8,   94.300000_r8,   94.300000_r8, &
                              94.300000_r8,   62.600000_r8,   62.600000_r8,   62.600000_r8,   62.600000_r8, &
                             392.190000_r8,  493.580000_r8,  101.390000_r8,  232.590000_r8,  203.400000_r8, &
                             226.580000_r8,  120.100000_r8,  194.710000_r8,  293.620000_r8,   67.670000_r8, &
                             165.300000_r8,  165.510000_r8,  229.610000_r8,  177.510000_r8,  313.750000_r8, &
                             133.750000_r8,  193.020000_r8,   34.470000_r8,  199.170000_r8,  483.390000_r8, &
                             321.300000_r8,  163.060000_r8,   82.389000_r8,  508.950000_r8,  354.830000_r8, &
                             339.590000_r8,   67.530000_r8,   95.550000_r8,  239.840000_r8,  646.280000_r8, &
                             406.160000_r8,  271.380000_r8,  105.040000_r8,  139.900000_r8,  150.110000_r8, &
                             319.370000_r8,  128.320000_r8,  319.360000_r8,  469.400000_r8,  163.060000_r8, &
                             482.430000_r8,  291.380000_r8,   67.540000_r8,  501.720000_r8 /)
      allocate( num_rnts(rxntot-phtcnt),stat=ios )
      if( ios /= 0 ) then
         write(iulog,*) 'set_sim_dat: failed to allocate num_rnts; error = ',ios
         call endrun
      end if
      num_rnts(:) = (/      2,     2,     3,     3,     2,     2,     2,     2,     2,     2, &
                            2,     2,     2,     3,     3,     2,     2,     2,     2,     2, &
                            2,     2,     2,     2,     2,     2,     2,     2,     3,     2, &
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
                            2,     2,     2,     2,     2,     2,     2,     2,     2,     2, &
                            2,     2,     2,     3,     3,     3,     3,     3,     3,     3, &
                            3,     3,     3,     3,     3,     3,     3,     2,     2,     2, &
                            2,     3,     2,     2,     2,     2,     2,     2,     2,     2, &
                            2,     2,     2,     3,     2,     2,     2,     2,     2,     3, &
                            2,     3,     2,     2,     2,     3,     2,     2,     2,     2, &
                            2,     2,     2,     2,     2,     2,     2,     2,     2,     2, &
                            2,     3,     2,     2,     2,     2,     2,     3,     2,     2, &
                            2,     2,     2,     2,     3,     2,     2,     2,     2,     2, &
                            2,     2,     2,     2,     2,     2,     2,     2,     3,     2, &
                            2,     2,     2,     2,     3,     2,     2,     3,     2,     2, &
                            2,     2,     2,     2,     2,     2,     2,     3,     2,     3, &
                            3,     3,     3,     3,     2,     3,     2,     3,     2,     3, &
                            2,     3,     2,     3,     2,     2,     2,     2,     2,     2, &
                            2,     2,     1,     1,     2,     2,     2,     2,     2,     2, &
                            2,     2,     2,     2,     2,     2,     2,     2,     2,     3, &
                            3,     2,     2,     2,     2,     2,     3,     2,     2,     2, &
                            2,     2,     2,     2,     2,     2,     3,     2,     2,     2, &
                            2,     2,     2,     2,     2,     2,     2,     2,     3,     2, &
                            2,     2,     2,     2,     2,     2,     3,     2,     2,     3, &
                            3,     3,     2,     2,     2,     2,     2,     2,     2,     2, &
                            2,     2,     2,     2,     2,     2,     2,     2,     2,     3, &
                            2,     2,     2,     2,     2,     2,     2,     2,     2,     2, &
                            2,     2,     2,     2,     2,     2,     2,     2,     2,     3, &
                            2,     2,     2,     2,     2,     2,     2,     2,     2,     2, &
                            3,     2,     2,     2,     2,     2,     2,     2,     2,     2, &
                            2,     2,     2,     2,     2,     2,     2,     2,     2,     2, &
                            2,     2,     2,     2,     2,     2,     2,     2,     2,     2, &
                            2,     2,     2,     2,     2,     2,     2,     2,     2,     2, &
                            2,     2,     3,     2,     2,     2,     2,     1,     1,     1, &
                            1,     1,     2,     1,     1,     1,     1,     2,     2,     2, &
                            1,     1,     2,     2,     2,     1,     1,     2,     1,     1, &
                            1,     2,     2,     2,     2,     2,     2,     2,     2,     2, &
                            2,     2,     2,     2,     2,     2,     2,     2,     2,     2, &
                            2,     2,     2,     2,     2,     2 /)

      end subroutine set_sim_dat

      end module mo_sim_dat
