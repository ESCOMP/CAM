
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

      clscnt(:) = (/     23,     0,     0,    50,     0 /)

      cls_rxt_cnt(:,1) = (/      3,    59,     0,    23 /)
      cls_rxt_cnt(:,4) = (/     30,   121,   138,    50 /)

      solsym(: 73) = (/ 'BRCL            ','BRO             ','BRONO2          ','BRY             ','CCL4            ', &
                        'CF2CLBR         ','CF3BR           ','CFC11           ','CFC113          ','CFC114          ', &
                        'CFC115          ','CFC12           ','CH2BR2          ','CH2O            ','CH3BR           ', &
                        'CH3CCL3         ','CH3CL           ','CH3O2           ','CH3OOH          ','CH4             ', &
                        'CHBR3           ','CL2             ','CL2O2           ','CLO             ','CLONO2          ', &
                        'CLY             ','CO              ','CO2             ','COF2            ','COFCL           ', &
                        'F               ','H               ','H2              ','H2402           ','H2O2            ', &
                        'HBR             ','HCFC141B        ','HCFC142B        ','HCFC22          ','HCL             ', &
                        'HF              ','HNO3            ','HO2             ','HO2NO2          ','HOBR            ', &
                        'HOCL            ','N               ','N2O             ','N2O5            ','NO              ', &
                        'NO2             ','NO3             ','O               ','O2              ','O3              ', &
                        'OCLO            ','SF6             ','BR              ','CL              ','e               ', &
                        'N2D             ','N2p             ','NOp             ','Np              ','O1D             ', &
                        'O2_1D           ','O2_1S           ','O2p             ','OH              ','Op              ', &
                        'Op2D            ','Op2P            ','H2O             ' /)

      adv_mass(: 73) = (/   115.356700_r8,    95.903400_r8,   141.908940_r8,    99.716850_r8,   153.821800_r8, &
                            165.364506_r8,   148.910210_r8,   137.367503_r8,   187.375310_r8,   170.921013_r8, &
                            154.466716_r8,   120.913206_r8,   173.833800_r8,    30.025200_r8,    94.937200_r8, &
                            133.402300_r8,    50.485900_r8,    47.032000_r8,    48.039400_r8,    16.040600_r8, &
                            252.730400_r8,    70.905400_r8,   102.904200_r8,    51.452100_r8,    97.457640_r8, &
                            100.916850_r8,    28.010400_r8,    44.009800_r8,    66.007206_r8,    82.461503_r8, &
                             18.998403_r8,     1.007400_r8,     2.014800_r8,   259.823613_r8,    34.013600_r8, &
                             80.911400_r8,   116.948003_r8,   100.493706_r8,    86.467906_r8,    36.460100_r8, &
                             20.005803_r8,    63.012340_r8,    33.006200_r8,    79.011740_r8,    96.910800_r8, &
                             52.459500_r8,    14.006740_r8,    44.012880_r8,   108.010480_r8,    30.006140_r8, &
                             46.005540_r8,    62.004940_r8,    15.999400_r8,    31.998800_r8,    47.998200_r8, &
                             67.451500_r8,   146.056419_r8,    79.904000_r8,    35.452700_r8, 0.548567E-03_r8, &
                             14.006740_r8,    28.013480_r8,    30.006140_r8,    14.006740_r8,    15.999400_r8, &
                             31.998800_r8,    31.998800_r8,    31.998800_r8,    17.006800_r8,    15.999400_r8, &
                             15.999400_r8,    15.999400_r8,    18.014200_r8 /)

      crb_mass(: 73) = (/     0.000000_r8,     0.000000_r8,     0.000000_r8,     0.000000_r8,    12.011000_r8, &
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
                              0.000000_r8,     0.000000_r8,     0.000000_r8,     0.000000_r8,     0.000000_r8, &
                              0.000000_r8,     0.000000_r8,     0.000000_r8,     0.000000_r8,     0.000000_r8, &
                              0.000000_r8,     0.000000_r8,     0.000000_r8 /)

      fix_mass(:  2) = (/ 0.00000000_r8, 28.0134800_r8 /)

      clsmap(: 23,1) = (/    4,   5,   6,   7,   8,   9,  10,  11,  12,  13, &
                            15,  16,  17,  20,  21,  26,  28,  34,  37,  38, &
                            39,  48,  57 /)
      clsmap(: 50,4) = (/    1,   2,   3,  14,  18,  19,  22,  23,  24,  25, &
                            27,  29,  30,  31,  32,  33,  35,  36,  40,  41, &
                            42,  43,  44,  45,  46,  47,  49,  50,  51,  52, &
                            53,  54,  55,  56,  58,  59,  60,  61,  62,  63, &
                            64,  65,  66,  67,  68,  69,  70,  71,  72,  73 /)

      permute(: 50,4) = (/    8,  39,  15,  35,  44,  11,   6,   1,  41,  30, &
                             16,   2,   7,  21,  36,  49,  17,  23,  34,   9, &
                             31,  47,  12,  20,  22,  29,  10,  38,  45,  37, &
                             48,  32,  50,   3,  43,  40,  27,  28,  18,  19, &
                             24,  42,   4,   5,  26,  46,  25,  14,  13,  33 /)

      diag_map(: 50) = (/    1,   4,   7,  10,  13,  15,  17,  21,  24,  27, &
                            33,  39,  46,  53,  59,  67,  71,  78,  87,  93, &
                           102, 111, 118, 127, 138, 149, 164, 178, 193, 204, &
                           216, 236, 248, 264, 282, 298, 316, 339, 363, 389, &
                           416, 442, 463, 483, 509, 541, 567, 607, 628, 647 /)

      extfrc_lst(: 11) = (/ 'NO2             ','CO              ','NO              ','O2p             ','Op              ', &
                            'N2p             ','N               ','OH              ','e               ','Np              ', &
                            'N2D             ' /)

      frc_from_dataset(: 11) = (/ .true., .true., .true., .false., .false., &
                                  .false., .false., .false., .false., .false., &
                                  .false. /)

      inv_lst(:  2) = (/ 'M               ', 'N2              ' /)

      slvd_lst(: 15) = (/ 'BR              ', 'CL              ', 'e               ', 'N2D             ', 'N2p             ', &
                          'NOp             ', 'Np              ', 'O1D             ', 'O2_1D           ', 'O2_1S           ', &
                          'O2p             ', 'OH              ', 'Op              ', 'Op2D            ', 'Op2P            ' /)

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
                                      'jho2no2_b                       ', 'jn2o                            ', &
                                      'jn2o5_a                         ', 'jn2o5_b                         ', &
                                      'jno                             ', 'jno_i                           ', &
                                      'jno2                            ', 'jno3_a                          ', &
                                      'jno3_b                          ', 'jch2o_a                         ', &
                                      'jch2o_b                         ', 'jch3ooh                         ', &
                                      'jch4_a                          ', 'jch4_b                          ', &
                                      'jco2                            ', 'jbrcl                           ', &
                                      'jbro                            ', 'jbrono2_b                       ', &
                                      'jbrono2_a                       ', 'jccl4                           ', &
                                      'jcf2clbr                        ', 'jcf3br                          ', &
                                      'jcfcl3                          ', 'jcfc113                         ', &
                                      'jcfc114                         ', 'jcfc115                         ', &
                                      'jcf2cl2                         ', 'jch2br2                         ', &
                                      'jch3br                          ', 'jch3ccl3                        ', &
                                      'jch3cl                          ', 'jchbr3                          ', &
                                      'jcl2                            ', 'jcl2o2                          ', &
                                      'jclo                            ', 'jclono2_a                       ', &
                                      'jclono2_b                       ', 'jcof2                           ', &
                                      'jcofcl                          ', 'jh2402                          ', &
                                      'jhbr                            ', 'jhcfc141b                       ', &
                                      'jhcfc142b                       ', 'jhcfc22                         ', &
                                      'jhcl                            ', 'jhf                             ', &
                                      'jhobr                           ', 'jhocl                           ', &
                                      'joclo                           ', 'jsf6                            ', &
                                      'jeuv_26                         ', 'jeuv_4                          ', &
                                      'jeuv_23                         ', 'jeuv_25                         ', &
                                      'jeuv_6                          ', 'jeuv_22                         ', &
                                      'jeuv_18                         ', 'jeuv_13                         ', &
                                      'jeuv_11                         ', 'jeuv_10                         ', &
                                      'jeuv_1                          ', 'jeuv_14                         ', &
                                      'jeuv_2                          ', 'jeuv_3                          ', &
                                      'jeuv_15                         ', 'jeuv_16                         ', &
                                      'jeuv_9                          ', 'jeuv_17                         ', &
                                      'jeuv_7                          ', 'jeuv_5                          ', &
                                      'jeuv_19                         ', 'jeuv_24                         ', &
                                      'jeuv_12                         ', 'jeuv_20                         ', &
                                      'jeuv_21                         ', 'jeuv_8                          ', &
                                      'ag1                             ', 'ag2                             ', &
                                      'O1D_H2                          ', 'O1D_H2O                         ', &
                                      'O1D_N2                          ', 'O1D_O2                          ', &
                                      'O1D_O2b                         ', 'O1D_O3                          ', &
                                      'O2_1D_N2                        ', 'O2_1D_O                         ', &
                                      'O2_1D_O2                        ', 'O2_1S_CO2                       ', &
                                      'O2_1S_N2                        ', 'O2_1S_O                         ', &
                                      'O2_1S_O2                        ', 'O2_1S_O3                        ', &
                                      'O_O3                            ', 'usr_O_O                         ', &
                                      'usr_O_O2                        ', 'H2_O                            ', &
                                      'H2O2_O                          ', 'H_HO2                           ', &
                                      'H_HO2a                          ', 'H_HO2b                          ', &
                                      'H_O2                            ', 'HO2_O                           ', &
                                      'HO2_O3                          ', 'H_O3                            ', &
                                      'OH_H2                           ', 'OH_H2O2                         ', &
                                      'OH_HO2                          ', 'OH_O                            ', &
                                      'OH_O3                           ', 'OH_OH                           ', &
                                      'OH_OH_M                         ', 'usr_HO2_HO2                     ', &
                                      'HO2NO2_OH                       ', 'N2D_O                           ', &
                                      'N2D_O2                          ', 'N_NO                            ', &
                                      'N_NO2a                          ', 'N_NO2b                          ', &
                                      'N_NO2c                          ', 'N_O2                            ', &
                                      'NO2_O                           ', 'NO2_O3                          ', &
                                      'NO2_O_M                         ', 'NO3_HO2                         ', &
                                      'NO3_NO                          ', 'NO3_O                           ', &
                                      'NO3_OH                          ', 'N_OH                            ', &
                                      'NO_HO2                          ', 'NO_O3                           ', &
                                      'NO_O_M                          ', 'O1D_N2Oa                        ', &
                                      'O1D_N2Ob                        ', 'tag_NO2_HO2                     ', &
                                      'tag_NO2_NO3                     ', 'tag_NO2_OH                      ', &
                                      'usr_HNO3_OH                     ', 'usr_HO2NO2_M                    ', &
                                      'usr_N2O5_M                      ', 'CL_CH2O                         ', &
                                      'CL_CH4                          ', 'CL_H2                           ', &
                                      'CL_H2O2                         ', 'CL_HO2a                         ', &
                                      'CL_HO2b                         ', 'CL_O3                           ', &
                                      'CLO_CH3O2                       ', 'CLO_CLOa                        ', &
                                      'CLO_CLOb                        ', 'CLO_CLOc                        ', &
                                      'CLO_HO2                         ', 'CLO_NO                          ', &
                                      'CLONO2_CL                       ', 'CLO_NO2_M                       ', &
                                      'CLONO2_O                        ', 'CLONO2_OH                       ', &
                                      'CLO_O                           ', 'CLO_OHa                         ', &
                                      'CLO_OHb                         ', 'HCL_O                           ', &
                                      'HCL_OH                          ', 'HOCL_CL                         ', &
                                      'HOCL_O                          ', 'HOCL_OH                         ', &
                                      'O1D_CCL4                        ', 'O1D_CF2CLBR                     ', &
                                      'O1D_CFC11                       ', 'O1D_CFC113                      ', &
                                      'O1D_CFC114                      ', 'O1D_CFC115                      ', &
                                      'O1D_CFC12                       ', 'O1D_HCLa                        ', &
                                      'O1D_HCLb                        ', 'tag_CLO_CLO_M                   ', &
                                      'usr_CL2O2_M                     ', 'BR_CH2O                         ', &
                                      'BR_HO2                          ', 'BR_O3                           ', &
                                      'BRO_BRO                         ', 'BRO_CLOa                        ', &
                                      'BRO_CLOb                        ', 'BRO_CLOc                        ', &
                                      'BRO_HO2                         ', 'BRO_NO                          ', &
                                      'BRO_NO2_M                       ', 'BRONO2_O                        ', &
                                      'BRO_O                           ', 'BRO_OH                          ', &
                                      'HBR_O                           ', 'HBR_OH                          ' /)
      rxt_tag_lst(   201:   290) = (/ 'HOBR_O                          ', 'O1D_CF3BR                       ', &
                                      'O1D_CHBR3                       ', 'O1D_H2402                       ', &
                                      'O1D_HBRa                        ', 'O1D_HBRb                        ', &
                                      'F_CH4                           ', 'F_H2                            ', &
                                      'F_H2O                           ', 'F_HNO3                          ', &
                                      'O1D_COF2                        ', 'O1D_COFCL                       ', &
                                      'CH2BR2_CL                       ', 'CH2BR2_OH                       ', &
                                      'CH3BR_CL                        ', 'CH3BR_OH                        ', &
                                      'CH3CCL3_OH                      ', 'CH3CL_CL                        ', &
                                      'CH3CL_OH                        ', 'CHBR3_CL                        ', &
                                      'CHBR3_OH                        ', 'HCFC141B_OH                     ', &
                                      'HCFC142B_OH                     ', 'HCFC22_OH                       ', &
                                      'O1D_CH2BR2                      ', 'O1D_CH3BR                       ', &
                                      'O1D_HCFC141B                    ', 'O1D_HCFC142B                    ', &
                                      'O1D_HCFC22                      ', 'CH2O_NO3                        ', &
                                      'CH2O_O                          ', 'CH2O_OH                         ', &
                                      'CH3O2_HO2                       ', 'CH3O2_NO                        ', &
                                      'CH3OOH_OH                       ', 'CH4_OH                          ', &
                                      'CO_OH_M                         ', 'O1D_CH4a                        ', &
                                      'O1D_CH4b                        ', 'O1D_CH4c                        ', &
                                      'usr_CO_OH_b                     ', 'usr_HO2_aer                     ', &
                                      'usr_N2O5_aer                    ', 'usr_NO2_aer                     ', &
                                      'usr_NO3_aer                     ', 'het1                            ', &
                                      'het10                           ', 'het11                           ', &
                                      'het12                           ', 'het13                           ', &
                                      'het14                           ', 'het15                           ', &
                                      'het16                           ', 'het17                           ', &
                                      'het2                            ', 'het3                            ', &
                                      'het4                            ', 'het5                            ', &
                                      'het6                            ', 'het7                            ', &
                                      'het8                            ', 'het9                            ', &
                                      'ag247nm                         ', 'ag373nm                         ', &
                                      'ag732nm                         ', 'elec1                           ', &
                                      'elec2                           ', 'elec3                           ', &
                                      'ion_N2p_O2                      ', 'ion_N2p_Oa                      ', &
                                      'ion_N2p_Ob                      ', 'ion_Np_O                        ', &
                                      'ion_Np_O2a                      ', 'ion_Np_O2b                      ', &
                                      'ion_O2p_N                       ', 'ion_O2p_N2                      ', &
                                      'ion_O2p_NO                      ', 'ion_Op_CO2                      ', &
                                      'ion_Op_N2                       ', 'ion_Op_N2D                      ', &
                                      'ion_Op_O2                       ', 'Op2D_e                          ', &
                                      'Op2D_N2                         ', 'Op2D_O                          ', &
                                      'Op2D_O2                         ', 'Op2P_ea                         ', &
                                      'Op2P_eb                         ', 'Op2P_N2a                        ', &
                                      'Op2P_N2b                        ', 'Op2P_O                          ' /)
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
                                      281, 282, 283, 284, 285, 286, 287, 288, 289, 290 /)
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
                              '                ', '                ', 'userdefined     ', '                ', &
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
                              'userdefined     ', 'userdefined     ', 'userdefined     ', 'userdefined     ', &
                              'userdefined     ', 'userdefined     ', 'userdefined     ', 'userdefined     ', &
                              'userdefined     ', 'userdefined     ', 'userdefined     ', 'userdefined     ', &
                              'userdefined     ', 'userdefined     ', 'userdefined     ', 'userdefined     ', &
                              'userdefined     ', 'userdefined     ', 'userdefined     ', 'userdefined     ', &
                              'userdefined     ', 'userdefined     ', 'userdefined     ', 'userdefined     ', &
                              'userdefined     ', 'userdefined     ' /)
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
                              'userdefined     ', 'userdefined     ', 'userdefined     ', 'userdefined     ', &
                              'userdefined     ', 'userdefined     ', 'userdefined     ', 'userdefined     ', &
                              'userdefined     ', 'userdefined     ', 'userdefined     ', 'userdefined     ', &
                              'userdefined     ', 'userdefined     ', 'userdefined     ', 'userdefined     ', &
                              'userdefined     ', 'userdefined     ', 'userdefined     ', 'userdefined     ', &
                              'userdefined     ', 'userdefined     ', 'userdefined     ', 'userdefined     ', &
                              'userdefined     ', 'userdefined     ' /)
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
                          1._r8 /)
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
                          1._r8 /)
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
      cph_rid(:)      = (/              91,             92,             93,             95,             96, &
                                        97,             99,            100,            101,            102, &
                                       103,            104,            105,            108,            111, &
                                       112,            113,            114,            117,            118, &
                                       119,            122,            124,            125,            126, &
                                       130,            131,            139,            140,            263, &
                                       264,            265,            266,            267,            268, &
                                       269,            270,            272,            273,            274, &
                                       275,            277,            279,            280,            281, &
                                       282,            283,            284,            285,            286, &
                                       287,            288,            289,            290 /)
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
      num_rnts(:) = (/      1,     1,     2,     2,     2,     2,     2,     2,     2,     2, &
                            2,     2,     2,     2,     2,     2,     2,     3,     3,     2, &
                            2,     2,     2,     2,     3,     2,     2,     2,     2,     2, &
                            2,     2,     2,     2,     3,     2,     2,     2,     2,     2, &
                            2,     2,     2,     2,     2,     2,     3,     2,     2,     2, &
                            2,     2,     2,     2,     3,     2,     2,     3,     3,     3, &
                            2,     2,     2,     2,     2,     2,     2,     2,     2,     2, &
                            2,     2,     2,     2,     2,     2,     2,     3,     2,     2, &
                            2,     2,     2,     2,     2,     2,     2,     2,     2,     2, &
                            2,     2,     2,     2,     2,     2,     2,     3,     2,     2, &
                            2,     2,     2,     2,     2,     2,     2,     2,     3,     2, &
                            2,     2,     2,     2,     2,     2,     2,     2,     2,     2, &
                            2,     2,     2,     2,     2,     2,     2,     2,     2,     2, &
                            2,     2,     2,     2,     2,     2,     2,     2,     2,     2, &
                            2,     2,     2,     2,     2,     2,     2,     2,     2,     2, &
                            3,     2,     2,     2,     2,     1,     1,     1,     1,     1, &
                            2,     1,     1,     1,     1,     2,     2,     2,     1,     1, &
                            2,     2,     2,     1,     1,     2,     1,     1,     1,     2, &
                            2,     2,     2,     2,     2,     2,     2,     2,     2,     2, &
                            2,     2,     2,     2,     2,     2,     2,     2,     2,     2, &
                            2,     2,     2,     2 /)

      end subroutine set_sim_dat

      end module mo_sim_dat
