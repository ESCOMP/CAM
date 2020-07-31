!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: getLandTypes.F90
!
! !DESCRIPTION: Subroutine getLandTypes converts the land types and leaf
!  area indices from the land model to the LandTypeFrac and XLAI_NATIVE
!  arrays in GEOS-Chem.
!
! !INTERFACE:
!
    SUBROUTINE getLandTypes( cam_in, nY, State_Met )
!
! !USES:
!
    USE camsrfexch,       ONLY : cam_in_t
    USE State_Met_Mod,    ONLY : MetState
    USE seq_drydep_mod,   ONLY : NPatch
    USE shr_kind_mod,     ONLY : r8 => shr_kind_r8
    USE PRECISION_MOD,    ONLY : fp, f4     ! Flexible precision
    USE CMN_SIZE_Mod,     ONLY : NSURFTYPE
    USE cam_abortutils,   ONLY : endrun

    IMPLICIT NONE
!
! !INPUT PARAMETERS:
!
    TYPE(cam_in_t),  INTENT(IN   ) :: cam_in    ! CAM
    INTEGER,         INTENT(IN   ) :: nY        ! Number of grid cells on chunk
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(MetState),  INTENT(INOUT) :: State_Met
!
! !REVISION HISTORY:
!  8 May 2020 - Thibaud M. Fritz - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER  :: J, T
    REAL(r8) :: waterFrac, landFrac

#if defined( CLM40 )

    ! Mapping for CLM4.0
    ! -----------------------------------|--------------------------------------
    !           Olson land type          |             CLM land type
    ! -----------------------------------|--------------------------------------
    ! Inland/sea water         (ID =  1) | Ocean fraction
    !                                    | Deeplake & Shallowlake    (LUID =3/4)
    ! Urban                    (ID =  2) | Urban - Not Applied       (LUID =  6)
    ! Low Sparse Grassland     (ID =  3) | 
    ! Coniferous Forest        (ID =  4) | 
    ! Deciduous Conifer Forest (ID =  5) | Needleleaf Deciduous Bor. (PAID =  3)
    ! Deciduous Broadleaf For. (ID =  6) |
    ! Evergreen Broadleaf For. (ID =  7) | Broadleaf Evergreen Temp. (PAID =  5)
    ! Tall Grasses and Shrubs  (ID =  8) | 
    ! Bare Desert              (ID =  9) | Not veg. \ Ice    (PAID = 0\LUID = 2)
    ! Upland Tundra            (ID = 10) | Broadleaf Deciduous Bore. (PAID = 11)
    ! Irrigated Grassland      (ID = 11) |
    ! Semi Desert              (ID = 12) | 
    ! Glacier ice              (ID = 13) | Land ice                  (LUID =  2)
    ! Wooded Wet Swamp         (ID = 14) | 
    ! -                        (ID = 15) | 
    ! -                        (ID = 16) | 
    ! Shrub Evergreen          (ID = 17) | Broadleaf Evergreen Shru. (PAID =  9)
    ! -                        (ID = 18) | 
    ! Shrub Deciduous          (ID = 19) | 
    ! Evergreen Forest and Fi. (ID = 20) | 
    ! Cool Rain Forest         (ID = 21) | 
    ! Conifer Boreal Forest    (ID = 22) | Needleleaf Evergreen Bor. (PAID =  2) 
    ! Cool Conifer Forest      (ID = 23) | 
    ! Cool Mixed Forest        (ID = 24) | Broadleaf Deciduous Bore. (PAID =  8)
    ! Mixed Forest             (ID = 25) | 
    ! Cool Broadleaf Forest    (ID = 26) | Broadleaf Deciduous Temp. (PAID =  7)
    ! Deciduous Broadleaf For. (ID = 27) | 
    ! Conifer Forest           (ID = 28) | Needleleaf Evergreen Tem. (PAID =  1)
    ! Montane Tropical Forests (ID = 29) | 
    ! Seasonal Tropical Fores. (ID = 30) | 
    ! Cool Crops and Towns     (ID = 31) | Winter Temp. Cereal       (PAID = 19) 
    ! Crops and Town           (ID = 32) | C3 Crop                   (PAID = 15)
    !                                    | C3 Irrigated              (PAID = 16)
    !                                    | Spring Temp. Cereal       (PAID = 18)
    ! Dry Tropical Woods       (ID = 33) | 
    ! Tropical Rainforest      (ID = 34) | Broadleaf Evergreen Trop. (PAID =  4)
    ! Tropical Degraded Forest (ID = 35) | 
    ! Corn and Beans Cropland  (ID = 36) | Corn                      (PAID = 17)
    !                                    | Soybean                   (PAID = 20)
    ! Rice Paddy and Field     (ID = 37) | 
    ! Hot Irrigated Cropland   (ID = 38) | 
    ! Cool Irrigated Cropland  (ID = 39) | 
    ! -                        (ID = 40) | 
    ! Cool Grasses and Shrubs  (ID = 41) | 
    ! Hot and Mild Grasses and (ID = 42) | C3 Non-Arctic Grass       (PAID = 13)
    ! Cold Grassland           (ID = 43) | C3 Arctic Grass           (PAID = 12)
    ! Savanna (Woods)          (ID = 44) | Broadleaf Deciduous Trop. (PAID =  6)
    !                                    | C4 Grass                  (PAID = 14)
    ! Mire, Bog, Fen           (ID = 45) | Wetland - Not Applied     (LUID =  5)
    ! Marsh Wetland            (ID = 46) |
    ! Mediterranean Scrub      (ID = 47) |
    ! Dry Woody Scrub          (ID = 48) | Broadleaf Deciduous Temp. (PAID = 10)
    ! -                        (ID = 49) | 
    ! -                        (ID = 50) | 
    ! -                        (ID = 51) | 
    ! Semi Desert Shrubs       (ID = 52) |
    ! Semi Desert Sage         (ID = 53) | 
    ! Barren Tundra            (ID = 54) | 
    ! Cool Southern Hemisphere (ID = 55) |
    ! Cool Fields and Woods    (ID = 56) | 
    ! Forest and Field         (ID = 57) | 
    ! Cool Forest and Field    (ID = 58) | 
    ! Fields and Woody Savanna (ID = 59) | 
    ! Succulent and Thorn Scr. (ID = 60) | 
    ! Small Leaf Mixed Woods   (ID = 61) | 
    ! Deciduous and Mixed Bor. (ID = 62) | 
    ! Narrow Conifers          (ID = 63) | 
    ! Wooded Tundra            (ID = 64) | 
    ! Heath Scrub              (ID = 65) | 
    ! -                        (ID = 66) | 
    ! -                        (ID = 67) | 
    ! -                        (ID = 68) | 
    ! -                        (ID = 69) | 
    ! Polar and Alpine Desert  (ID = 70) | 
    ! -                        (ID = 71) | 
    ! -                        (ID = 72) | 
    ! Mangrove                 (ID = 73) | 

    !==================================================================
    ! The urban and wetland land unit types seem to be already
    ! accounted for in patches, as it introduces total land fractions
    ! (summed over all types) greater than 100%.
    ! Thibaud M. Fritz - 06 May 2020
    !==================================================================

    DO J = 1, nY
        waterFrac = cam_in%ocnFrac(J)    + cam_in%iceFrac(J)    &
                  + cam_in%lwtgcell(J,3) + cam_in%lwtgcell(J,4)
        landFrac  = 1.0e+0_fp - waterFrac

        ! Initialize fraction land for this grid cell
        State_Met%LandTypeFrac(1,J, 1) = waterFrac
        !State_Met%LandTypeFrac(1,J, 2) = cam_in%lwtgcell(J, 6)
        State_Met%LandTypeFrac(1,J, 5) = cam_in%pwtgcell(J, 4)
        State_Met%LandTypeFrac(1,J, 7) = cam_in%pwtgcell(J, 6)
        State_Met%LandTypeFrac(1,J, 9) = cam_in%pwtgcell(J, 1) &
                                       - cam_in%lwtgcell(J, 2)
        State_Met%LandTypeFrac(1,J,10) = cam_in%pwtgcell(J,12)
        State_Met%LandTypeFrac(1,J,13) = cam_in%lwtgcell(J, 2)
        State_Met%LandTypeFrac(1,J,17) = cam_in%pwtgcell(J,10)
        State_Met%LandTypeFrac(1,J,22) = cam_in%pwtgcell(J, 3)
        State_Met%LandTypeFrac(1,J,24) = cam_in%pwtgcell(J, 9)
        State_Met%LandTypeFrac(1,J,26) = cam_in%pwtgcell(J, 8)
        State_Met%LandTypeFrac(1,J,28) = cam_in%pwtgcell(J, 2)
        State_Met%LandTypeFrac(1,J,31) = cam_in%pwtgcell(J,20)
        State_Met%LandTypeFrac(1,J,32) = cam_in%pwtgcell(J,16) &
                                       + cam_in%pwtgcell(J,17) &
                                       + cam_in%pwtgcell(J,19)
        State_Met%LandTypeFrac(1,J,34) = cam_in%pwtgcell(J, 5)
        State_Met%LandTypeFrac(1,J,36) = cam_in%pwtgcell(J,18) &
                                       + cam_in%pwtgcell(J,21)
        State_Met%LandTypeFrac(1,J,42) = cam_in%pwtgcell(J,14)
        State_Met%LandTypeFrac(1,J,43) = cam_in%pwtgcell(J,13)
        State_Met%LandTypeFrac(1,J,44) = cam_in%pwtgcell(J, 7) & 
                                       + cam_in%pwtgcell(J,15)
        !State_Met%LandTypeFrac(1,J,45) = cam_in%lwtgcell(J, 5)
        State_Met%LandTypeFrac(1,J,48) = cam_in%pwtgcell(J,11)

        State_Met%XLAI_NATIVE(1,J, 5)  = cam_in%lai(J, 4)
        State_Met%XLAI_NATIVE(1,J, 7)  = cam_in%lai(J, 6)
        State_Met%XLAI_NATIVE(1,J,10)  = cam_in%lai(J,12)
        State_Met%XLAI_NATIVE(1,J,17)  = cam_in%lai(J,10)
        State_Met%XLAI_NATIVE(1,J,22)  = cam_in%lai(J, 3)
        State_Met%XLAI_NATIVE(1,J,24)  = cam_in%lai(J, 9)
        State_Met%XLAI_NATIVE(1,J,26)  = cam_in%lai(J, 8)
        State_Met%XLAI_NATIVE(1,J,28)  = cam_in%lai(J, 2)
        State_Met%XLAI_NATIVE(1,J,31)  = cam_in%lai(J,20)
        State_Met%XLAI_NATIVE(1,J,32)  = cam_in%lai(J,16) &
                                       + cam_in%lai(J,17) &
                                       + cam_in%lai(J,19)
        State_Met%XLAI_NATIVE(1,J,34)  = cam_in%lai(J, 5)
        State_Met%XLAI_NATIVE(1,J,36)  = cam_in%lai(J,18) &
                                       + cam_in%lai(J,21)
        State_Met%XLAI_NATIVE(1,J,42)  = cam_in%lai(J,14)
        State_Met%XLAI_NATIVE(1,J,43)  = cam_in%lai(J,13)
        State_Met%XLAI_NATIVE(1,J,44)  = cam_in%lai(J, 7) & 
                                       + cam_in%lai(J,15)
        State_Met%XLAI_NATIVE(1,J,48)  = cam_in%lai(J,11)

        DO T = 2, NSURFTYPE
            State_Met%LandTypeFrac(1,J,T) = &
            State_Met%LandTypeFrac(1,J,T) * landFrac

            State_Met%XLAI_NATIVE(1,J,T) = &
            State_Met%XLAI_NATIVE(1,J,T) * landFrac

            ! Make sure that the land type fractions do not exceed 1
            IF ( State_Met%LandTypeFrac(1,J,T) > 1.0e+0_fp ) THEN
                State_Met%LandTypeFrac(1,J,T) = 1.0e+0_fp
            ELSEIF ( State_Met%LandTypeFrac(1,J,T) < 0.0e+0_fp ) THEN
                State_Met%LandTypeFrac(1,J,T) = 0.0e+0_fp
            ENDIF
        ENDDO

    ENDDO

#elif defined( CLM45 ) || defined( CLM50 )

    ! Mapping for CLM4.5/CLM5.0
    ! -----------------------------------|--------------------------------------
    !           Olson land type          |             CLM land type
    ! -----------------------------------|--------------------------------------
    ! Inland/sea water         (ID =  1) | Ocean fraction
    !                                    | Deeplake                  (LUID =  5)
    ! Urban                    (ID =  2) | Urban - Not Applied       (LUID =7-9)
    ! Low Sparse Grassland     (ID =  3) | 
    ! Coniferous Forest        (ID =  4) | 
    ! Deciduous Conifer Forest (ID =  5) | Needleleaf Deciduous Bor. (PAID =  3)
    ! Deciduous Broadleaf For. (ID =  6) |
    ! Evergreen Broadleaf For. (ID =  7) | Broadleaf Evergreen Temp. (PAID =  5)
    ! Tall Grasses and Shrubs  (ID =  8) | 
    ! Bare Desert              (ID =  9) | Not veg. \ Ice    (PAID = 0\LUID = 4)
    ! Upland Tundra            (ID = 10) | Broadleaf Deciduous Bore. (PAID = 11)
    ! Irrigated Grassland      (ID = 11) |
    ! Semi Desert              (ID = 12) | 
    ! Glacier ice              (ID = 13) | Land ice                  (LUID =  4)
    ! Wooded Wet Swamp         (ID = 14) | 
    ! -                        (ID = 15) | 
    ! -                        (ID = 16) | 
    ! Shrub Evergreen          (ID = 17) | Broadleaf Evergreen Shru. (PAID =  9)
    ! -                        (ID = 18) | 
    ! Shrub Deciduous          (ID = 19) | 
    ! Evergreen Forest and Fi. (ID = 20) | 
    ! Cool Rain Forest         (ID = 21) | 
    ! Conifer Boreal Forest    (ID = 22) | Needleleaf Evergreen Bor. (PAID =  2) 
    ! Cool Conifer Forest      (ID = 23) | 
    ! Cool Mixed Forest        (ID = 24) | Broadleaf Deciduous Bore. (PAID =  8)
    ! Mixed Forest             (ID = 25) | 
    ! Cool Broadleaf Forest    (ID = 26) | Broadleaf Deciduous Temp. (PAID =  7)
    ! Deciduous Broadleaf For. (ID = 27) | 
    ! Conifer Forest           (ID = 28) | Needleleaf Evergreen Tem. (PAID =  1)
    ! Montane Tropical Forests (ID = 29) | 
    ! Seasonal Tropical Fores. (ID = 30) | 
    ! Cool Crops and Towns     (ID = 31) |
    ! Crops and Town           (ID = 32) | C3 Crop                   (PAID = 15)
    !                                    | C3 Irrigated              (PAID = 16)
    ! Dry Tropical Woods       (ID = 33) | 
    ! Tropical Rainforest      (ID = 34) | Broadleaf Evergreen Trop. (PAID =  4)
    ! Tropical Degraded Forest (ID = 35) | 
    ! Corn and Beans Cropland  (ID = 36) | Corn                      (PAID = 17)
    !                                    | Irrigated Temperate Corn  (PAID = 18)
    !                                    | Spring Wheat              (PAID = 19)
    !                                    | Irrigated Spring Wheat    (PAID = 20)
    !                                    | Winter Wheat              (PAID = 21)
    !                                    | Irrigated Winter Wheat    (PAID = 22)
    !                                    | Temperated Soybean        (PAID = 23)
    !                                    | Irrigated Temperate Soyb. (PAID = 24)
    !                                    | Barley                    (PAID = 25)
    !                                    | Irrigated Barley          (PAID = 26)
    !                                    | Winter Barley             (PAID = 27)
    !                                    | Irrigated Winter Barley   (PAID = 28)
    !                                    | Rye                       (PAID = 29)
    !                                    | Irrigated Rye             (PAID = 30)
    !                                    | Winter Rye                (PAID = 31)
    !                                    | Irrigated Winter Rye      (PAID = 32)
    !                                    | Cassava                   (PAID = 33)
    !                                    | Irrigated Cassava         (PAID = 34)
    !                                    | Citrus                    (PAID = 35)
    !                                    | Irrigated Citrus          (PAID = 36)
    !                                    | Cocoa                     (PAID = 37)
    !                                    | Irrigated Cocoa           (PAID = 38)
    !                                    | Coffee                    (PAID = 39)
    !                                    | Irrigated Coffee          (PAID = 40)
    !                                    | Cotton                    (PAID = 41)
    !                                    | Irrigated Cotton          (PAID = 42)
    !                                    | Datepalm                  (PAID = 43)
    !                                    | Irrigated Datepalm        (PAID = 44)
    !                                    | Foddergrass               (PAID = 45)
    !                                    | Irrigated Foddergrass     (PAID = 46)
    !                                    | Grapes                    (PAID = 47)
    !                                    | Irrigated Grapes          (PAID = 48)
    !                                    | Groundnuts                (PAID = 49)
    !                                    | Irrigated Groundnuts      (PAID = 50)
    !                                    | Millet                    (PAID = 51)
    !                                    | Irrigated Millet          (PAID = 52)
    !                                    | Oilpalm                   (PAID = 53)
    !                                    | Irrigated Oilpalm         (PAID = 54)
    !                                    | Potatoes                  (PAID = 55)
    !                                    | Irrigated Potatoes        (PAID = 56)
    !                                    | Pulses                    (PAID = 57)
    !                                    | Irrigated Pulses          (PAID = 58)
    !                                    | Rapeseed                  (PAID = 59)
    !                                    | Irrigated Rapeseed        (PAID = 60)
    !                                    | Rice                      (PAID = 61)
    !                                    | Irrigated Rice            (PAID = 62)
    !                                    | Sorghum                   (PAID = 63)
    !                                    | Irrigated Sorghum         (PAID = 64)
    !                                    | Sugarbeet                 (PAID = 65)
    !                                    | Irrigated Sugarbeet       (PAID = 66)
    !                                    | Sugarcane                 (PAID = 67)
    !                                    | Irrigated Sugarcane       (PAID = 68)
    !                                    | Sunflower                 (PAID = 69)
    !                                    | Irrigated Sunflower       (PAID = 70)
    !                                    | Miscanthus                (PAID = 71)
    !                                    | Irrigated Miscanthus      (PAID = 72)
    !                                    | Switchgrass               (PAID = 73)
    !                                    | Irrigated Switchgrass     (PAID = 74)
    !                                    | Tropical Corn             (PAID = 75)
    !                                    | Irrigated Tropical Corn   (PAID = 76)
    !                                    | Tropical Soybean          (PAID = 77)
    !                                    | Irrigated Tropical Soybe. (PAID = 78)
    ! Rice Paddy and Field     (ID = 37) | 
    ! Hot Irrigated Cropland   (ID = 38) | 
    ! Cool Irrigated Cropland  (ID = 39) | 
    ! -                        (ID = 40) | 
    ! Cool Grasses and Shrubs  (ID = 41) | 
    ! Hot and Mild Grasses and (ID = 42) | C3 Non-Arctic Grass       (PAID = 13)
    ! Cold Grassland           (ID = 43) | C3 Arctic Grass           (PAID = 12)
    ! Savanna (Woods)          (ID = 44) | Broadleaf Deciduous Trop. (PAID =  6)
    !                                    | C4 Grass                  (PAID = 14)
    ! Mire, Bog, Fen           (ID = 45) | Wetland - Not Applied     (LUID =  6)
    ! Marsh Wetland            (ID = 46) |
    ! Mediterranean Scrub      (ID = 47) |
    ! Dry Woody Scrub          (ID = 48) | Broadleaf Deciduous Temp. (PAID = 10)
    ! -                        (ID = 49) | 
    ! -                        (ID = 50) | 
    ! -                        (ID = 51) | 
    ! Semi Desert Shrubs       (ID = 52) |
    ! Semi Desert Sage         (ID = 53) | 
    ! Barren Tundra            (ID = 54) | 
    ! Cool Southern Hemisphere (ID = 55) |
    ! Cool Fields and Woods    (ID = 56) | 
    ! Forest and Field         (ID = 57) | 
    ! Cool Forest and Field    (ID = 58) | 
    ! Fields and Woody Savanna (ID = 59) | 
    ! Succulent and Thorn Scr. (ID = 60) | 
    ! Small Leaf Mixed Woods   (ID = 61) | 
    ! Deciduous and Mixed Bor. (ID = 62) | 
    ! Narrow Conifers          (ID = 63) | 
    ! Wooded Tundra            (ID = 64) | 
    ! Heath Scrub              (ID = 65) | 
    ! -                        (ID = 66) | 
    ! -                        (ID = 67) | 
    ! -                        (ID = 68) | 
    ! -                        (ID = 69) | 
    ! Polar and Alpine Desert  (ID = 70) | 
    ! -                        (ID = 71) | 
    ! -                        (ID = 72) | 
    ! Mangrove                 (ID = 73) | 

    State_Met%LandTypeFrac(:,:,:) = 0.0e+0_fp
    State_Met%XLAI_NATIVE(:,:,:)  = 0.0e+0_fp

    DO J = 1, nY
        waterFrac = cam_in%ocnFrac(J)    + cam_in%iceFrac(J)    &
                  + cam_in%lwtgcell(J,5)
        landFrac  = 1.0e+0_fp - waterFrac

        ! Initialize fraction land for this grid cell
        State_Met%LandTypeFrac(1,J, 1) = waterFrac
        !State_Met%LandTypeFrac(1,J, 2) = cam_in%lwtgcell(J, 7) &
        !                               + cam_in%lwtgcell(J, 8) &
        !                               + cam_in%lwtgcell(J, 9)
        State_Met%LandTypeFrac(1,J, 5) = cam_in%pwtgcell(J, 4)
        State_Met%LandTypeFrac(1,J, 7) = cam_in%pwtgcell(J, 6)
        State_Met%LandTypeFrac(1,J, 9) = cam_in%pwtgcell(J, 1) &
                       * ( 1.0e+0_fp - cam_in%lwtgcell(J, 4) )
        State_Met%LandTypeFrac(1,J,10) = cam_in%pwtgcell(J,12)
        State_Met%LandTypeFrac(1,J,13) = cam_in%lwtgcell(J, 4)
        State_Met%LandTypeFrac(1,J,17) = cam_in%pwtgcell(J,10)
        State_Met%LandTypeFrac(1,J,22) = cam_in%pwtgcell(J, 3)
        State_Met%LandTypeFrac(1,J,24) = cam_in%pwtgcell(J, 9)
        State_Met%LandTypeFrac(1,J,26) = cam_in%pwtgcell(J, 8)
        State_Met%LandTypeFrac(1,J,28) = cam_in%pwtgcell(J, 2)
        State_Met%LandTypeFrac(1,J,32) = cam_in%pwtgcell(J,16) &
                                       + cam_in%pwtgcell(J,17)
        State_Met%LandTypeFrac(1,J,34) = cam_in%pwtgcell(J, 5)
        State_Met%LandTypeFrac(1,J,36) = cam_in%pwtgcell(J,18) &
                                       + cam_in%pwtgcell(J,19) &
                                       + cam_in%pwtgcell(J,20) &
                                       + cam_in%pwtgcell(J,21) &
                                       + cam_in%pwtgcell(J,22) &
                                       + cam_in%pwtgcell(J,23) &
                                       + cam_in%pwtgcell(J,24) &
                                       + cam_in%pwtgcell(J,25) &
                                       + cam_in%pwtgcell(J,26) &
                                       + cam_in%pwtgcell(J,27) &
                                       + cam_in%pwtgcell(J,28) &
                                       + cam_in%pwtgcell(J,29) &
                                       + cam_in%pwtgcell(J,30) &
                                       + cam_in%pwtgcell(J,31) &
                                       + cam_in%pwtgcell(J,32) &
                                       + cam_in%pwtgcell(J,33) &
                                       + cam_in%pwtgcell(J,34) &
                                       + cam_in%pwtgcell(J,35) &
                                       + cam_in%pwtgcell(J,36) &
                                       + cam_in%pwtgcell(J,37) &
                                       + cam_in%pwtgcell(J,38) &
                                       + cam_in%pwtgcell(J,39) &
                                       + cam_in%pwtgcell(J,40) &
                                       + cam_in%pwtgcell(J,41) &
                                       + cam_in%pwtgcell(J,42) &
                                       + cam_in%pwtgcell(J,43) &
                                       + cam_in%pwtgcell(J,44) &
                                       + cam_in%pwtgcell(J,45) &
                                       + cam_in%pwtgcell(J,46) &
                                       + cam_in%pwtgcell(J,47) &
                                       + cam_in%pwtgcell(J,48) &
                                       + cam_in%pwtgcell(J,49) &
                                       + cam_in%pwtgcell(J,50) &
                                       + cam_in%pwtgcell(J,51) &
                                       + cam_in%pwtgcell(J,52) &
                                       + cam_in%pwtgcell(J,53) &
                                       + cam_in%pwtgcell(J,54) &
                                       + cam_in%pwtgcell(J,55) &
                                       + cam_in%pwtgcell(J,56) &
                                       + cam_in%pwtgcell(J,57) &
                                       + cam_in%pwtgcell(J,58) &
                                       + cam_in%pwtgcell(J,59) &
                                       + cam_in%pwtgcell(J,60) &
                                       + cam_in%pwtgcell(J,61) &
                                       + cam_in%pwtgcell(J,62) &
                                       + cam_in%pwtgcell(J,63) &
                                       + cam_in%pwtgcell(J,64) &
                                       + cam_in%pwtgcell(J,65) &
                                       + cam_in%pwtgcell(J,66) &
                                       + cam_in%pwtgcell(J,67) &
                                       + cam_in%pwtgcell(J,68) &
                                       + cam_in%pwtgcell(J,69) &
                                       + cam_in%pwtgcell(J,70) &
                                       + cam_in%pwtgcell(J,71) &
                                       + cam_in%pwtgcell(J,72) &
                                       + cam_in%pwtgcell(J,73) &
                                       + cam_in%pwtgcell(J,74) &
                                       + cam_in%pwtgcell(J,75) &
                                       + cam_in%pwtgcell(J,76) &
                                       + cam_in%pwtgcell(J,77) &
                                       + cam_in%pwtgcell(J,78) &
                                       + cam_in%pwtgcell(J,79)
        State_Met%LandTypeFrac(1,J,42) = cam_in%pwtgcell(J,14)
        State_Met%LandTypeFrac(1,J,43) = cam_in%pwtgcell(J,13)
        State_Met%LandTypeFrac(1,J,44) = cam_in%pwtgcell(J, 7) &
                                       + cam_in%pwtgcell(J,15)
        !State_Met%LandTypeFrac(1,J,45) = cam_in%lwtgcell(J, 6)
        State_Met%LandTypeFrac(1,J,48) = cam_in%pwtgcell(J,11)

        State_Met%XLAI_NATIVE(1,J, 5)  = cam_in%lai(J, 4)
        State_Met%XLAI_NATIVE(1,J, 7)  = cam_in%lai(J, 6)
        State_Met%XLAI_NATIVE(1,J,10)  = cam_in%lai(J,12)
        State_Met%XLAI_NATIVE(1,J,17)  = cam_in%lai(J,10)
        State_Met%XLAI_NATIVE(1,J,22)  = cam_in%lai(J, 3)
        State_Met%XLAI_NATIVE(1,J,24)  = cam_in%lai(J, 9)
        State_Met%XLAI_NATIVE(1,J,26)  = cam_in%lai(J, 8)
        State_Met%XLAI_NATIVE(1,J,28)  = cam_in%lai(J, 2)
        State_Met%XLAI_NATIVE(1,J,32)  = cam_in%lai(J,16) &
                                       + cam_in%lai(J,17)
        State_Met%XLAI_NATIVE(1,J,34)  = cam_in%lai(J, 5)
        State_Met%XLAI_NATIVE(1,J,36)  = cam_in%lai(J,18) &
                                       + cam_in%lai(J,19) &
                                       + cam_in%lai(J,20) &
                                       + cam_in%lai(J,21) &
                                       + cam_in%lai(J,22) &
                                       + cam_in%lai(J,23) &
                                       + cam_in%lai(J,24) &
                                       + cam_in%lai(J,25) &
                                       + cam_in%lai(J,26) &
                                       + cam_in%lai(J,27) &
                                       + cam_in%lai(J,28) &
                                       + cam_in%lai(J,29) &
                                       + cam_in%lai(J,30) &
                                       + cam_in%lai(J,31) &
                                       + cam_in%lai(J,32) &
                                       + cam_in%lai(J,33) &
                                       + cam_in%lai(J,34) &
                                       + cam_in%lai(J,35) &
                                       + cam_in%lai(J,36) &
                                       + cam_in%lai(J,37) &
                                       + cam_in%lai(J,38) &
                                       + cam_in%lai(J,39) &
                                       + cam_in%lai(J,40) &
                                       + cam_in%lai(J,41) &
                                       + cam_in%lai(J,42) &
                                       + cam_in%lai(J,43) &
                                       + cam_in%lai(J,44) &
                                       + cam_in%lai(J,45) &
                                       + cam_in%lai(J,46) &
                                       + cam_in%lai(J,47) &
                                       + cam_in%lai(J,48) &
                                       + cam_in%lai(J,49) &
                                       + cam_in%lai(J,50) &
                                       + cam_in%lai(J,51) &
                                       + cam_in%lai(J,52) &
                                       + cam_in%lai(J,53) &
                                       + cam_in%lai(J,54) &
                                       + cam_in%lai(J,55) &
                                       + cam_in%lai(J,56) &
                                       + cam_in%lai(J,57) &
                                       + cam_in%lai(J,58) &
                                       + cam_in%lai(J,59) &
                                       + cam_in%lai(J,60) &
                                       + cam_in%lai(J,61) &
                                       + cam_in%lai(J,62) &
                                       + cam_in%lai(J,63) &
                                       + cam_in%lai(J,64) &
                                       + cam_in%lai(J,65) &
                                       + cam_in%lai(J,66) &
                                       + cam_in%lai(J,67) &
                                       + cam_in%lai(J,68) &
                                       + cam_in%lai(J,69) &
                                       + cam_in%lai(J,70) &
                                       + cam_in%lai(J,71) &
                                       + cam_in%lai(J,72) &
                                       + cam_in%lai(J,73) &
                                       + cam_in%lai(J,74) &
                                       + cam_in%lai(J,75) &
                                       + cam_in%lai(J,76) &
                                       + cam_in%lai(J,77) &
                                       + cam_in%lai(J,78) &
                                       + cam_in%lai(J,79)
        State_Met%XLAI_NATIVE(1,J,42)  = cam_in%lai(J,14)
        State_Met%XLAI_NATIVE(1,J,43)  = cam_in%lai(J,13)
        State_Met%XLAI_NATIVE(1,J,44)  = cam_in%lai(J, 7) &
                                       + cam_in%lai(J,15)
        State_Met%XLAI_NATIVE(1,J,48)  = cam_in%lai(J,11)

        DO T = 2, NSURFTYPE
            State_Met%LandTypeFrac(1,J,T) = &
            State_Met%LandTypeFrac(1,J,T) * landFrac

            State_Met%XLAI_NATIVE(1,J,T) = &
            State_Met%XLAI_NATIVE(1,J,T) * landFrac

            ! Make sure that the land type fractions do not exceed 1
            IF ( State_Met%LandTypeFrac(1,J,T) > 1.0e+0_fp ) THEN
                State_Met%LandTypeFrac(1,J,T) = 1.0e+0_fp
            ELSEIF ( State_Met%LandTypeFrac(1,J,T) < 0.0e+0_fp ) THEN
                State_Met%LandTypeFrac(1,J,T) = 0.0e+0_fp
            ENDIF
        ENDDO

    ENDDO

#else
    CALL endrun('Cannot figure out which version of CLM')
#endif

    END SUBROUTINE getLandTypes
!EOC
