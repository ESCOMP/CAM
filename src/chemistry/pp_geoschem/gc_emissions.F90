!================================================================================================
! This is the "GEOS-Chem" chemistry emissions interface
!================================================================================================
module GC_Emissions_Mod

  use Shr_kind_mod,        only : r8 => shr_kind_r8
  use Spmd_utils,          only : MasterProc, myCPU=>iam, nCPUs=>npes
  use Cam_logfile,         only : iulog
  use Cam_abortutils,      only : endrun

  use Chem_mods,           only : NTracers
  use Chem_mods,           only : TracerNames
  use Chem_mods,           only : Map2GC

  use Tracer_data,         only : trfld,trfile

  IMPLICIT NONE

  TYPE :: Emission
      INTEGER                    :: Spc_Ndx
      REAL(r8)                   :: MW
      REAL(r8)                   :: Scalefactor
      CHARACTER(LEN=256)         :: Filename
      CHARACTER(LEN=16)          :: Species
      CHARACTER(LEN=8)           :: Units
      INTEGER                    :: Nsectors
      CHARACTER(LEN=32), POINTER :: Sectors(:)
      TYPE(trfld),       POINTER :: Fields(:)
      TYPE(trfile)               :: File
  ENDTYPE Emission

  PRIVATE

  PUBLIC :: GC_Emissions_Init
  PUBLIC :: GC_Emissions_Calc
  PUBLIC :: GC_Emissions_Final

  ! Stand-in: emissions
  TYPE(Emission), ALLOCATABLE :: Emissions(:)
  INTEGER                     :: N_Emis_Files

!================================================================================================
contains
!================================================================================================

  subroutine GC_Emissions_Init

  INTEGER :: Ierr

  N_Emis_Files=1
  ALLOCATE(Emissions(N_Emis_Files), STAT=IERR)
  IF (IERR.NE.0) CALL ENDRUN('Could not allocate GC emissions')

  end subroutine GC_Emissions_Init

  subroutine GC_Emissions_Calc(Eflx)

  ! Emissions in kg/m2/s
  ! Dimensions: [N columns x K levels x C constituents ]
  REAL(r8), INTENT(OUT) :: EFlx(:,:,:)
  INTEGER               :: I_Trc, I_Emis

  EFlx(:,:,:) = 0.0e+0_r8
  DO I_Emis = 1, N_Emis_Files
      ! Read emissions file
      DO I_Trc = 1, NTracers
      ENDDO
  ENDDO

  end subroutine GC_Emissions_Calc

  subroutine GC_Emissions_Final
  IF (ALLOCATED(Emissions)) DEALLOCATE(Emissions)
  end subroutine GC_Emissions_Final

  end module GC_Emissions_Mod
