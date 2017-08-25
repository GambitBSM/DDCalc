MODULE DDExHelp

USE DDExperiments
USE DDCommandLine
USE DDDetectors

IMPLICIT NONE
PRIVATE

INTERFACE DDCalc_ChooseAnalysisCommandLine
  MODULE PROCEDURE ChooseAnalysisCommandLine
END INTERFACE

PUBLIC :: DDCalc_ChooseAnalysisCommandLine
PUBLIC :: DDCalc_InitDetector, C_DDCalc_InitDetector

CONTAINS

!-----------------------------------------------------------------------
! Initialise default detector
! 
FUNCTION DDCalc_InitDetector(intervals) RESULT(Detector)

  IMPLICIT NONE
  LOGICAL, INTENT(IN), OPTIONAL :: intervals
  TYPE(DetectorStruct) :: Detector
  LOGICAL :: intervals0
  
  ! Include sub-intervals?
  intervals0 = .TRUE.
  IF (PRESENT(intervals)) intervals0 = intervals
  IF (GetLongArg('no-intervals')) intervals0 = .FALSE.

  Detector = DDCalc_InitDefaultDetector(intervals0)

END FUNCTION

!-----------------------------------------------------------------------
! C/C++ wrapper for DDCalc_InitDetector
!
INTEGER(KIND=C_INT) FUNCTION C_DDCalc_InitDetector() &
 BIND(C,NAME='C_DDExperiments_ddcalc_initdetector') 
  USE ISO_C_BINDING, only: C_INT
  IMPLICIT NONE
  N_Detectors = N_Detectors + 1
  IF (N_Detectors .GT. Max_Detectors) stop 'DDCalc: Max_Detectors exceeded.&
   Please run FreeDetectors or modify Max_Detectors in DDTypes.f90.'
  ALLOCATE(Detectors(N_Detectors)%p)
  Detectors(N_Detectors)%p = DDCalc_InitDetector()
  C_DDCalc_InitDetector = N_Detectors
END FUNCTION


!-----------------------------------------------------------------------
! Parses command-line arguments and chooses the requested analysis.
!
! Optional input arguments:
!   eff_file    File from which efficiencies shoud be read.  If not
!               given or set to '', an internal default will be used
!               (the Xenon1T 2017 result).  Non-empty value takes precedence
!               over --file option.
!   intervals   Specify if sub-intervals should be loaded from the
!               efficiency file (if available); otherwise, only the
!               full interval is used for calculations (default: true).
!               Superceded by --no-intervals setting.
! 
FUNCTION ChooseAnalysisCommandLine(eff_file,intervals) RESULT(Detector)

  CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: eff_file
  LOGICAL, INTENT(IN), OPTIONAL :: intervals
  TYPE(DetectorStruct) :: Detector
  CHARACTER(LEN=1024) :: eff_file0
  LOGICAL :: intervals0

  ! File name ('' for internal default)
  IF (.NOT. GetLongArgString('file',eff_file0)) eff_file0 = ''
  IF (PRESENT(eff_file)) THEN
    IF (TRIM(eff_file) .NE. '') eff_file0 = eff_file
  END IF
  
  ! Include sub-intervals?
  intervals0 = .TRUE.
  IF (PRESENT(intervals)) intervals0 = intervals
  IF (GetLongArg('no-intervals')) intervals0 = .FALSE.

  ! Experiment-specific settings
  Detector = AvailableAnalyses(intervals0)
  
  ! If given, load efficiency file
  IF (TRIM(eff_file0) .NE. '') THEN
    CALL SetDetector(Detector,eff_file=eff_file0,intervals=intervals0)
  END IF

END FUNCTION


END MODULE
