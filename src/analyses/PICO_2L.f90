MODULE PICO_2L

!=======================================================================
! PICO_2L ANALYSIS ROUTINES
! Based upon http://arxiv.org/pdf/1601.03729v2.pdf .  
!=======================================================================

USE DDTypes
USE DDDetectors

IMPLICIT NONE

CONTAINS


!-----------------------------------------------------------------------
! Initializes a DetectorStruct to the PICO_2L analysis.
! 
! This analysis includes only fluorine!
!
FUNCTION PICO_2L_Init() RESULT(D)

  IMPLICIT NONE
  TYPE(DetectorStruct) :: D
  INTEGER, PARAMETER :: NE = 101
  INTEGER, PARAMETER :: NBINS = 0
  REAL*8, PARAMETER :: EMIN = 3.3d0
  ! Efficiency curves energy tabulation points
  REAL*8, PARAMETER :: E(NE)                                            &
      =       (/ 3.3d0,     3.39712d0, 3.49709d0, 3.6d0,     3.70595d0, &
      3.81501d0, 3.92728d0, 4.04286d0, 4.16184d0, 4.28431d0, 4.4104d0,  &
      4.54019d0, 4.6738d0,  4.81135d0, 4.95294d0, 5.0987d0,  5.24875d0, &
      5.40322d0, 5.56223d0, 5.72592d0, 5.89442d0, 6.06789d0, 6.24646d0, &
      6.43029d0, 6.61953d0, 6.81433d0, 7.01487d0, 7.22131d0, 7.43383d0, &
      7.6526d0,  7.8778d0,  8.10964d0, 8.3483d0,  8.59398d0, 8.84689d0, &
      9.10724d0, 9.37526d0, 9.65117d0, 9.93519d0, 10.2276d0, 10.5286d0, &
      10.8384d0, 11.1574d0, 11.4857d0, 11.8237d0, 12.1717d0, 12.5299d0, &
      12.8986d0, 13.2782d0, 13.669d0,  14.0712d0, 14.4853d0, 14.9116d0, &
      15.3505d0, 15.8022d0, 16.2673d0, 16.746d0,  17.2388d0, 17.7461d0, &
      18.2684d0, 18.806d0,  19.3594d0, 19.9292d0, 20.5157d0, 21.1194d0, &
      21.7409d0, 22.3807d0, 23.0394d0, 23.7174d0, 24.4154d0, 25.1339d0, &
      25.8736d0, 26.635d0,  27.4188d0, 28.2258d0, 29.0564d0, 29.9115d0, &
      30.7918d0, 31.6979d0, 32.6308d0, 33.5911d0, 34.5796d0, 35.5973d0, &
      36.6448d0, 37.7233d0, 38.8334d0, 39.9762d0, 41.1527d0, 42.3638d0, &
      43.6105d0, 44.8939d0, 46.2151d0, 47.5752d0, 48.9752d0, 50.4165d0, &
      51.9002d0, 53.4276d0, 54.9999d0, 56.6185d0, 58.2847d0, 60.d0 /)
  ! Efficiency (total)
  REAL*8, PARAMETER :: EFF0(NE)                                         &
      =       (/ 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, &
      0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, &
      0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, &
      0.00000d0, 0.01260d0, 0.18145d0, 0.32697d0, 0.44480d0, 0.57556d0, &
      0.75532d0, 0.90435d0, 1.00000d0, 1.00000d0, 1.00000d0, 1.00000d0, &
      1.00000d0, 1.00000d0, 1.00000d0, 1.00000d0, 1.00000d0, 1.00000d0, &
      1.00000d0, 1.00000d0, 1.00000d0, 1.00000d0, 1.00000d0, 1.00000d0, &
      1.00000d0, 1.00000d0, 1.00000d0, 1.00000d0, 1.00000d0, 1.00000d0, &
      1.00000d0, 1.00000d0, 1.00000d0, 1.00000d0, 1.00000d0, 1.00000d0, &
      1.00000d0, 1.00000d0, 1.00000d0, 1.00000d0, 1.00000d0, 1.00000d0, &
      1.00000d0, 1.00000d0, 1.00000d0, 1.00000d0, 1.00000d0, 1.00000d0, &
      1.00000d0, 1.00000d0, 1.00000d0, 1.00000d0, 1.00000d0, 1.00000d0, &
      1.00000d0, 1.00000d0, 1.00000d0, 1.00000d0, 1.00000d0, 1.00000d0, &
      1.00000d0, 1.00000d0, 1.00000d0, 1.00000d0, 1.00000d0, 1.00000d0, &
      1.00000d0, 1.00000d0, 1.00000d0, 1.00000d0, 1.00000d0, 1.00000d0, &
      1.00000d0, 1.00000d0, 1.00000d0, 1.00000d0, 1.00000d0, 1.00000d0, &
      1.00000d0, 1.00000d0, 1.00000d0, 1.00000d0, 1.00000d0, 1.00000d0 /)
  ! Efficiencies array (2D)
  REAL*8, PARAMETER :: EFF(NE,0:NBINS)                                   &
      = RESHAPE( (/ EFF0(:) /),SHAPE(EFF))

  ! The fiducial mass is reduced to account for the fluorine fraction.
  ! NOTE: PICO_2L does not attempt background subtraction

  CALL SetDetector(D,mass=1.57d0,time=66.3d0,Nevents_tot=1,              &
                   Backgr_tot=0.0d0,Nelem=1,Zelem=(/9/),                 &
                   NE=NE,E=E,Nbins=NBINS,eff_all=EFF,                    &
                   Emin=EMIN)
  D%eff_file = '[PICO_2L]'
  
END FUNCTION


! C++ interface wrapper
INTEGER(KIND=C_INT) FUNCTION C_PICO_2L_Init() &
 BIND(C,NAME='C_DDCalc_pico_2l_init') 
  USE ISO_C_BINDING, only: C_BOOL, C_INT
  IMPLICIT NONE
  N_Detectors = N_Detectors + 1
  ALLOCATE(Detectors(N_Detectors)%p)
  Detectors(N_Detectors)%p = PICO_2L_Init()
  C_PICO_2L_Init = N_Detectors
END FUNCTION


END MODULE
