MODULE PICO_60_F

!=======================================================================
! PICO_60 ANALYSIS ROUTINES
! Based upon https://arxiv.org/pdf/1510.07754v3.pdf .  
!=======================================================================

USE DDTypes
USE DDDetectors

IMPLICIT NONE

CONTAINS


!-----------------------------------------------------------------------
! Initializes a DetectorStruct to the PICO_60_F analysis.
!
! This analysis includes only fluorine, which is sufficient for SD scattering.
!
FUNCTION PICO_60_F_Init(intervals) RESULT(D)

  IMPLICIT NONE
  TYPE(DetectorStruct) :: D
  LOGICAL, INTENT(IN) :: intervals
  INTEGER, PARAMETER :: NE = 101
  INTEGER, PARAMETER :: NEFF = 1
  REAL*8, PARAMETER :: EMIN = 5.0d0
  ! Efficiency curves energy tabulation points
  REAL*8, PARAMETER :: E(NE)                                            &
      =       (/ 5.d0,      5.15205d0, 5.30873d0, 5.47017d0, 5.63652d0, &
      5.80793d0, 5.98455d0, 6.16655d0, 6.35408d0, 6.54731d0, 6.74641d0, &
      6.95158d0, 7.16298d0, 7.38081d0, 7.60526d0, 7.83654d0, 8.07486d0, &
      8.32042d0, 8.57345d0, 8.83417d0, 9.10282d0, 9.37964d0, 9.66488d0, &
      9.9588d0,  10.2617d0, 10.5737d0, 10.8953d0, 11.2266d0, 11.568d0,  &
      11.9198d0, 12.2823d0, 12.6558d0, 13.0407d0, 13.4372d0, 13.8459d0, &
      14.2669d0, 14.7008d0, 15.1479d0, 15.6085d0, 16.0832d0, 16.5723d0, &
      17.0762d0, 17.5955d0, 18.1306d0, 18.682d0,  19.2501d0, 19.8355d0, &
      20.4387d0, 21.0603d0, 21.7007d0, 22.3607d0, 23.0407d0, 23.7414d0, &
      24.4633d0, 25.2073d0, 25.9739d0, 26.7637d0, 27.5776d0, 28.4163d0, &
      29.2804d0, 30.1709d0, 31.0884d0, 32.0338d0, 33.008d0,  34.0118d0, &
      35.0461d0, 36.1119d0, 37.21d0,   38.3416d0, 39.5076d0, 40.7091d0, &
      41.947d0,  43.2227d0, 44.5371d0, 45.8915d0, 47.2871d0, 48.7251d0, &
      50.2069d0, 51.7337d0, 53.3069d0, 54.928d0,  56.5984d0, 58.3196d0, &
      60.0931d0, 61.9206d0, 63.8036d0, 65.744d0,  67.7433d0, 69.8034d0, &
      71.9261d0, 74.1134d0, 76.3673d0, 78.6896d0, 81.0826d0, 83.5484d0, &
      86.0892d0, 88.7072d0, 91.4048d0, 94.1845d0, 97.0487d0, 100.d0 /)
  ! The efficiency is obtained by convoluting the (time-dependent) energy
  ! threshold with the efficiency given for 13.6 keV.
  ! Efficiency (total)
  REAL*8, PARAMETER :: EFF0(NE)                                         &
      =       (/ 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, &
      0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, &
      0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, &
      0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, &
      0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, &
      0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, &
      0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00005d0, 0.00029d0, &
      0.00079d0, 0.00159d0, 0.00251d0, 0.00340d0, 0.00434d0, 0.00553d0, &
      0.00729d0, 0.01005d0, 0.01368d0, 0.01815d0, 0.02356d0, 0.03136d0, &
      0.04299d0, 0.05765d0, 0.07300d0, 0.08921d0, 0.10916d0, 0.13664d0, &
      0.17110d0, 0.21335d0, 0.26507d0, 0.32146d0, 0.37453d0, 0.42395d0, &
      0.47787d0, 0.54257d0, 0.60605d0, 0.66718d0, 0.73220d0, 0.78894d0, &
      0.84033d0, 0.87927d0, 0.90937d0, 0.93311d0, 0.95289d0, 0.96765d0, &
      0.97782d0, 0.98745d0, 0.99396d0, 0.99725d0, 0.99938d0, 0.99998d0, &
      1.00000d0, 1.00000d0, 1.00000d0, 1.00000d0, 1.00000d0, 1.00000d0, &
      1.00000d0, 1.00000d0, 1.00000d0, 1.00000d0, 1.00000d0, 1.00000d0, &
      1.00000d0, 1.00000d0, 1.00000d0, 1.00000d0, 1.00000d0, 1.00000d0 /)
  ! Efficiency (first and only interval)
  REAL*8, PARAMETER :: EFF1(NE) = EFF0
  ! Efficiencies array (2D)
  REAL*8, PARAMETER :: EFF(NE,0:NEFF)                                   &
      = RESHAPE( (/ EFF0(:), EFF1(:) /) ,SHAPE(EFF))

  ! The fiducial mass is reduced to account for the fluorine fraction.
  ! It is furthermore reduced by a trial factor of 1.8.
  ! NOTE: PICO-2L does not attempt background subtraction (to be implemented)

  CALL SetDetector(D,mass=2.32d0,time=92.8d0,Nevents=0,                 &
                   background=0.0d0,Nelem=1,Zelem=(/9/),                &
                   NEeff=NE,Eeff=E,Neff=NEFF,eff=EFF,                   &
                   intervals=intervals,Emin=EMIN)
  D%eff_file = '[PICO_60 F]'
  
END FUNCTION


! C++ interface wrapper
INTEGER(KIND=C_INT) FUNCTION C_PICO_60_F_Init(intervals) &
 BIND(C,NAME='C_DDCalc_pico_60_f_init') 
  USE ISO_C_BINDING, only: C_BOOL, C_INT
  IMPLICIT NONE
  LOGICAL(KIND=C_BOOL), INTENT(IN) :: intervals
  N_Detectors = N_Detectors + 1
  ALLOCATE(Detectors(N_Detectors)%p)
  Detectors(N_Detectors)%p = PICO_60_F_Init(LOGICAL(intervals))
  C_PICO_60_F_Init = N_Detectors
END FUNCTION


END MODULE
