MODULE LUX_2016_prelim

!=======================================================================
! LUX 2016 prelim ANALYSIS ROUTINES
! Based upon the IDM talk on the LUX 2016 analysis.  
!=======================================================================

USE DDTypes
USE DDDetectors

IMPLICIT NONE

CONTAINS


!-----------------------------------------------------------------------
! Initializes a DetectorStruct to the LUX 2016 analysis.
! 
! The following assumptions and simplifications are made:
! 1. We consider only events below the mean of the nuclear recoil band.
! 2. Since most events that fall below the nuclear recoil band are observed either at very energies (S1 < 3) 
! or at very large energies (S1 > 33), we limit the search window to 3 < S1 < 33, leading to only one observed event.
! 3. To estimate the signal acceptance in this window, we perform a Monte Carlo simulation of nuclear recoils, using 
! all available information on the photon and electron yield, the conversion of primary photons/electrons to S1 and S2
! and the S2 detection efficiency. We also include a very simple model for fluctuations of primary photons/electrons 
! and of the detection probability.
! 4. To estimate the background in this window, we assume that background events due to leakage from the electron recoil band
! are distributed uniformly in S1 and that all other backgrounds fall off exponentially in the same way as the spectrum from
! calibration with neutrons. This leads to a prediction of 2.3 background events in the search window.
! 
FUNCTION LUX_2016_prelim_Init(intervals) RESULT(D)

  IMPLICIT NONE
  TYPE(DetectorStruct) :: D
  LOGICAL, INTENT(IN) :: intervals
  INTEGER, PARAMETER :: NE = 101
  INTEGER, PARAMETER :: NEFF = 1
  REAL*8, PARAMETER :: EMIN = 1.1d0
  ! Efficiency curves energy tabulation points
  REAL*8, PARAMETER :: E(NE)                                            &
      =       (/ 1.1d0,     1.1428d0,  1.18726d0, 1.23345d0, 1.28143d0, &
      1.33129d0, 1.38308d0, 1.43689d0, 1.49279d0, 1.55087d0, 1.6112d0,  &
      1.67389d0, 1.73901d0, 1.80667d0, 1.87695d0, 1.94998d0, 2.02584d0, &
      2.10465d0, 2.18654d0, 2.2716d0,  2.35998d0, 2.45179d0, 2.54718d0, &
      2.64628d0, 2.74923d0, 2.85619d0, 2.96731d0, 3.08275d0, 3.20269d0, &
      3.32729d0, 3.45673d0, 3.59122d0, 3.73093d0, 3.87608d0, 4.02688d0, &
      4.18355d0, 4.34631d0, 4.5154d0,  4.69107d0, 4.87358d0, 5.06318d0, &
      5.26017d0, 5.46481d0, 5.67742d0, 5.8983d0,  6.12777d0, 6.36617d0, &
      6.61385d0, 6.87116d0, 7.13848d0, 7.4162d0,  7.70472d0, 8.00448d0, &
      8.31589d0, 8.63942d0, 8.97553d0, 9.32472d0, 9.6875d0,  10.0644d0, &
      10.4559d0, 10.8627d0, 11.2853d0, 11.7244d0, 12.1805d0, 12.6544d0, &
      13.1467d0, 13.6582d0, 14.1896d0, 14.7416d0, 15.3151d0, 15.911d0,  &
      16.53d0,   17.1731d0, 17.8412d0, 18.5353d0, 19.2564d0, 20.0056d0, &
      20.7839d0, 21.5925d0, 22.4326d0, 23.3053d0, 24.212d0,  25.1539d0, &
      26.1326d0, 27.1492d0, 28.2055d0, 29.3028d0, 30.4428d0, 31.6272d0, &
      32.8577d0, 34.136d0,  35.464d0,  36.8437d0, 38.2771d0, 39.7663d0, &
      41.3134d0, 42.9207d0, 44.5905d0, 46.3253d0, 48.1276d0, 50.d0 /)
  ! LOWER 50% NR BAND >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! Efficiency (total)
  REAL*8, PARAMETER :: EFF0(NE)                                         &
      =       (/ 0.00108d0, 0.00185d0, 0.00195d0, 0.00299d0, 0.00344d0, &
      0.00424d0, 0.00489d0, 0.00592d0, 0.00791d0, 0.00881d0, 0.01080d0, &
      0.01261d0, 0.01491d0, 0.01694d0, 0.02067d0, 0.02445d0, 0.02769d0, &
      0.03315d0, 0.03869d0, 0.04019d0, 0.04754d0, 0.04883d0, 0.05939d0, &
      0.06460d0, 0.07848d0, 0.07967d0, 0.08888d0, 0.10340d0, 0.10810d0, &
      0.12270d0, 0.13540d0, 0.14270d0, 0.15580d0, 0.17690d0, 0.19550d0, &
      0.20230d0, 0.21600d0, 0.22520d0, 0.23640d0, 0.25800d0, 0.27290d0, &
      0.27200d0, 0.28550d0, 0.28930d0, 0.29520d0, 0.30530d0, 0.31920d0, &
      0.31550d0, 0.32850d0, 0.34090d0, 0.34900d0, 0.36890d0, 0.37210d0, &
      0.39450d0, 0.41610d0, 0.43530d0, 0.44670d0, 0.44380d0, 0.45980d0, &
      0.46240d0, 0.46800d0, 0.47240d0, 0.47230d0, 0.47730d0, 0.46940d0, &
      0.46300d0, 0.46760d0, 0.46760d0, 0.46870d0, 0.46120d0, 0.47040d0, &
      0.47270d0, 0.47650d0, 0.46340d0, 0.45050d0, 0.46170d0, 0.45180d0, &
      0.45280d0, 0.44910d0, 0.44680d0, 0.43770d0, 0.43300d0, 0.41960d0, &
      0.39410d0, 0.37610d0, 0.31950d0, 0.27650d0, 0.22310d0, 0.16820d0, &
      0.12230d0, 0.06659d0, 0.03409d0, 0.01530d0, 0.00730d0, 0.00280d0, &
      0.00080d0, 0.00000d0, 0.00010d0, 0.00000d0, 0.00000d0, 0.00000d0 /)
  ! Efficiency (first and only interval)
  REAL*8, PARAMETER :: EFF1(NE) = EFF0
  ! Efficiencies array (2D)
  REAL*8, PARAMETER :: EFF(NE,0:NEFF)                                   &
      = RESHAPE( (/ EFF0(:), EFF1(:) /) ,SHAPE(EFF))

  ! Note that LUX 2016 doesn't have a well-defined fiducial mass (it varies with time).
  ! We calculate the average mass from the known duration (332 days) and the known exposure (33500 kg-days).

  CALL SetDetector(D,mass=100.9d0,time=332.0d0,Nevents=1,               &
                   background=2.3d0,Nelem=1,Zelem=(/54/),               &
                   NEeff=NE,Eeff=E,Neff=NEFF,eff=EFF,                   &
                   intervals=intervals,Emin=EMIN)
  D%eff_file = '[LUX 2016 prelim]'
  
END FUNCTION


! C++ interface wrapper
INTEGER(KIND=C_INT) FUNCTION C_LUX_2016_prelim_Init(intervals) &
 BIND(C,NAME='C_DDCalc_lux_2016_prelim_init') 
  USE ISO_C_BINDING, only: C_BOOL, C_INT
  IMPLICIT NONE
  LOGICAL(KIND=C_BOOL), INTENT(IN) :: intervals
  N_Detectors = N_Detectors + 1
  ALLOCATE(Detectors(N_Detectors)%p)
  Detectors(N_Detectors)%p = LUX_2016_prelim_Init(LOGICAL(intervals))
  C_LUX_2016_prelim_Init = N_Detectors
END FUNCTION


END MODULE
