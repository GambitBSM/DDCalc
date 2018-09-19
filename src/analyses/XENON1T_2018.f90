MODULE XENON1T_2018

!=======================================================================
! XENON1T 2018 ANALYSIS ROUTINES
! Based upon arXiv:1805.12562.  
!=======================================================================

USE DDTypes
USE DDDetectors

IMPLICIT NONE

CONTAINS


!-----------------------------------------------------------------------
! Initializes a DetectorStruct to the XENON1T 2018 analysis.
! 
FUNCTION XENON1T_2018_Init() RESULT(D)

  IMPLICIT NONE
  TYPE(DetectorStruct) :: D
  INTEGER, PARAMETER :: NE = 101
  INTEGER, PARAMETER :: NBINS = 2
  REAL*8, PARAMETER :: EMIN = 1.5d0
  ! Efficiency curves energy tabulation points

  REAL*8, PARAMETER :: E(NE)                                            &
      =       (/ 1.50000d0, 1.55637d0, 1.61485d0, 1.67553d0, 1.73850d0, 1.80382d0, &
      1.87161d0, 1.94194d0, 2.01491d0, 2.09063d0, 2.16919d0, 2.25070d0, &
      2.33528d0, 2.42303d0, 2.51408d0, 2.60856d0, 2.70658d0, 2.80829d0, &
      2.91382d0, 3.02331d0, 3.13692d0, 3.25480d0, 3.37710d0, 3.50401d0, &
      3.63568d0, 3.77230d0, 3.91405d0, 4.06114d0, 4.21374d0, 4.37209d0, &
      4.53638d0, 4.70684d0, 4.88372d0, 5.06724d0, 5.25765d0, 5.45522d0, &
      5.66021d0, 5.87291d0, 6.09360d0, 6.32258d0, 6.56017d0, 6.80669d0, &
      7.06247d0, 7.32786d0, 7.60322d0, 7.88893d0, 8.18538d0, 8.49297d0, &
      8.81211d0, 9.14325d0, 9.48683d0, 9.84333d0, 10.21320d0, 10.59700d0, &
      10.99520d0, 11.40840d0, 11.83710d0, 12.28190d0, 12.74340d0, &
      13.22230d0, 13.71920d0, 14.23470d0, 14.76960d0, 15.32460d0, &
      15.90050d0, 16.49800d0, 17.11790d0, 17.76120d0, 18.42860d0, &
      19.12110d0, 19.83960d0, 20.58510d0, 21.35870d0, 22.16130d0, &
      22.99410d0, 23.85810d0, 24.75470d0, 25.68490d0, 26.65010d0, &
      27.65150d0, 28.69060d0, 29.76870d0, 30.88730d0, 32.04800d0, &
      33.25230d0, 34.50180d0, 35.79830d0, 37.14360d0, 38.53930d0, &
      39.98750d0, 41.49020d0, 43.04930d0, 44.66700d0, 46.34540d0, &
      48.08700d0, 49.89400d0, 51.76890d0, 53.71420d0, 55.73270d0, &
      57.82700d0, 60.00000d0 /)

  ! The efficiency is based on Fig. 1 of arXiv:1805.12562 with two additional correction factors:
  ! 1. A factor of 0.475 accounts for the fraction of nuclear recoil events that fall into the reference region.
  ! 2. A factor of 0.5 accounts for the fraction of WIMP events that fall into the inner 0.65 t of the fiducial volume.

  ! Efficiency (total)
  REAL*8, PARAMETER :: EFF0(NE)                                         &
      =       (/ 0.00000d0, 0.00000d0, 0.00000d0, 0.00010d0, 0.00023d0, 0.00036d0, &
      0.00050d0, 0.00092d0, 0.00137d0, 0.00184d0, 0.00245d0, 0.00321d0, &
      0.00398d0, 0.00494d0, 0.00605d0, 0.00721d0, 0.00866d0, 0.01022d0, &
      0.01195d0, 0.01395d0, 0.01607d0, 0.01857d0, 0.02116d0, 0.02422d0, &
      0.02742d0, 0.03109d0, 0.03501d0, 0.03936d0, 0.04405d0, 0.04907d0, &
      0.05453d0, 0.06027d0, 0.06633d0, 0.07267d0, 0.07925d0, 0.08601d0, &
      0.09301d0, 0.10007d0, 0.10721d0, 0.11439d0, 0.12132d0, 0.12827d0, &
      0.13463d0, 0.14090d0, 0.14661d0, 0.15203d0, 0.15709d0, 0.16157d0, &
      0.16585d0, 0.16957d0, 0.17294d0, 0.17604d0, 0.17875d0, 0.18142d0, &
      0.18378d0, 0.18622d0, 0.18834d0, 0.19042d0, 0.19231d0, 0.19403d0, &
      0.19561d0, 0.19699d0, 0.19821d0, 0.19927d0, 0.20013d0, 0.20091d0, &
      0.20144d0, 0.20193d0, 0.20232d0, 0.20263d0, 0.20292d0, 0.20317d0, &
      0.20339d0, 0.20357d0, 0.20375d0, 0.20392d0, 0.20411d0, 0.20432d0, &
      0.20456d0, 0.20476d0, 0.20481d0, 0.20462d0, 0.20403d0, 0.20288d0, &
      0.20099d0, 0.19814d0, 0.19368d0, 0.18726d0, 0.17826d0, 0.16540d0, &
      0.14854d0, 0.12729d0, 0.10362d0, 0.07914d0, 0.05637d0, 0.03763d0, &
      0.02364d0, 0.01386d0, 0.00744d0, 0.00346d0, 0.00135d0 /)

  ! The energy range is divided into two bins: [3, 35] and [35, 70]

  REAL*8, PARAMETER :: EFF1(NE)                                         &
      =       (/ 0.00000d0, 0.00000d0, 0.00000d0, 0.00010d0, 0.00023d0, 0.00036d0, &
      0.00050d0, 0.00092d0, 0.00137d0, 0.00184d0, 0.00245d0, 0.00321d0, &
      0.00398d0, 0.00494d0, 0.00605d0, 0.00721d0, 0.00866d0, 0.01022d0, &
      0.01195d0, 0.01395d0, 0.01607d0, 0.01857d0, 0.02116d0, 0.02422d0, &
      0.02742d0, 0.03109d0, 0.03501d0, 0.03936d0, 0.04405d0, 0.04907d0, &
      0.05453d0, 0.06027d0, 0.06633d0, 0.07267d0, 0.07925d0, 0.08601d0, &
      0.09301d0, 0.10007d0, 0.10721d0, 0.11439d0, 0.12132d0, 0.12827d0, &
      0.13463d0, 0.14090d0, 0.14661d0, 0.15203d0, 0.15709d0, 0.16157d0, &
      0.16585d0, 0.16957d0, 0.17294d0, 0.17604d0, 0.17875d0, 0.18142d0, &
      0.18378d0, 0.18622d0, 0.18834d0, 0.19042d0, 0.19231d0, 0.19403d0, &
      0.19560d0, 0.19697d0, 0.19815d0, 0.19913d0, 0.19981d0, 0.20020d0, &
      0.19997d0, 0.19898d0, 0.19673d0, 0.19261d0, 0.18591d0, 0.17592d0, &
      0.16210d0, 0.14441d0, 0.12342d0, 0.10043d0, 0.07722d0, 0.05586d0, &
      0.03769d0, 0.02350d0, 0.01341d0, 0.00694d0, 0.00314d0, 0.00125d0, &
      0.00044d0, 0.00014d0, 0.00004d0, 0.00001d0, 0.00000d0, 0.00000d0, &
      0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, &
      0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0 /)

  ! Efficiencies array (2D)
  REAL*8, PARAMETER :: EFF(NE,0:NBINS)                                   &
      = RESHAPE( (/ EFF0(:), EFF1(:), EFF0(:) - EFF1(:) /), SHAPE(EFF))

! The background distributions are estimates based on arXiv:1512.07501
! Zero events were observed in the first bin, 2 events in the second bin

  CALL SetDetector(D,mass=1300.0d0,time=278.8d0,Nevents_bin=[0,0],  &
                   Backgr_bin=[0.46d0,0.34d0],Nelem=1,Zelem=(/54/), &
                   NE=NE,E=E,Nbins=NBINS,eff_all=EFF,               &
                   Emin=EMIN)

!  Uncomment the lines below for an unbinned analysis of XENON1T
!  CALL SetDetector(D,mass=1300.0d0,time=278.8d0,Nevents_tot=2,     &
!                   Backgr_tot=0.8d0,Nelem=1,Zelem=(/54/),          &
!                   NE=NE,E=E,Nbins=NBINS,eff_all=EFF,              &
!                   Emin=EMIN)
   D%eff_file = '[XENON1T 2018]'
  
END FUNCTION  


! C++ interface wrapper
INTEGER(KIND=C_INT) FUNCTION C_XENON1T_2018_Init() &
 BIND(C,NAME='C_DDCalc_xenon1t_2018_init') 
  USE ISO_C_BINDING, only: C_BOOL, C_INT
  IMPLICIT NONE
  N_Detectors = N_Detectors + 1
  ALLOCATE(Detectors(N_Detectors)%p)
  Detectors(N_Detectors)%p = XENON1T_2018_Init()
  C_XENON1T_2018_Init = N_Detectors
END FUNCTION


END MODULE

