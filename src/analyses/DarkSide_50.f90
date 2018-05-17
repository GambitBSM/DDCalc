MODULE DarkSide_50

!=======================================================================
! DarkSide-50 ANALYSIS ROUTINES
! Based upon arXiv:1802.07198.  
!=======================================================================

USE DDTypes
USE DDDetectors

IMPLICIT NONE

CONTAINS


!-----------------------------------------------------------------------
! Initializes a DetectorStruct to the DarkSide-50 analysis.
! 
FUNCTION DarkSide_50_Init() RESULT(D)

  IMPLICIT NONE
  TYPE(DetectorStruct) :: D
  INTEGER, PARAMETER :: NE = 71
  INTEGER, PARAMETER :: NBINS = 0
  REAL*8, PARAMETER :: EMIN = 40d0
  ! Efficiency curves energy tabulation points
  REAL*8, PARAMETER :: E(NE)                                            &
      =       (/ 40d0, 40.738d0, 41.6869d0, 42.658d0, 43.6516d0,        &
      44.6684d0, 45.7088d0, 46.7735d0, 47.863d0, 48.9779d0, 50.1187d0,  &
      51.2861d0, 52.4807d0, 53.7032d0, 54.9541d0, 56.2341d0, 57.544d0,  &
      58.8844d0, 60.256d0, 61.6595d0, 63.0957d0, 64.5654d0, 66.0693d0,  &
      67.6083d0, 69.1831d0, 70.7946d0, 72.4436d0, 74.131d0, 75.8578d0,  &
      77.6247d0, 79.4328d0, 81.2831d0, 83.1764d0, 85.1138d0, 87.0964d0, &
      89.1251d0, 91.2011d0, 93.3254d0, 95.4993d0, 97.7237d0, 100.d0,    &
      102.329d0, 104.713d0, 107.152d0, 109.648d0, 112.202d0, 114.815d0, &
      117.49d0, 120.226d0, 123.027d0, 125.893d0, 128.825d0, 131.826d0,  &
      134.896d0, 138.038d0, 141.254d0, 144.544d0, 147.911d0, 151.356d0, &
      154.882d0, 158.489d0, 162.181d0, 165.959d0, 169.824d0, 173.78d0,  &
      177.828d0, 181.97d0, 186.209d0, 190.546d0, 194.984d0, 200d0 /)
  ! Efficiency (total)
  REAL*8, PARAMETER :: EFF0(NE)                                         &
      =       (/ 0d0, 0d0, 0d0, 0d0, 0d0,                               &
      0.00133863d0, 0.00355493d0, 0.00582285d0, 0.00995072d0, 0.01431d0,&
      0.0286665d0, 0.0468955d0, 0.0906758d0, 0.135789d0, 0.185331d0,    &
      0.230552d0, 0.271346d0, 0.311883d0, 0.352708d0, 0.413193d0,       &
      0.46071d0, 0.488894d0, 0.521126d0, 0.536742d0, 0.548441d0,        &
      0.571357d0, 0.579175d0, 0.585033d0, 0.588919d0, 0.592042d0,       &
      0.594637d0, 0.594953d0, 0.595276d0, 0.595606d0, 0.595944d0,       &
      0.59599d0, 0.595956d0, 0.595921d0, 0.595884d0, 0.595848d0,        &
      0.59581d0, 0.595772d0, 0.595732d0, 0.595692d0, 0.595651d0,        &
      0.595609d0, 0.595426d0, 0.594963d0, 0.594488d0, 0.594003d0,       &
      0.593506d0, 0.592997d0, 0.592477d0, 0.592367d0, 0.593462d0,       &
      0.594582d0, 0.595728d0, 0.596901d0, 0.598101d0, 0.598067d0,       &
      0.598007d0, 0.597946d0, 0.597884d0, 0.597133d0, 0.596241d0,       &
      0.595328d0, 0.594383d0, 0.593275d0, 0.592141d0, 0.590462d0,       &
      0.571295d0 /)
  ! Efficiencies array (2D)
  REAL*8, PARAMETER :: EFF(NE,0:NBINS)                                   &
      = RESHAPE( (/ EFF0(:) /), SHAPE(EFF))

  CALL SetDetector(D,mass=36.9d0,time=532.4d0,Nevents_tot=0,            &
                   Backgr_tot=0.09d0,Nelem=1,Zelem=(/18/),              &
                   NE=NE,E=E,Nbins=NBINS,eff_all=EFF,                   &
                   Emin=EMIN)
  
END FUNCTION  


! C++ interface wrapper
INTEGER(KIND=C_INT) FUNCTION C_DarkSide_50_Init() &
 BIND(C,NAME='C_DDCalc_darkside_50_init') 
  USE ISO_C_BINDING, only: C_BOOL, C_INT
  IMPLICIT NONE
  N_Detectors = N_Detectors + 1
  ALLOCATE(Detectors(N_Detectors)%p)
  Detectors(N_Detectors)%p = DarkSide_50_Init()
  C_DarkSide_50_Init = N_Detectors
END FUNCTION


END MODULE
