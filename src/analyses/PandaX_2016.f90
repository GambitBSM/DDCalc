MODULE PandaX_2016

!=======================================================================
! PandaX 2016 ANALYSIS ROUTINES
! Based upon arXiv:1607.07400.  
!=======================================================================

USE DDTypes
USE DDDetectors

IMPLICIT NONE

CONTAINS


!-----------------------------------------------------------------------
! Initializes a DetectorStruct to the PandaX 2016 analysis.
! 
FUNCTION PandaX_2016_Init(intervals) RESULT(D)

  IMPLICIT NONE
  TYPE(DetectorStruct) :: D
  LOGICAL, INTENT(IN) :: intervals
  INTEGER, PARAMETER :: NE = 101
  INTEGER, PARAMETER :: NEFF = 1
  REAL*8, PARAMETER :: EMIN = 1.1d0
  ! Efficiency curves energy tabulation points
  REAL*8, PARAMETER :: E(NE)                                            &
      =       (/ 1.1d0,     1.14488d0, 1.19159d0, 1.24021d0, 1.29081d0, & 
      1.34348d0, 1.39829d0, 1.45534d0, 1.51472d0, 1.57653d0, 1.64085d0, &
      1.7078d0,  1.77748d0, 1.85d0,    1.92548d0, 2.00404d0, 2.08581d0, &
      2.17091d0, 2.25948d0, 2.35167d0, 2.44762d0, 2.54749d0, 2.65143d0, &
      2.75961d0, 2.8722d0,  2.98939d0, 3.11136d0, 3.2383d0,  3.37043d0, &
      3.50794d0, 3.65107d0, 3.80004d0, 3.95508d0, 4.11645d0, 4.28441d0, &
      4.45921d0, 4.64115d0, 4.83051d0, 5.0276d0,  5.23273d0, 5.44623d0, &
      5.66844d0, 5.89972d0, 6.14043d0, 6.39097d0, 6.65172d0, 6.92312d0, &
      7.20558d0, 7.49958d0, 7.80557d0, 8.12404d0, 8.45551d0, 8.8005d0,  &
      9.15956d0, 9.53328d0, 9.92224d0, 10.3271d0, 10.7484d0, 11.187d0,  &
      11.6434d0, 12.1185d0, 12.6129d0, 13.1275d0, 13.6631d0, 14.2206d0, &
      14.8008d0, 15.4047d0, 16.0332d0, 16.6874d0, 17.3683d0, 18.0769d0, &
      18.8144d0, 19.5821d0, 20.381d0,  21.2126d0, 22.0781d0, 22.9789d0, &
      23.9165d0, 24.8923d0, 25.9079d0, 26.9649d0, 28.0651d0, 29.2102d0, &
      30.402d0,  31.6424d0, 32.9335d0, 34.2772d0, 35.6757d0, 37.1313d0, &
      38.6463d0, 40.2231d0, 41.8642d0, 43.5723d0, 45.3501d0, 47.2004d0, &
      49.1262d0, 51.1306d0, 53.2168d0, 55.3881d0, 57.6479d0, 60.d0 /)
  ! LOWER 50% NR BAND >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! Efficiency (total)
  REAL*8, PARAMETER :: EFF0(NE)                                         &
      =       (/ 0.01266d0, 0.01454d0, 0.01650d0, 0.01853d0, 0.02065d0, &
      0.02286d0, 0.02515d0, 0.02754d0, 0.03002d0, 0.03392d0, 0.03800d0, &
      0.04224d0, 0.04666d0, 0.05125d0, 0.05604d0, 0.06179d0, 0.07072d0, &
      0.08002d0, 0.08970d0, 0.09977d0, 0.10814d0, 0.11521d0, 0.12256d0, &
      0.13021d0, 0.13818d0, 0.14635d0, 0.15484d0, 0.16368d0, 0.17288d0, &
      0.18244d0, 0.19182d0, 0.20158d0, 0.21174d0, 0.22231d0, 0.23273d0, &
      0.24349d0, 0.25469d0, 0.26634d0, 0.27605d0, 0.28537d0, 0.29508d0, &
      0.30518d0, 0.31569d0, 0.32495d0, 0.33433d0, 0.34408d0, 0.35423d0, &
      0.36480d0, 0.37152d0, 0.37841d0, 0.38558d0, 0.39304d0, 0.39901d0, &
      0.40495d0, 0.41114d0, 0.41600d0, 0.41963d0, 0.42342d0, 0.42736d0, &
      0.43185d0, 0.43670d0, 0.44175d0, 0.44662d0, 0.44820d0, 0.44985d0, &
      0.45156d0, 0.45334d0, 0.45386d0, 0.45302d0, 0.45215d0, 0.45124d0, &
      0.45008d0, 0.44876d0, 0.44739d0, 0.44427d0, 0.43766d0, 0.43079d0, &
      0.41966d0, 0.40455d0, 0.38633d0, 0.36645d0, 0.34565d0, 0.32399d0, &
      0.30592d0, 0.28741d0, 0.26493d0, 0.23642d0, 0.20674d0, 0.17510d0, &
      0.14163d0, 0.11482d0, 0.08883d0, 0.06538d0, 0.04839d0, 0.03733d0, &
      0.02654d0, 0.01851d0, 0.01182d0, 0.00747d0, 0.00483d0, 0.00338d0 /)
  ! Efficiency (first and only interval)
  REAL*8, PARAMETER :: EFF1(NE) = EFF0
  ! Efficiencies array (2D)
  REAL*8, PARAMETER :: EFF(NE,0:NEFF)                                   &
      = RESHAPE( (/ EFF0(:), EFF1(:) /) ,SHAPE(EFF))

  CALL SetDetector(D,mass=334.3d0,time=98.7d0,Nevents=3,                &
                   background=4.8d0,Nelem=1,Zelem=(/54/),               &
                   NEeff=NE,Eeff=E,Neff=NEFF,eff=EFF,                   &
                   intervals=intervals,Emin=EMIN)
  D%eff_file = '[PandaX 2016]'
  
END FUNCTION


! C++ interface wrapper
INTEGER(KIND=C_INT) FUNCTION C_PandaX_2016_Init(intervals) &
 BIND(C,NAME='C_DDCalc_pandax_2016_init') 
  USE ISO_C_BINDING, only: C_BOOL, C_INT
  IMPLICIT NONE
  LOGICAL(KIND=C_BOOL), INTENT(IN) :: intervals
  N_Detectors = N_Detectors + 1
  ALLOCATE(Detectors(N_Detectors)%p)
  Detectors(N_Detectors)%p = PandaX_2016_Init(LOGICAL(intervals))
  C_PandaX_2016_Init = N_Detectors
END FUNCTION


END MODULE
