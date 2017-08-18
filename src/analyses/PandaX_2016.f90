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
  INTEGER, PARAMETER :: NBINS = 0
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
      =       (/ 0.00033d0, 0.00049d0, 0.00072d0, 0.00101d0, 0.00137d0, &
      0.00178d0, 0.00222d0, 0.00271d0, 0.00331d0, 0.00406d0, 0.00501d0, &
      0.00620d0, 0.00728d0, 0.00887d0, 0.01055d0, 0.01221d0, 0.01422d0, &
      0.01656d0, 0.01918d0, 0.02220d0, 0.02566d0, 0.02961d0, 0.03413d0, &
      0.03927d0, 0.04510d0, 0.05165d0, 0.05894d0, 0.06693d0, 0.07554d0, &
      0.08470d0, 0.09448d0, 0.10476d0, 0.11551d0, 0.12680d0, 0.13870d0, &
      0.15114d0, 0.16408d0, 0.17732d0, 0.19094d0, 0.20497d0, 0.21931d0, &
      0.23382d0, 0.24853d0, 0.26330d0, 0.27781d0, 0.29203d0, 0.30584d0, &
      0.31923d0, 0.33211d0, 0.34437d0, 0.35565d0, 0.36588d0, 0.37560d0, &
      0.38519d0, 0.39462d0, 0.40300d0, 0.40944d0, 0.41476d0, 0.41965d0, &
      0.42474d0, 0.42916d0, 0.43282d0, 0.43588d0, 0.43901d0, 0.44137d0, &
      0.44284d0, 0.44379d0, 0.44442d0, 0.44466d0, 0.44568d0, 0.44602d0, &
      0.44486d0, 0.44321d0, 0.44032d0, 0.43611d0, 0.43047d0, 0.42314d0, &
      0.41384d0, 0.40229d0, 0.38809d0, 0.37134d0, 0.35267d0, 0.33291d0, &
      0.31204d0, 0.28997d0, 0.26706d0, 0.24325d0, 0.21851d0, 0.19279d0, &
      0.16606d0, 0.13956d0, 0.12110d0, 0.10394d0, 0.08823d0, 0.07409d0, &
      0.06164d0, 0.05107d0, 0.04250d0, 0.03549d0, 0.02964d0, 0.02427d0 /)
  ! Efficiencies array (2D)
  INTEGER, PARAMETER :: NELEM=1
  REAL*8, PARAMETER :: EFF(NELEM,NE,0:NBINS)                                   &
      = RESHAPE( (/ (/ EFF0(:) /) /),SHAPE(EFF))

  CALL SetDetector(D,mass=334.3d0,time=98.7d0,Nevents=(/3/),                &
                   background=(/4.8d0/),Nelem=NELEM,Zelem=(/54/),               &
                   NE=NE,E=E,Nbins=NBINS,eff=EFF,                   &
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
