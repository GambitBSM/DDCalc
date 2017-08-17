MODULE LUX_2015

!=======================================================================
! LUX 2015 ANALYSIS ROUTINES
! Based upon http://arxiv.org/pdf/1512.03506v3.pdf .  
!=======================================================================

USE DDTypes
USE DDDetectors

IMPLICIT NONE

CONTAINS


!-----------------------------------------------------------------------
! Initializes a DetectorStruct to the LUX 2015 analysis.
! 
! LUX observes one event at high energies, which reduces the sensitivity
! to heavy WIMPs but plays no role for light WIMPs.
! To reflect this, we reduce the LUX search window to recoil energies
! below 30 keV, such that the observed event is outside the search window
! We assume that this has negligible impact on the expected background,
! so we take 0.64 events, as quoted in http://arxiv.org/pdf/1310.8214v2.pdf .
!
FUNCTION LUX_2015_Init(intervals) RESULT(D)

  IMPLICIT NONE
  TYPE(DetectorStruct) :: D
  LOGICAL, INTENT(IN) :: intervals
  INTEGER, PARAMETER :: NE = 101
  INTEGER, PARAMETER :: NEFF = 1
  REAL*8, PARAMETER :: EMIN = 1.1d0
  ! Efficiency curves energy tabulation points
  REAL*8, PARAMETER :: E(NE)                                            &
      =       (/ 1.1d0,     1.13697d0, 1.17519d0, 1.21469d0, 1.25551d0, &
      1.29771d0, 1.34133d0, 1.38642d0, 1.43302d0, 1.48118d0, 1.53097d0, &
      1.58242d0, 1.63561d0, 1.69059d0, 1.74741d0, 1.80614d0, 1.86685d0, &
      1.9296d0,  1.99445d0, 2.06149d0, 2.13078d0, 2.2024d0,  2.27642d0, &
      2.35294d0, 2.43202d0, 2.51377d0, 2.59826d0, 2.68559d0, 2.77585d0, &
      2.86915d0, 2.96559d0, 3.06527d0, 3.1683d0,  3.27479d0, 3.38486d0, &
      3.49863d0, 3.61622d0, 3.73777d0, 3.8634d0,  3.99325d0, 4.12747d0, &
      4.2662d0,  4.40959d0, 4.55781d0, 4.711d0,   4.86934d0, 5.03301d0, &
      5.20218d0, 5.37703d0, 5.55776d0, 5.74456d0, 5.93765d0, 6.13722d0, &
      6.3435d0,  6.55671d0, 6.77709d0, 7.00488d0, 7.24032d0, 7.48368d0, &
      7.73522d0, 7.99521d0, 8.26394d0, 8.5417d0,  8.8288d0,  9.12555d0, &
      9.43227d0, 9.7493d0,  10.077d0,  10.4157d0, 10.7658d0, 11.1276d0, &
      11.5016d0, 11.8882d0, 12.2878d0, 12.7008d0, 13.1277d0, 13.569d0,  &
      14.025d0,  14.4964d0, 14.9837d0, 15.4873d0, 16.0078d0, 16.5459d0, &
      17.102d0,  17.6768d0, 18.271d0,  18.8851d0, 19.5199d0, 20.1759d0, &
      20.8541d0, 21.555d0,  22.2795d0, 23.0284d0, 23.8024d0, 24.6024d0, &
      25.4293d0, 26.284d0,  27.1675d0, 28.0806d0, 29.0244d0, 30.d0 /)
  ! For the efficiency there is an extra factor 1/2 to consider only events
  ! below the nuclear recoil band and an extra factor 18^2/20^2 to consider
  ! only events in the inner fiducial volume.
  ! Efficiency (total)
  REAL*8, PARAMETER :: EFF0(NE)                                         &
      =       (/ 0.00150d0, 0.00188d0, 0.00233d0, 0.00293d0, 0.00355d0, &
      0.00420d0, 0.00528d0, 0.00642d0, 0.00760d0, 0.00881d0, 0.01100d0, &
      0.01330d0, 0.01597d0, 0.01969d0, 0.02355d0, 0.02745d0, 0.03129d0, &
      0.03526d0, 0.03936d0, 0.04504d0, 0.05242d0, 0.06005d0, 0.06806d0, &
      0.07659d0, 0.08541d0, 0.09452d0, 0.10394d0, 0.11548d0, 0.12758d0, &
      0.14010d0, 0.15336d0, 0.16775d0, 0.18263d0, 0.19801d0, 0.21391d0, &
      0.22284d0, 0.23120d0, 0.23983d0, 0.24876d0, 0.25799d0, 0.26909d0, &
      0.28090d0, 0.29310d0, 0.30571d0, 0.31875d0, 0.32667d0, 0.33469d0, &
      0.34299d0, 0.35156d0, 0.36042d0, 0.36265d0, 0.36419d0, 0.36578d0, &
      0.36742d0, 0.36922d0, 0.37352d0, 0.37796d0, 0.38256d0, 0.38731d0, &
      0.39222d0, 0.39729d0, 0.40109d0, 0.40114d0, 0.40119d0, 0.40124d0, &
      0.40129d0, 0.40134d0, 0.40140d0, 0.40145d0, 0.40150d0, 0.40155d0, &
      0.40160d0, 0.40165d0, 0.40170d0, 0.40175d0, 0.40181d0, 0.40187d0, &
      0.40215d0, 0.40249d0, 0.40284d0, 0.40320d0, 0.40357d0, 0.40396d0, &
      0.40436d0, 0.40477d0, 0.40460d0, 0.40370d0, 0.40276d0, 0.40180d0, &
      0.40080d0, 0.39977d0, 0.39870d0, 0.39760d0, 0.39646d0, 0.39528d0, &
      0.39522d0, 0.39662d0, 0.39806d0, 0.39955d0, 0.40110d0, 0.40269d0 /)
  ! Efficiency (first and only interval)
  REAL*8, PARAMETER :: EFF1(NE) = EFF0
  ! Efficiencies array (2D)
  REAL*8, PARAMETER :: EFF(NE,0:NEFF)                                   &
      = RESHAPE( (/ EFF0(:), EFF1(:) /) ,SHAPE(EFF))

  CALL SetDetector(D,mass=118.0d0,time=85.3d0,Nevents=0,                &
                   background=0.64d0,Nelem=1,Zelem=(/54/),              &
                   NE=NE,E=E,Neff=NEFF,eff=EFF,                   &
                   intervals=intervals,Emin=EMIN)
  D%eff_file = '[LUX 2015]'
  
END FUNCTION


! C++ interface wrapper
INTEGER(KIND=C_INT) FUNCTION C_LUX_2015_Init(intervals) &
 BIND(C,NAME='C_DDCalc_lux_2015_init') 
  USE ISO_C_BINDING, only: C_BOOL, C_INT
  IMPLICIT NONE
  LOGICAL(KIND=C_BOOL), INTENT(IN) :: intervals
  N_Detectors = N_Detectors + 1
  ALLOCATE(Detectors(N_Detectors)%p)
  Detectors(N_Detectors)%p = LUX_2015_Init(LOGICAL(intervals))
  C_LUX_2015_Init = N_Detectors
END FUNCTION


END MODULE
