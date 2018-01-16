MODULE DummyExp

!=======================================================================
! Dummy experiment ANALYSIS ROUTINES
!=======================================================================

USE DDTypes
USE DDDetectors

IMPLICIT NONE

CONTAINS


!-----------------------------------------------------------------------
! Initializes a DetectorStruct of a dummy experiment
!
FUNCTION DummyExp_Init() RESULT(D)

  IMPLICIT NONE
  TYPE(DetectorStruct) :: D
  INTEGER, PARAMETER :: NE = 151
  INTEGER, PARAMETER :: NEeff = 8
  INTEGER, PARAMETER :: NBINS = 2
  REAL*8, ALLOCATABLE ::  EFF(:,:)
  REAL*8, PARAMETER :: Eeff(NEeff)                                      &
      =       (/ 0.1d0, 1.0d0,  2.0d0,  3.0d0,  8.0d0,  &
      9.0d0,  10.0d0, 100.00d0 /)
  REAL*8, PARAMETER :: EFF0(NEeff)                                      &
      =       (/ 0.d0, 0.d0, 0.5d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0 /)
  REAL*8, PARAMETER :: EFF1(NEeff)                                      &
      =       (/ 0.d0, 0.d0, 0.5d0, 1.0d0, 1.0d0, 0.5d0, 0.0d0, 0.0d0 /)
  REAL*8, PARAMETER :: EFF2(NEeff)                                      &
      =       (/ 0.d0, 0.d0, 0.0d0, 0.0d0, 0.0d0, 0.5d0, 1.0d0, 1.0d0 /)

  REAL*8, PARAMETER :: EFFeff(NEeff,0:NBINS)                            &
      = RESHAPE( (/ (/ EFF0(:), EFF1(:), EFF2(:) /) /) ,SHAPE(EFFeff))

  REAL*8, PARAMETER :: E(NE)                                            &
      =       (/ 0.10000d0, 0.10471d0, 0.10965d0, 0.11482d0, 0.12023d0, &
      0.12589d0, 0.13183d0, 0.13804d0, 0.14454d0, 0.15136d0, 0.15849d0, &
      0.16596d0, 0.17378d0, 0.18197d0, 0.19055d0, 0.19953d0, 0.20893d0, &
      0.21878d0, 0.22909d0, 0.23988d0, 0.25119d0, 0.26303d0, 0.27542d0, &
      0.28840d0, 0.30200d0, 0.31623d0, 0.33113d0, 0.34674d0, 0.36308d0, &
      0.38019d0, 0.39811d0, 0.41687d0, 0.43652d0, 0.45709d0, 0.47863d0, &
      0.50119d0, 0.52481d0, 0.54954d0, 0.57544d0, 0.60256d0, 0.63096d0, &
      0.66069d0, 0.69183d0, 0.72444d0, 0.75858d0, 0.79433d0, 0.83176d0, &
      0.87096d0, 0.91201d0, 0.95499d0, 1.0000d0,  1.0471d0,  1.0965d0,  &
      1.1482d0,  1.2023d0,  1.2589d0,  1.3183d0,  1.3804d0,  1.4454d0,  &
      1.5136d0,  1.5849d0,  1.6596d0,  1.7378d0,  1.8197d0,  1.9055d0,  &
      1.9953d0,  2.0893d0,  2.1878d0,  2.2909d0,  2.3988d0,  2.5119d0,  &
      2.6303d0,  2.7542d0,  2.8840d0,  3.0200d0,  3.1623d0,  3.3113d0,  &
      3.4674d0,  3.6308d0,  3.8019d0,  3.9811d0,  4.1687d0,  4.3652d0,  &
      4.5709d0,  4.7863d0,  5.0119d0,  5.2481d0,  5.4954d0,  5.7544d0,  &
      6.0256d0,  6.3096d0,  6.6069d0,  6.9183d0,  7.2444d0,  7.5858d0,  &
      7.9433d0,  8.3176d0,  8.7096d0,  9.1201d0,  9.5499d0, 10.000d0,   &
     10.471d0,  10.965d0,  11.482d0,  12.023d0,  12.589d0,  13.183d0,   &
     13.804d0,  14.454d0,  15.136d0,  15.849d0,  16.596d0,  17.378d0,   &
     18.197d0,  19.055d0,  19.953d0,  20.893d0,  21.878d0,  22.909d0,   &
     23.988d0,  25.119d0,  26.303d0,  27.542d0,  28.840d0,  30.200d0,   &
     31.623d0,  33.113d0,  34.674d0,  36.308d0,  38.019d0,  39.811d0,   &
     41.687d0,  43.652d0,  45.709d0,  47.863d0,  50.119d0,  52.481d0,   &
     54.954d0,  57.544d0,  60.256d0,  63.096d0,  66.069d0,  69.183d0,   &
     72.444d0,  75.858d0,  79.433d0,  83.176d0,  87.096d0,  91.201d0,   &
     95.499d0, 100.00d0 /)

  CALL RetabulateEfficiency(NEeff,Eeff,NBINS,EFFeff,NE,E,EFF)

  WRITE(*,*) EFF(:,0)
  WRITE(*,*) EFF(:,1)
  WRITE(*,*) EFF(:,2)

  CALL SetDetector(D,mass=118d0,time=85.3d0,Nevents_tot=11,               &
                   Backgr_tot = 1.3d0, &
                   Nelem=1,Zelem=(/54/),               &
                   NE=NE,E=E,Nbins=NBINS,eff_all=EFF)
  D%eff_file = '[DummyExp]'
  
END FUNCTION


! C++ interface wrapper
INTEGER(KIND=C_INT) FUNCTION C_DummyExp_Init() &
 BIND(C,NAME='C_DDCalc_dummyexp_init') 
  USE ISO_C_BINDING, only: C_BOOL, C_INT
  IMPLICIT NONE
  N_Detectors = N_Detectors + 1
  ALLOCATE(Detectors(N_Detectors)%p)
  Detectors(N_Detectors)%p = DummyExp_Init()
  C_DummyExp_Init = N_Detectors
END FUNCTION


END MODULE
