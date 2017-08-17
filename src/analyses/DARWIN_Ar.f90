MODULE DARWIN_Ar

!=======================================================================
! DARWIN ARGON ANALYSIS ROUTINES
! Based upon a DARWIN argon-based analysis (2015 estimated parameters)
! [15MM.NNNNN].  Zero events assumed in the analysis region.
!=======================================================================

USE DDTypes
USE DDDetectors

IMPLICIT NONE

CONTAINS


!-----------------------------------------------------------------------
! Initializes a DetectorStruct to the DARWIN Argon projection.
! 
! The efficiencies used here were generated using TPCMC.
! 
! Required input arguments:
!     intervals   Indicates if sub-intervals should be included
! 
FUNCTION DARWIN_Ar_Init(intervals) RESULT(D)

  IMPLICIT NONE
  TYPE(DetectorStruct) :: D
  LOGICAL, INTENT(IN) :: intervals
  INTEGER, PARAMETER :: NE = 151
  INTEGER, PARAMETER :: NEFF = 1
  ! Efficiency curves energy tabulation points
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
  ! Efficiency (total)
  REAL*8, PARAMETER :: EFF0(NE)                                         &
      =       (/ 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, &
      0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, &
      0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, &
      0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, &
      1.00000d-6,1.00000d-6,3.00000d-6,5.00000d-6,4.00000d-6,4.00000d-6,&
      6.00000d-6,9.00000d-6,1.20000d-5,1.70000d-5,3.60000d-5,5.40000d-5,&
      7.00000d-5,1.06000d-4,1.48000d-4,2.48000d-4,3.74000d-4,5.11000d-4,&
      6.61000d-4,9.53000d-4,1.32400d-3,1.91900d-3,2.57200d-3,3.42800d-3,&
      4.78200d-3,6.31100d-3,8.37200d-3,1.11320d-2,1.43600d-2,1.88560d-2,&
      2.39690d-2,3.04640d-2,3.85180d-2,4.79960d-2,5.92610d-2,7.20370d-2,&
      8.75540d-2,1.04670d-1,1.24360d-1,1.44260d-1,1.67600d-1,1.92550d-1,&
      2.19830d-1,2.45220d-1,2.71260d-1,2.95940d-1,3.20590d-1,3.44590d-1,&
      3.66380d-1,3.85290d-1,4.02490d-1,4.17410d-1,4.30070d-1,4.40120d-1,&
      4.47510d-1,4.55340d-1,4.59120d-1,4.63940d-1,4.66210d-1,4.67680d-1,&
      4.70120d-1,4.70570d-1,4.72180d-1,4.73500d-1,4.75110d-1,4.75610d-1,&
      4.76660d-1,4.77150d-1,4.77350d-1,4.77890d-1,4.78970d-1,4.79310d-1,&
      4.80410d-1,4.81470d-1,4.81580d-1,4.83040d-1,4.82920d-1,4.82570d-1,&
      4.83820d-1,4.84220d-1,4.85050d-1,4.84690d-1,4.85190d-1,4.84650d-1,&
      4.84290d-1,4.83720d-1,4.84160d-1,4.83160d-1,4.81950d-1,4.79380d-1,&
      4.76860d-1,4.73220d-1,4.62790d-1,4.43970d-1,4.14520d-1,3.73500d-1,&
      3.19880d-1,2.59970d-1,1.96280d-1,1.36860d-1,8.72900d-2,5.09020d-2,&
      2.65530d-2,1.23380d-2,5.26700d-3,1.86600d-3,6.02000d-4,1.50000d-4,&
      4.60000d-5,3.00000d-6,3.00000d-6,0.00000d0, 0.00000d0, 0.00000d0, &
      0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, &
      0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, &
      0.00000d0, 0.00000d0 /)
  ! Efficiency (first and only interval)
  REAL*8, PARAMETER :: EFF1(NE) = EFF0
  ! Efficiencies array (2D)
  REAL*8, PARAMETER :: EFF(NE,0:NEFF)                                   &
      = RESHAPE( (/ EFF0(:), EFF1(:) /) ,SHAPE(EFF))
  
  ! One call for all settings.
  ! Most of these _must_ be there to ensure everything get initialized.
  CALL SetDetector(D,mass=20d3,time=2d0*365d0,Nevents=0,                &
                   background=0.5d0,Nelem=1,Zelem=(/18/),               &
                   NE=NE,E=E,Neff=NEFF,eff=EFF,                   &
                   intervals=intervals)
  D%eff_file = '[DARWIN Ar 2015]'
  
END FUNCTION


! C++ interface wrapper
INTEGER(KIND=C_INT) FUNCTION C_DARWIN_Ar_Init(intervals) &
 BIND(C,NAME='C_DDCalc_darwin_ar_init') 
  USE ISO_C_BINDING, only: C_BOOL, C_INT
  IMPLICIT NONE
  LOGICAL(KIND=C_BOOL), INTENT(IN) :: intervals
  N_Detectors = N_Detectors + 1
  ALLOCATE(Detectors(N_Detectors)%p)
  Detectors(N_Detectors)%p = DARWIN_Ar_Init(LOGICAL(intervals))
  C_DARWIN_Ar_Init = N_Detectors
END FUNCTION


END MODULE
