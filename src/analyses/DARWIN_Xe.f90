MODULE DARWIN_Xe

!=======================================================================
! DARWIN XENON ANALYSIS ROUTINES
! Based upon a DARWIN xenon-based analysis (2015 estimated parameters)
! [15MM.NNNNN].  Zero events assumed in the analysis region.
!=======================================================================

USE DDTypes
USE DDDetectors

IMPLICIT NONE

CONTAINS


!-----------------------------------------------------------------------
! Initializes a DetectorStruct to the DARWIN Xenon projection.
! 
! The efficiencies used here were generated using TPCMC.
! 
! Required input arguments:
!     intervals   Indicates if sub-intervals should be included
! 
FUNCTION DARWIN_Xe_Init(intervals) RESULT(D)

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
      0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, &
      0.00000d0, 0.00000d0, 0.00000d0, 1.00000d-5,2.00000d-5,0.00000d0, &
      1.00000d-5,2.00000d-5,1.00000d-5,3.00000d-5,6.00000d-5,4.00000d-5,&
      1.60000d-4,2.50000d-4,3.30000d-4,4.30000d-4,5.80000d-4,8.30000d-4,&
      1.29000d-3,1.73000d-3,2.61000d-3,3.38000d-3,4.71000d-3,6.67000d-3,&
      8.93000d-3,1.18500d-2,1.66700d-2,2.11800d-2,2.76600d-2,3.70100d-2,&
      4.62800d-2,5.57800d-2,7.03000d-2,8.48500d-2,1.00920d-1,1.19770d-1,&
      1.40370d-1,1.59900d-1,1.83670d-1,2.13730d-1,2.34620d-1,2.60470d-1,&
      2.90490d-1,3.10700d-1,3.37130d-1,3.59870d-1,3.81220d-1,3.97590d-1,&
      4.13050d-1,4.29350d-1,4.33960d-1,4.43460d-1,4.50340d-1,4.54710d-1,&
      4.58700d-1,4.64210d-1,4.63840d-1,4.65210d-1,4.67910d-1,4.69700d-1,&
      4.71550d-1,4.70310d-1,4.72880d-1,4.73510d-1,4.77040d-1,4.76320d-1,&
      4.78180d-1,4.74830d-1,4.77390d-1,4.79300d-1,4.80260d-1,4.81980d-1,&
      4.79980d-1,4.83400d-1,4.85070d-1,4.83660d-1,4.86320d-1,4.83460d-1,&
      4.84910d-1,4.85100d-1,4.81940d-1,4.83840d-1,4.83400d-1,4.78140d-1,&
      4.77990d-1,4.70180d-1,4.48760d-1,4.18240d-1,3.67580d-1,3.03320d-1,&
      2.31100d-1,1.54890d-1,9.42000d-2,5.02400d-2,2.18100d-2,8.18000d-3,&
      2.44000d-3,7.40000d-4,1.70000d-4,1.00000d-5,1.00000d-5,0.00000d0, &
      0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, &
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
  CALL SetDetector(D,mass=12d3,time=2d0*365d0,Nevents=0,                &
                   background=0.5d0,Nelem=1,Zelem=(/54/),               &
                   NEeff=NE,Eeff=E,Neff=NEFF,eff=EFF,                   &
                   intervals=intervals)
  D%eff_file = '[DARWIN Xe 2015]'
  
END FUNCTION


! C++ interface wrapper
INTEGER(KIND=C_INT) FUNCTION C_DARWIN_Xe_Init(intervals) &
 BIND(C,NAME='C_DDCalc_darwin_xe_init') 
  USE ISO_C_BINDING, only: C_BOOL, C_INT
  IMPLICIT NONE
  LOGICAL(KIND=C_BOOL), INTENT(IN) :: intervals
  N_Detectors = N_Detectors + 1
  ALLOCATE(Detectors(N_Detectors)%p)
  Detectors(N_Detectors)%p = DARWIN_Xe_Init(LOGICAL(intervals))
  C_DARWIN_Xe_Init = N_Detectors
END FUNCTION


END MODULE
