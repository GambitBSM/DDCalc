MODULE LZ_2022

!=======================================================================
! LZ_2022 ANALYSIS ROUTINES
!=======================================================================

USE DDTypes
USE DDDetectors

IMPLICIT NONE

CONTAINS


!-----------------------------------------------------------------------
! Initializes a DetectorStruct to the LZ_2022 analysis.
!
FUNCTION LZ_2022_Init() RESULT(D)

  IMPLICIT NONE
  TYPE(DetectorStruct) :: D
  INTEGER, PARAMETER :: NE = 200
  INTEGER, PARAMETER :: NBINS = 1
  REAL*8, PARAMETER :: EMIN = 2.5d-1


  CALL SetDetector(D,exposure=1.65d5,Nevents_tot=0,           &
                   Backgr_tot=0d0,Nelem=1,Zelem=(/54/),               &
                   Nbins=NBINS,E_file='LZ_2022/energies.dat',          &
                   eff_file_all='LZ_2022/efficiencies.dat')

END FUNCTION

! C++ interface wrapper
INTEGER(KIND=C_INT) FUNCTION C_LZ_2022_Init() &
 BIND(C,NAME='C_DDCalc_lz_2022_init') 
  USE ISO_C_BINDING, only: C_BOOL, C_INT
  IMPLICIT NONE
  N_Detectors = N_Detectors + 1
  ALLOCATE(Detectors(N_Detectors)%p)
  Detectors(N_Detectors)%p = LZ_2022_Init()
  C_LZ_2022_Init = N_Detectors
END FUNCTION


END MODULE
