MODULE CDMSlite

!=======================================================================
! CDMSlite ANALYSIS ROUTINES
! Based upon arXiv:1509.02448.  
!=======================================================================

USE DDTypes
USE DDDetectors
USE DDInput

IMPLICIT NONE

CONTAINS


!-----------------------------------------------------------------------
! Initializes a DetectorStruct to the CDMSlite analysis.
! 
FUNCTION CDMSlite_Init() RESULT(D)

  IMPLICIT NONE
  TYPE(DetectorStruct) :: D
  INTEGER, PARAMETER :: NBINS = 10

  INTEGER, PARAMETER :: Nev_bin(Nbins) = (/5, 4, 2, 3, 21, 6, 3, 5, 3, 3/)
  REAL*8, PARAMETER :: Bg_bin(Nbins) = (/0,0,0,0,0,0,0,0,0,0/)

  CALL SetDetector(D,exposure=70.1d0,Nevents_bin=Nev_bin,               &
                   Backgr_bin=Bg_bin,Nelem=1,Zelem=(/32/),              &
                   Nbins=NBINS,E_file='CDMSlite/energies.dat',          &
                   eff_file_all='CDMSlite/efficiencies.dat')
  D%eff_file = '[CDMSlite]'
  
END FUNCTION  


! C++ interface wrapper
INTEGER(KIND=C_INT) FUNCTION C_CDMSlite_Init() &
 BIND(C,NAME='C_DDCalc_cdmslite_init') 
  USE ISO_C_BINDING, only: C_BOOL, C_INT
  IMPLICIT NONE
  N_Detectors = N_Detectors + 1
  ALLOCATE(Detectors(N_Detectors)%p)
  Detectors(N_Detectors)%p = CDMSlite_Init()
  C_CDMSlite_Init = N_Detectors
END FUNCTION


END MODULE
