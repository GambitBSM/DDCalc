MODULE DARWIN

!=======================================================================
! DARWIN ANALYSIS ROUTINES
! Based upon arXiv:1506.08309.  
!=======================================================================

USE DDTypes
USE DDDetectors

IMPLICIT NONE

CONTAINS


!-----------------------------------------------------------------------
! Initializes a DetectorStruct to the DARWIN analysis.
! 
FUNCTION DARWIN_Init() RESULT(D)

  IMPLICIT NONE
  TYPE(DetectorStruct) :: D
! Binning and resolution is taken from arXiv:1506.08309
  INTEGER, PARAMETER :: NBINS = 5

  INTEGER, PARAMETER :: Nev_bin(Nbins) = (/1,0,0,0,0/)
! We assume that the total background of 2.37 events is uniformly distributed
  REAL*8, PARAMETER :: Bg_bin(Nbins) = (/0.83,0.41,0.29,0.24,0.19/)

  CALL SetDetector(D,exposure=7.3d7,Nevents_bin=Nev_bin,                &
                   Backgr_bin=Bg_bin,Nelem=1,Zelem=(/54/),              &
                   Nbins=NBINS,E_file='data/DARWIN/energies.dat',       &
                   eff_file_all='data/DARWIN/efficiencies.dat')
  
END FUNCTION  


! C++ interface wrapper
INTEGER(KIND=C_INT) FUNCTION C_DARWIN_Init() &
 BIND(C,NAME='C_DDCalc_darwin_init') 
  USE ISO_C_BINDING, only: C_BOOL, C_INT
  IMPLICIT NONE
  N_Detectors = N_Detectors + 1
  ALLOCATE(Detectors(N_Detectors)%p)
  Detectors(N_Detectors)%p = DARWIN_Init()
  C_DARWIN_Init = N_Detectors
END FUNCTION


END MODULE
