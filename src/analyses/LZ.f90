MODULE LZ

!=======================================================================
! LZ ANALYSIS ROUTINES
! Based upon arXiv:1509.02910.  
!=======================================================================

USE DDTypes
USE DDDetectors

IMPLICIT NONE

CONTAINS


!-----------------------------------------------------------------------
! Initializes a DetectorStruct to the LZ analysis.
! 
FUNCTION LZ_Init() RESULT(D)

  IMPLICIT NONE
  TYPE(DetectorStruct) :: D
! Binning and resolution is taken from arXiv:1712.04793
  INTEGER, PARAMETER :: NBINS = 6

  INTEGER, PARAMETER :: Nev_bin(Nbins) = (/0,0,0,0,0,0/)
! We assume that the total background of 2.37 events is uniformly distributed
  REAL*8, PARAMETER :: Bg_bin(Nbins) = (/0.394,0.394,0.394,0.394,0.394,0.394/)

  CALL SetDetector(D,exposure=5.6d6,Nevents_bin=Nev_bin,                &
                   Backgr_bin=Bg_bin,Nelem=1,Zelem=(/54/),              &
                   Nbins=NBINS,E_file='data/LZ/energies.dat',           &
                   eff_file_all='data/LZ/efficiencies.dat')
  
END FUNCTION  


! C++ interface wrapper
INTEGER(KIND=C_INT) FUNCTION C_LZ_Init() &
 BIND(C,NAME='C_DDCalc_lz_init') 
  USE ISO_C_BINDING, only: C_BOOL, C_INT
  IMPLICIT NONE
  N_Detectors = N_Detectors + 1
  ALLOCATE(Detectors(N_Detectors)%p)
  Detectors(N_Detectors)%p = LZ_Init()
  C_LZ_Init = N_Detectors
END FUNCTION


END MODULE
