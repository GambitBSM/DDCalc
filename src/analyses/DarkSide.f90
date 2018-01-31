MODULE DarkSide

!=======================================================================
! DarkSide ANALYSIS ROUTINES
! Based upon arXiv:1707.08145.  
!=======================================================================

USE DDTypes
USE DDDetectors

IMPLICIT NONE

CONTAINS


!-----------------------------------------------------------------------
! Initializes a DetectorStruct to the DarkSide analysis.
! 
FUNCTION DarkSide_Init() RESULT(D)

  IMPLICIT NONE
  TYPE(DetectorStruct) :: D
! We approximate the energy resolution shown in Fig. 77 by sigma(E) = 0.05 * (E/keV) + 2
! We consider 10 bins of size 5 keV in the range from 30 to 80 keV
  INTEGER, PARAMETER :: NBINS = 10

  INTEGER, PARAMETER :: Nev_bin(Nbins) = (/0,0,0,0,0,0,0,0,0,0/)
! DarkSide claims to be virtually background-free. Here we assume 0.4 events, evenly distributed across all bins
  REAL*8, PARAMETER :: Bg_bin(Nbins) = (/0.04,0.04,0.04,0.04,0.04,0.04,0.04,0.04,0.04,0.04/)

! The efficiencies include the f_90 cut shown in Fig. 92
  CALL SetDetector(D,exposure=3.65d7,Nevents_bin=Nev_bin,               &
                   Backgr_bin=Bg_bin,Nelem=1,Zelem=(/18/),              &
                   Nbins=NBINS,E_file='data/DarkSide/energies.dat',     &
                   eff_file_all='data/DarkSide/efficiencies.dat')
  
END FUNCTION  


! C++ interface wrapper
INTEGER(KIND=C_INT) FUNCTION C_DarkSide_Init() &
 BIND(C,NAME='C_DDCalc_darkside_init') 
  USE ISO_C_BINDING, only: C_BOOL, C_INT
  IMPLICIT NONE
  N_Detectors = N_Detectors + 1
  ALLOCATE(Detectors(N_Detectors)%p)
  Detectors(N_Detectors)%p = DarkSide_Init()
  C_DarkSide_Init = N_Detectors
END FUNCTION


END MODULE
