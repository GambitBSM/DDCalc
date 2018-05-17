MODULE DarkSide_20k

!=======================================================================
! DarkSide_20k ANALYSIS ROUTINES
! Based upon arXiv:1707.08145.  
!=======================================================================

USE DDTypes
USE DDDetectors

IMPLICIT NONE

CONTAINS


!-----------------------------------------------------------------------
! Initializes a DetectorStruct to the DarkSide_20k analysis.
! 
FUNCTION DarkSide_20k_Init() RESULT(D)

  IMPLICIT NONE
  TYPE(DetectorStruct) :: D
! We approximate the energy resolution shown in Fig. 77 by sigma(E) = 0.05 * (E/keV) + 2
! We consider 10 bins of size 5 keV in the range from 30 to 80 keV
  INTEGER, PARAMETER :: NBINS = 10

! The dominant background is coherent neutrino scattering.
! We take the curves from https://arxiv.org/pdf/1307.5458.pdf and rescale them to
! 1.6 events in the ROI as stated in https://arxiv.org/pdf/1707.08145.pdf
! In addition, we add 0.04 events per bin from instrumental backgrounds
  REAL*8, PARAMETER :: Bg_bin(Nbins) = (/0.077, 0.117, 0.149, 0.156,    &
                       0.149, 0.138, 0.128, 0.119, 0.111, 0.102/)

  INTEGER, PARAMETER :: Nev_bin(Nbins) = (/0,0,0,0,0,0,0,0,0,0/)

! The efficiencies include the f_200 cut shown in Fig. 92
  CALL SetDetector(D,exposure=3.65d7,Nevents_bin=Nev_bin,               &
                   Backgr_bin=Bg_bin,Nelem=1,Zelem=(/18/),              &
                   Nbins=NBINS,E_file='DarkSide/energies.dat',          &
                   eff_file_all='DarkSide/efficiencies.dat')
  
END FUNCTION  


! C++ interface wrapper
INTEGER(KIND=C_INT) FUNCTION C_DarkSide_20k_Init() &
 BIND(C,NAME='C_DDCalc_darkside_20k_init') 
  USE ISO_C_BINDING, only: C_BOOL, C_INT
  IMPLICIT NONE
  N_Detectors = N_Detectors + 1
  ALLOCATE(Detectors(N_Detectors)%p)
  Detectors(N_Detectors)%p = DarkSide_20k_Init()
  C_DarkSide_20k_Init = N_Detectors
END FUNCTION


END MODULE
