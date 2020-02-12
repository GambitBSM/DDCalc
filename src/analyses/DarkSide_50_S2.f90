MODULE DarkSide_50_S2

!=======================================================================
! DarkSide-50 S2-ONLY ANALYSIS ROUTINES
! Based upon arXiv:1802.06994.  
!=======================================================================

USE DDTypes
USE DDDetectors

IMPLICIT NONE

CONTAINS


!-----------------------------------------------------------------------
! Initializes a DetectorStruct to the DarkSide-50 S2-only analysis.
! 
FUNCTION DarkSide_50_S2_Init() RESULT(D)

  IMPLICIT NONE
  TYPE(DetectorStruct) :: D
  INTEGER, PARAMETER :: NBINS = 3

  INTEGER, PARAMETER :: Nev_bin(Nbins) = (/680,630,4130/)
  REAL*8, PARAMETER :: Bg_bin(Nbins) = (/0.,0.,4080./)

! DarkSide specifies mass=46.4d0 and a total exposure of 6786 kg*days
! including an acceptance of 0.43. Since the acceptance is included separately
! we assume a run-time of 340 days to obtain a total exposure of 15800 kg*days  
  CALL SetDetector(D,mass=46.4d0,time=340.0d0,Nevents_bin=Nev_bin,      &
                   Backgr_bin=Bg_bin,Nelem=1,Zelem=(/18/),              &
                   Nbins=NBINS,E_file='DarkSide_50_S2/energies.dat',    &
                   eff_file=(/'DarkSide_50_S2/efficiencies.dat'/))
  
END FUNCTION  


! C++ interface wrapper
INTEGER(KIND=C_INT) FUNCTION C_DarkSide_50_S2_Init() &
 BIND(C,NAME='C_DDCalc_darkside_50_s2_init') 
  USE ISO_C_BINDING, only: C_BOOL, C_INT
  IMPLICIT NONE
  N_Detectors = N_Detectors + 1
  ALLOCATE(Detectors(N_Detectors)%p)
  Detectors(N_Detectors)%p = DarkSide_50_S2_Init()
  C_DarkSide_50_S2_Init = N_Detectors
END FUNCTION


END MODULE

