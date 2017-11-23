MODULE CRESST_II

!=======================================================================
! CRESST_III ANALYSIS ROUTINES
! Based upon https://arxiv.org/pdf/1510.07754v3.pdf .  
!=======================================================================

USE DDTypes
USE DDDetectors

IMPLICIT NONE

CONTAINS


!-----------------------------------------------------------------------
! Initializes a DetectorStruct to the PICO_60 analysis.
!
FUNCTION CRESST_II_Init(intervals) RESULT(D)

  IMPLICIT NONE
  TYPE(DetectorStruct) :: D
  LOGICAL, INTENT(IN) :: intervals
  INTEGER, PARAMETER :: NBINS = 10
  INTEGER, PARAMETER :: NELEM = 3

  INTEGER, PARAMETER :: Nev_bin(Nbins) = (/6,11,13,11,9,12,5,19,8,19/)
  REAL*8, PARAMETER :: Bg_bin(Nbins) = (/0,0,0,0,0,0,0,0,0,0/)

  CALL SetDetector(D,mass=0.3d0,time=174.d0,Nevents_bin=Nev_bin,         &
                   Backgr_bin=Bg_bin,Nelem=NELEM,Nbins=NBINS,            &
                   Zelem=(/8,20,74/),stoich=(/4,1,1/),                   &
                   E_file='data/CRESST-II/energies.dat',                 &
                   eff_file=(/'data/CRESST-II/08.dat',                   &
                   'data/CRESST-II/20.dat',                              &
                   'data/CRESST-II/74.dat'/),                            &
                   intervals=intervals)
  D%eff_file = '[CRESST_II]'
  
END FUNCTION


! C++ interface wrapper
INTEGER(KIND=C_INT) FUNCTION C_CRESST_II_Init(intervals) &
 BIND(C,NAME='C_DDCalc_cresst_ii_init') 
  USE ISO_C_BINDING, only: C_BOOL, C_INT
  IMPLICIT NONE
  LOGICAL(KIND=C_BOOL), INTENT(IN) :: intervals
  N_Detectors = N_Detectors + 1
  ALLOCATE(Detectors(N_Detectors)%p)
  Detectors(N_Detectors)%p = CRESST_II_Init(LOGICAL(intervals))
  C_CRESST_II_Init = N_Detectors
END FUNCTION


END MODULE
