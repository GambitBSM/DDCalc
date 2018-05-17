MODULE CRESST_II

!=======================================================================
! CRESST_II ANALYSIS ROUTINES
! Based upon arXiv:1509.01515 and arXiv:1701.08157 .  
!=======================================================================

USE DDTypes
USE DDDetectors

IMPLICIT NONE

CONTAINS


!-----------------------------------------------------------------------
! Initializes a DetectorStruct to the CRESST_II analysis.
!
FUNCTION CRESST_II_Init() RESULT(D)

  IMPLICIT NONE
  TYPE(DetectorStruct) :: D
  INTEGER, PARAMETER :: NBINS = 10
  INTEGER, PARAMETER :: NELEM = 3

  INTEGER, PARAMETER :: Nev_bin(Nbins) = (/6, 17, 18, 21, 36, 57, 75, 97, 89, 84/)
  REAL*8, PARAMETER :: Bg_bin(Nbins) = (/0,0,0,0,0,0,0,0,0,0/)

  CALL SetDetector(D,mass=0.3d0,time=174.d0,Nevents_bin=Nev_bin,         &
                   Backgr_bin=Bg_bin,Nelem=NELEM,Nbins=NBINS,            &
                   Zelem=(/8,20,74/),stoich=(/4,1,1/),                   &
                   E_file='CRESST-II/energies.dat',                      &
                   eff_file=(/'CRESST-II/08.dat',                        &
                   'CRESST-II/20.dat',                                   &
                   'CRESST-II/74.dat'/))
  D%eff_file = '[CRESST_II]'
  
END FUNCTION


! C++ interface wrapper
INTEGER(KIND=C_INT) FUNCTION C_CRESST_II_Init() &
 BIND(C,NAME='C_DDCalc_cresst_ii_init') 
  USE ISO_C_BINDING, only: C_BOOL, C_INT
  IMPLICIT NONE
  N_Detectors = N_Detectors + 1
  ALLOCATE(Detectors(N_Detectors)%p)
  Detectors(N_Detectors)%p = CRESST_II_Init()
  C_CRESST_II_Init = N_Detectors
END FUNCTION


END MODULE
