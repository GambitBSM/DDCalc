MODULE CRESST_III

!=======================================================================
! CRESST_III ANALYSIS ROUTINES
! Based upon arXiv:1509.01515 and arXiv:1701.08157 .
!=======================================================================

USE DDTypes
USE DDDetectors

IMPLICIT NONE

CONTAINS


!-----------------------------------------------------------------------
! Initializes a DetectorStruct to the CRESST_III analysis.
!
FUNCTION CRESST_III_Init() RESULT(D)

  IMPLICIT NONE
  TYPE(DetectorStruct) :: D
  INTEGER, PARAMETER :: NBINS = 10
  INTEGER, PARAMETER :: NELEM = 3

  INTEGER, PARAMETER :: Nev_bin(Nbins) = (/355, 219, 103, 57, 21, 7, 9, 11, 18, 74/)
  REAL*8, PARAMETER :: Bg_bin(Nbins) = (/0,0,0,0,0,0,0,0,0,0/)

  CALL SetDetector(D,mass=0.024d0,time=237.d0,Nevents_bin=Nev_bin,         &
                   Backgr_bin=Bg_bin,Nelem=NELEM,Nbins=NBINS,            &
                   Zelem=(/8,20,74/),stoich=(/4,1,1/),                   &
                   E_file='CRESST-III/energies.dat',                      &
                   eff_file=(/'CRESST-III/08.dat',                        &
                   'CRESST-III/20.dat',                                   &
                   'CRESST-III/74.dat'/))
  D%eff_file = '[CRESST_III]'

END FUNCTION


! C++ interface wrapper
INTEGER(KIND=C_INT) FUNCTION C_CRESST_III_Init() &
 BIND(C,NAME='C_DDCalc_cresst_iii_init')
  USE ISO_C_BINDING, only: C_BOOL, C_INT
  IMPLICIT NONE
  N_Detectors = N_Detectors + 1
  ALLOCATE(Detectors(N_Detectors)%p)
  Detectors(N_Detectors)%p = CRESST_III_Init()
  C_CRESST_III_Init = N_Detectors
END FUNCTION


END MODULE
