MODULE PandaX_4T

!=======================================================================
! PandaX_4T ANALYSIS ROUTINES
!=======================================================================

USE DDTypes
USE DDDetectors

IMPLICIT NONE

CONTAINS


!-----------------------------------------------------------------------
! Initializes a DetectorStruct to the PandaX_4T analysis.
!
FUNCTION PandaX_4T_Init() RESULT(D)

  IMPLICIT NONE
  TYPE(DetectorStruct) :: D
  INTEGER, PARAMETER :: NE = 201
  INTEGER, PARAMETER :: NBINS = 1
  REAL*8, PARAMETER :: EMIN = 1.0d0


  CALL SetDetector(D,exposure=1.15d5,Nevents_tot=6,  &
                   Backgr_tot=9.8d0,Nelem=1,Zelem=(/54/),               &
                   Nbins=NBINS,E_file='PandaX_4T/energies.dat',          &
                   eff_file_all='PandaX_4T/efficiencies.dat')

END FUNCTION

! C++ interface wrapper
INTEGER(KIND=C_INT) FUNCTION C_PandaX_4T_Init() &
 BIND(C,NAME='C_DDCalc_pandax_4t_init') 
  USE ISO_C_BINDING, only: C_BOOL, C_INT
  IMPLICIT NONE
  N_Detectors = N_Detectors + 1
  ALLOCATE(Detectors(N_Detectors)%p)
  Detectors(N_Detectors)%p = PandaX_4T_Init()
  C_PandaX_4T_Init = N_Detectors
END FUNCTION


END MODULE
