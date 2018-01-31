MODULE PICO_500

!=======================================================================
! PICO_500 ANALYSIS ROUTINES
! Based upon https://indico.cern.ch/event/606690/contributions/2623446/attachments/1497228/2330240/Fallows_2017_07_24__TAUP__PICO-40L_v1.2.pdf.  
!=======================================================================

USE DDTypes
USE DDDetectors

IMPLICIT NONE

CONTAINS


!-----------------------------------------------------------------------
! Initializes a DetectorStruct to the PICO_500 analysis.
! 
FUNCTION PICO_500_Init() RESULT(D)

  IMPLICIT NONE
  TYPE(DetectorStruct) :: D
! We consider the two modes with different threshold as separate bins.
  INTEGER, PARAMETER :: NBINS = 2

! The background in the lower bin is dominated by solar neutrinos, the background in the higher bin by neutrons.
  REAL*8, PARAMETER :: Bg_bin(Nbins) = (/3.0,0.9/)
! We assume the number of observed events closest to the number of expected events.
  INTEGER, PARAMETER :: Nev_bin(Nbins) = (/3,1/)

! The shorter run-time with low threshold is reflected by a reduced efficiency for the first bin.
! The total exposure is corrected for the mass fraction of fluorine.

  CALL SetDetector(D,exposure=1.03d5,Nevents_bin=Nev_bin,               &
                   Backgr_bin=Bg_bin,Nelem=1,Zelem=(/9/),               &
                   Nbins=NBINS,E_file='data/PICO-500/energies.dat',     &
                   eff_file_all='data/PICO-500/efficiencies.dat')
  
END FUNCTION  


! C++ interface wrapper
INTEGER(KIND=C_INT) FUNCTION C_PICO_500_Init() &
 BIND(C,NAME='C_DDCalc_pico_500_init') 
  USE ISO_C_BINDING, only: C_BOOL, C_INT
  IMPLICIT NONE
  N_Detectors = N_Detectors + 1
  ALLOCATE(Detectors(N_Detectors)%p)
  Detectors(N_Detectors)%p = PICO_500_Init()
  C_PICO_500_Init = N_Detectors
END FUNCTION


END MODULE
