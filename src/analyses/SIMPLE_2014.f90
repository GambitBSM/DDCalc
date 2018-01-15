MODULE SIMPLE_2014

!=======================================================================
! SIMPLE 2014 ANALYSIS ROUTINES
! Based upon the SIMPLE Phase II WIMP search: PRD 89, 072013 (2014)
! [1404.4309].  7+1 candidate events seen with 10.8+1.9 expected
! background events.
!=======================================================================

USE DDTypes
USE DDDetectors

IMPLICIT NONE

CONTAINS


!-----------------------------------------------------------------------
! Initializes a DetectorStruct to the SIMPLE 2014 analysis.
! 
FUNCTION SIMPLE_2014_Init() RESULT(D)

  IMPLICIT NONE
  TYPE(DetectorStruct) :: D
  REAL*8, ALLOCATABLE :: EFF_AllIso(:,:,:)
  INTEGER :: Kiso
  INTEGER :: K,Kmin,Kmax,Niso,Niso0
  INTEGER, ALLOCATABLE :: Ziso(:),Ziso0(:),Aiso(:),Aiso0(:)
  REAL*8, ALLOCATABLE :: fiso(:),fiso0(:),Miso(:),Miso0(:)
  
  ! Will build array of efficiencies below
  INTEGER :: NE
  REAL*8 :: Ethresh,Gamma
  REAL*8, ALLOCATABLE :: E(:),eff(:,:)
  
  ! Uses C2ClF5, but C threshold is ~ 100 keV, so we will drop
  ! the C components.  First get full set of isotopes to get mass
  ! fractions correct.  C is placed last in the stoichiometry to
  ! allow easier extraction of F & Cl parts.
  CALL CompoundIsotopeList(3,(/17,9,6/),(/1,5,2/),                      &
                           Niso0,Ziso0,Aiso0,fiso0,Miso0)
  Niso = COUNT(Ziso0 .NE. 6)
  ALLOCATE(Ziso(Niso),Aiso(Niso),fiso(Niso),Miso(Niso))
  Ziso = Ziso0(1:Niso)
  Aiso = Aiso0(1:Niso)
  fiso = fiso0(1:Niso)
  Miso = Miso0(1:Niso)
  
  ! Set efficiencies.
  ! Tabulation set at 100 per decade, up to 1000 keV.
  Ethresh = 8d0      ! F & Cl only, C is ~ 100 keV
  Gamma   = 4.2d0    ! 4.2 +/- 0.3
  Kmin = INT(100*LOG10(Ethresh))
  Kmax = NINT(100*LOG10(1d3))
  NE = (Kmax-Kmin)+1
  ALLOCATE(E(1:NE),eff(1:NE,0:0))
  DO K=Kmin,Kmax
    E(K-Kmin+1) = 10**(K/100d0)
  END DO
  eff(:,0) = 1d0 - EXP(-Gamma*(1d0-Ethresh/MAX(E,Ethresh)))
  
  ! One call for all settings.
  ! Most of these _must_ be there to ensure everything get initialized.
  CALL SetDetector(D,exposure=18.24d0,Nevents_tot=8,Backgr_tot=12.7d0,      &
                   Niso=Niso,Ziso=Ziso,Aiso=Aiso,fiso=fiso,                 &
                   NE=NE,E=E,Nbins=0,eff_all=eff)
  D%eff_file = '[SIMPLE 2014]'
  
END FUNCTION


! C++ interface wrapper
INTEGER(KIND=C_INT) FUNCTION C_SIMPLE_2014_Init() &
 BIND(C,NAME='C_DDCalc_simple_2014_init') 
  USE ISO_C_BINDING, only: C_BOOL, C_INT
  IMPLICIT NONE
  N_Detectors = N_Detectors + 1
  ALLOCATE(Detectors(N_Detectors)%p)
  Detectors(N_Detectors)%p = SIMPLE_2014_Init()
  C_SIMPLE_2014_Init = N_Detectors
END FUNCTION


END MODULE
