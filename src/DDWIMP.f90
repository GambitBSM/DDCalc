MODULE DDWIMP

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! DDWIMP
!    Routines for initializing, modifying, or viewing WIMP properties,
!    mainly the mass and couplings.
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

USE DDTypes
USE DDUtils
USE DDCouplings

IMPLICIT NONE
PRIVATE

! WIMP mass & couplings routines
PUBLIC :: DDCalc_GetWIMP,DDCalc_SetWIMP,     &
          DDCalc_InitWIMP,                   &
          C_DDCalc_InitWIMP
INTERFACE DDCalc_GetWIMP
  MODULE PROCEDURE GetWIMP
END INTERFACE
INTERFACE DDCalc_SetWIMP
  MODULE PROCEDURE SetWIMP
END INTERFACE
INTERFACE DDCalc_InitWIMP
  MODULE PROCEDURE InitWIMP
END INTERFACE

CONTAINS


! ----------------------------------------------------------------------
! Get various WIMP quantities.
! 
! Optional output arguments:
!   m           WIMP mass [GeV].
!   DMtype      WIMP type (see definition of WIMP stucture)
!   params      WIMP parameters
!   Nparams     Number of WIMP parameters
SUBROUTINE GetWIMP(WIMP, m, DMtype, params, Nparams)

  IMPLICIT NONE
  TYPE(WIMPStruct), INTENT(IN) :: WIMP
  REAL*8, INTENT(OUT), OPTIONAL :: m
  CHARACTER(LEN=24), INTENT(OUT), OPTIONAL :: DMtype
  REAL*8, ALLOCATABLE, INTENT(OUT), OPTIONAL :: params(:)
  INTEGER, INTENT(OUT), OPTIONAL :: Nparams

  IF (PRESENT(m))       m       = WIMP%m
  IF (PRESENT(DMtype))  DMtype  = WIMP%DMtype
  IF (PRESENT(Nparams)) Nparams  = WIMP%Nparams
  IF (PRESENT(params)) THEN
    ALLOCATE(params(WIMP%Nparams))
    params = WIMP%params
  END IF

END SUBROUTINE


! ----------------------------------------------------------------------
! Set various WIMP quantities.
! 
! Optional input arguments:
!   m           WIMP mass [GeV].
!   DMtype      WIMP type (see definition of WIMP stucture)
!   params      WIMP parameters (must have the correct length)
! Note that both DMtype and params must be specified at the same time
SUBROUTINE SetWIMP(WIMP, m, DMtype, params)

  IMPLICIT NONE
  TYPE(WIMPStruct), INTENT(OUT) :: WIMP
  REAL*8, INTENT(IN), OPTIONAL :: m
  CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: DMtype
  REAL*8, INTENT(IN), OPTIONAL :: params(:)
  LOGICAL :: updated  

  IF (PRESENT(m))       WIMP%m       = m
  IF (PRESENT(DMtype) .AND. PRESENT(params)) THEN
    updated = .FALSE.
    IF ( DMtype .EQ. 'SIonly' ) THEN
      WIMP%DMtype  = DMtype
      WIMP%Nparams = 2
      updated = .TRUE.
    END IF
    IF ( DMtype .EQ. 'SDonly' ) THEN
      WIMP%DMtype  = DMtype
      WIMP%Nparams = 2
      updated = .TRUE.
    END IF
    IF ( DMtype .EQ. 'SISD' ) THEN
      WIMP%DMtype  = DMtype
      WIMP%Nparams = 4
      updated = .TRUE.
    END IF
    IF ( DMtype .EQ. 'SILR' ) THEN
      WIMP%DMtype  = DMtype
      WIMP%Nparams = 3
      updated = .TRUE.
    END IF
    IF ( DMtype .EQ. 'NREffectiveTheory' ) THEN
      WIMP%DMtype  = DMtype
      WIMP%Nparams = 45
      updated = .TRUE.
    END IF
    IF ( DMtype .EQ. 'NREFT_CPT' ) THEN
      WIMP%DMtype  = DMtype
      WIMP%Nparams = 51
      updated = .TRUE.
    END IF
    IF ( updated ) THEN
      ALLOCATE(WIMP%params(WIMP%Nparams))
      WIMP%params = params(1:WIMP%Nparams)
    ELSE
      WRITE(0,*) 'WARNING: WIMP type not recognized by SetWIMP.'
    END IF
  ELSE
    IF (PRESENT(DMtype) .OR. PRESENT(params)) THEN
      WRITE(0,*) 'WARNING: Incomplete information provided to SetWIMP.'
    END IF
  END IF

END SUBROUTINE

!-----------------------------------------------------------------------
! Initializes WIMP.
! Simply sets some default values for the WIMP parameters.
! 
FUNCTION InitWIMP() RESULT(WIMP)

  IMPLICIT NONE
  TYPE(WIMPStruct) :: WIMP

  CALL SetWIMP(WIMP,m=100d0,DMtype='SIonly',params=[1d-9,1d-9])
  
END FUNCTION

!-----------------------------------------------------------------------
! C/C++ wrapper for InitWIMP.
!
INTEGER(KIND=C_INT) FUNCTION C_DDCalc_InitWIMP() &
 BIND(C,NAME='C_DDWIMP_ddcalc_initwimp') 
  USE ISO_C_BINDING, only: C_INT
  IMPLICIT NONE
  N_WIMPs = N_WIMPs + 1
  IF (N_WIMPs .GT. Max_WIMPs) stop 'DDCalc: Max_WIMPs exceeded.&
   Please run FreeWIMPs or modify Max_WIMPs in DDTypes.f90.'
  ALLOCATE(WIMPs(N_WIMPs)%p)
  WIMPs(N_WIMPs)%p = InitWIMP()
  C_DDCalc_InitWIMP = N_WIMPs
END FUNCTION

END MODULE
