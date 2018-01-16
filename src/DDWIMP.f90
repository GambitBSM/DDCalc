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
          ParseWIMPInput, C_DDCalc_InitWIMP
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
    IF ( DMtype .EQ. 'HiggsPortal' ) THEN
      WIMP%DMtype  = DMtype
      WIMP%Nparams = 4
      updated = .TRUE.
    END IF
    IF ( DMtype .EQ. 'NREffectiveTheory' ) THEN
      WIMP%DMtype  = DMtype
      WIMP%Nparams = 37
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



!-----------------------------------------------------------------------
! Routine to read a line containing WIMP parameters from standard input
! (if line not given explicitly), parse it, and set WIMP parameters
! accordingly.  Returns true if a non-empty line was found and was
! parsable (false if no line found (EOF), line was empty, or line
! could not be parsed).
! 
! The following forms are allowed:
!   m
!   m sigmaSI
!   m sigmaSI sigmaSD
!   m sigmaSI sigmanSD sigmanSD
!   m sigmapSI sigmanSI sigmanSD sigmanSD
!   (formerly:  m GpSI GnSI GpSD GnSD)
! In the first case, all couplings are set to 1 pb, while in the second
! case, the SD couplings are set to zero.
! 
! Optional input argument:
!   line        String containing WIMP parameters.  If not given,
!               a line is taken from standard output.
! 
FUNCTION ParseWIMPInput(WIMP,line) RESULT(status)
  IMPLICIT NONE
  TYPE(WIMPStruct), INTENT(INOUT) :: WIMP
  LOGICAL :: status
  CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: line
  CHARACTER(LEN=1024) :: line0
  INTEGER :: ios,Np
  REAL*8 :: params(6)
  
  status = .FALSE.
  
  ! Use given line or read from standard input
  IF (PRESENT(line)) THEN
    line0 = line
  ELSE
    READ(*,'(A)',IOSTAT=ios) line0
    IF (ios .NE. 0) RETURN
  END IF
  
  ! Check if empty string
  IF (TRIM(line0) .EQ. '') RETURN
  
  ! Determine number of parameters and read them in
  Np = MIN(NumberOfFields(line0),6)
  IF (Np .LE. 0) RETURN
  READ(line0,*,IOSTAT=ios) params(1:Np)
  IF (ios .NE. 0) RETURN
  
  ! Now parse parameters
  status = ParseWIMPParameters(Np,params,WIMP)
  
END FUNCTION


!-----------------------------------------------------------------------
! Routine to take an array containing WIMP parameters and set the
! internal WIMP parameters to that.  The meaning of the given
! parameters depends upon the number of parameters and is the same
! as expected for the commandline.  Returns false if an invalid
! number of parameters or invalid value is found.
! 
! The following array lengths and parameters are allowed:
!   N=1:  m
!   N=2:  m sigmaSI
!   N=3:  m sigmaSI sigmaSD
!   N=4:  m sigmaSI sigmanSD sigmanSD
!   N=5:  m sigmapSI sigmanSI sigmanSD sigmanSD
!         (formerly:  m GpSI GnSI GpSD GnSD)
! In the first case, all couplings are set to 1 pb, while in the second
! case, the SD couplings are set to zero.
! 
! Input arguments:
!   N           Number of parameters
!   p           Array of parameters of size [1:N]
! 
FUNCTION ParseWIMPParameters(N,p,WIMP) RESULT(status)

  IMPLICIT NONE
  TYPE(WIMPStruct), INTENT(INOUT) :: WIMP
  LOGICAL :: status
  INTEGER, INTENT(IN) :: N
  REAL*8, INTENT(IN) :: p(N)
  REAL*8 :: m, fp, fn, ap, an
  
  ! Check for bad cases (e.g. non-positive mass)
  status = .FALSE.
  IF (N .LT. 1) RETURN
  IF (p(1) .LE. 0d0) RETURN
  
  status = .TRUE.
  m=p(1)
  
  ! Meaning of parameters depends on number of parameters.

  SELECT CASE (N)

  CASE (1)
    fp = 1d-9
    fn = 1d-9
    ap = 0d0
    an = 0d0
  CASE (2)
    ! Form: m sigmaSI
    ! Set SD couplings to zero
    fp = SigmapSItoFp(m,p(2))
    fn = SigmanSItoFn(m,p(2))
    ap = 0d0
    an = 0d0
  CASE (3)
    ! Form: m sigmaSI sigmaSD
    fp = SigmapSItoFp(m,p(2))
    fn = SigmanSItoFn(m,p(2))
    ap = SigmapSDtoAp(m,p(3))
    an = SigmanSDtoAn(m,p(3))
  CASE (4)
    fp = SigmapSItoFp(m,p(2))
    fn = SigmanSItoFn(m,p(2))
    ap = SigmapSDtoAp(m,p(3))
    an = SigmanSDtoAn(m,p(4))
  CASE (5)
    fp = SigmapSItoFp(m,p(2))
    fn = SigmanSItoFn(m,p(3))
    ap = SigmapSDtoAp(m,p(4))
    an = SigmanSDtoAn(m,p(5))
  CASE (6:)
    status = .FALSE.
  END SELECT

  if ( status ) CALL SetWIMP(WIMP,m=m, DMtype='SISD', params=[fp,fp,ap,an])
  
END FUNCTION


END MODULE
