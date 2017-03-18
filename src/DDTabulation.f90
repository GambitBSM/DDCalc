MODULE DDTabulation

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!  DDTabulation
!    Routines for determining tabulation points in DDCalc.
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

USE DDTypes
USE DDCommandLine

IMPLICIT NONE
PRIVATE

PUBLIC :: GetTabulationArgs, InitTabulation, TabulationValue, &
          TabulationInterval, TabulationBin

! Structure to contain fixed spacing linear or logarithmic tabulation
! parametrization.
TYPE, PUBLIC :: TabulationStruct
  ! Linear or logarithmic spacing?
  LOGICAL :: logarithmic = .FALSE.
  ! Number of intervals between tabulation points
  INTEGER :: N
  ! Range of tabulation
  REAL*8 :: xmin,xmax
  ! Logarithm of tabulation range values
  REAL*8 :: lnxmin,lnxmax
  ! Spacing between tabulation points:
  !   linear:       delta = x_{k+1} - x_k
  !   logarithmic:  delta = ln(x_{k+1}) - ln(x_k) = ln(x_{k+1}/x_k)
  REAL*8 :: delta
END TYPE

CONTAINS


! ----------------------------------------------------------------------
! Extracts parameters from arguments of the form:
!   --<akey>=<min>,<max>,<N>,<use_log>
! Does not change the existing value if the parameter is not given
! or cannot be parsed.  The integer argument N can be given as a
! floating-point number and the logical argument use_log can be
! as 'T'/'F' or numerically (zero is false, all other values true).
! 
SUBROUTINE GetTabulationArgs(akey,xmin,xmax,N,use_log)
  IMPLICIT NONE
  CHARACTER*(*), INTENT(IN) :: akey
  REAL*8, INTENT(INOUT) :: xmin,xmax
  INTEGER, INTENT(INOUT) :: N
  LOGICAL, INTENT(INOUT) :: use_log
  LOGICAL :: status,Ltmp
  INTEGER :: Nval,ios,Itmp
  REAL*8 :: Rtmp
  ! Older compiler compatibility
  INTEGER, PARAMETER :: NCHAR = 32
  CHARACTER(LEN=NCHAR), DIMENSION(:), ALLOCATABLE :: aval
  ! ...but this would be better (needs gfortran 4.6+)
  !CHARACTER(LEN=:), DIMENSION(:), ALLOCATABLE :: aval
  
  IF (.NOT. GetLongArgStrings(akey,NCHAR,aval,Nval)) RETURN
  
  IF (Nval .GE. 1) THEN
    READ(UNIT=aval(1),FMT=*,IOSTAT=ios) Rtmp
    IF (ios .EQ. 0) xmin = Rtmp
  END IF
  
  IF (Nval .GE. 2) THEN
    READ(UNIT=aval(2),FMT=*,IOSTAT=ios) Rtmp
    IF (ios .EQ. 0) xmax = Rtmp
  END IF
  
  IF (Nval .GE. 3) THEN
    READ(UNIT=aval(3),FMT=*,IOSTAT=ios) Rtmp  ! read as real
    IF (ios .EQ. 0) N = Rtmp
  END IF
  
  IF (Nval .GE. 4) THEN
    READ(UNIT=aval(4),FMT=*,IOSTAT=ios) Ltmp
    IF (ios .EQ. 0) THEN
      use_log = Ltmp
    ELSE
      READ(UNIT=aval(4),FMT=*,IOSTAT=ios) Rtmp
      IF (ios .EQ. 0) use_log = .NOT. (Rtmp .EQ. 0)
    END IF
  END IF
  
END SUBROUTINE


! ----------------------------------------------------------------------
! Initializes the given TabulationStruct to use a tabulation specified
! by the given parameters.  Some sanity checks are performed.
! 
! Output argument:
!   TS              TabulationStruct to be filled with tabulation
!                   parameterization data.
! Input arguments:
!   xmin,xmax       Tabulation range.  Might be adjusted by up to 1/2
!                   bin size to ensure even sized bins in logarithmic
!                   case with N < 0.
!   N               Number of intervals between tabulation points.
!                   A negative number indicates intervals per decade.
!   use_log         If logarithmic rather than linear spacing is to
!                   be used.
! 
SUBROUTINE InitTabulation(TS,xmin,xmax,N,use_log)
  IMPLICIT NONE
  TYPE(TabulationStruct), INTENT(OUT) :: TS
  REAL*8, INTENT(IN) :: xmin,xmax
  INTEGER, INTENT(IN) :: N
  LOGICAL, INTENT(IN) :: use_log
  
  TS%logarithmic = use_log
  
  IF (use_log) THEN
    IF (xmin .LE. 0d0) THEN
      IF (xmax .LE. 0d0) THEN
        TS%xmin =   1d0
        TS%xmax = 100d0
      ELSE
        TS%xmin = xmax/100
        TS%xmax = xmax
      END IF
    ELSE
      IF (xmax .LE. 0d0) THEN
        TS%xmin = xmin
        TS%xmax = 100*xmin
      ELSE
        TS%xmin = xmin
        TS%xmax = xmax
      END IF
    END IF
    IF (N .EQ. 0) THEN
      TS%N = -10
    ELSE
      TS%N = N
    END IF
    TS%lnxmin = LOG(TS%xmin)
    TS%lnxmax = LOG(TS%xmax)
  ELSE
    TS%xmin = MIN(xmin,xmax)
    TS%xmax = MAX(xmin,xmax)
    IF (N .EQ. 0) THEN
      TS%N = 100
    ELSE
      TS%N = ABS(N)
    END IF
    TS%lnxmin = 0d0
    TS%lnxmax = 0d0
  END IF
  
  ! Special case (single point)
  IF (TS%xmax .EQ. TS%xmin) THEN
    TS%N     = 0
    TS%delta = 0d0
    RETURN
  END IF
  
  ! Logarithmic case
  IF (use_log) THEN
    IF (TS%N .LT. 0) THEN
      TS%delta = LOG(10d0) / ABS(TS%N)
      TS%N     = NINT(ABS(TS%N) * LOG10(TS%xmax/TS%xmin))
      TS%xmax  = TS%xmin*EXP(TS%N*TS%delta)
      TS%lnxmax = LOG(TS%xmax)
    ELSE
      TS%delta = (TS%lnxmax - TS%lnxmin) / TS%N
    END IF
  ! Linear case
  ELSE
    TS%delta = (TS%xmax - TS%xmin) / TS%N
  END IF
  
END SUBROUTINE


! ----------------------------------------------------------------------
! Returns the Kth tabulation value, with K=0 (TS%N) corresponding to
! the minimum (maximum) of the tabulation range.
! 
PURE FUNCTION TabulationValue(TS,K) RESULT (x)
  IMPLICIT NONE
  REAL*8 :: x
  TYPE(TabulationStruct), INTENT(IN) :: TS
  INTEGER, INTENT(IN) :: K
  
  IF (TS%logarithmic) THEN
    x = EXP(TS%lnxmin + K*TS%delta)
  ELSE
    x = TS%xmin + K*TS%delta
  END IF
  
END FUNCTION


! ----------------------------------------------------------------------
! Returns the tabulation index K such that x_K <= x < x_{K+1}, or -1 if
! x does not fall within the tabulation range.  Note tabulation points
! are x_K with K=0,...,TS%N.
! 
PURE FUNCTION TabulationInterval(TS,x) RESULT (K)
  IMPLICIT NONE
  INTEGER :: K
  TYPE(TabulationStruct), INTENT(IN) :: TS
  REAL*8, INTENT(IN) :: x
  
  IF (TS%N .EQ. 0) THEN
    IF (x .EQ. TS%xmin) THEN
      K = 0
    ELSE
      K = -1
    END IF
    RETURN
  END IF
  
  IF (x .LT. TS%xmin) THEN
    K = -1
  ELSE IF (TS%logarithmic) THEN
    K = INT((LOG(x)-TS%lnxmin)/TS%delta)
  ELSE
    K = INT((x-TS%xmin)/TS%delta)
  END IF
  IF ((K .LT. 0) .OR. (K .GT. TS%N-1)) K = -1
  
END FUNCTION


! ----------------------------------------------------------------------
! Returns the tabulation index K such that x falls within a bin centered
! on x_K and with a width equal to the tabulation spacing.  A value of
! -1 is returned if x falls outside the tabulation bins.  For
! logarithmic tabulation, the bin center and size are taken in
! logarithmic space.  Note tabulation points are x_K with K=0,...,TS%N.
! 
PURE FUNCTION TabulationBin(TS,x) RESULT (K)
  IMPLICIT NONE
  INTEGER :: K
  TYPE(TabulationStruct), INTENT(IN) :: TS
  REAL*8, INTENT(IN) :: x
  
  IF (TS%N .EQ. 0) THEN
    IF (x .EQ. TS%xmin) THEN
      K = 0
    ELSE
      K = -1
    END IF
    RETURN
  END IF
  
  IF (TS%logarithmic) THEN
    K = NINT((LOG(x)-TS%lnxmin)/TS%delta)
  ELSE
    K = NINT((x-TS%xmin)/TS%delta)
  END IF
  IF ((K .LT. 0) .OR. (K .GT. TS%N)) K = -1
  
END FUNCTION


END MODULE
