MODULE DDUtils

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! DDUtils
!    DDCalc utility functions
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

USE DDTypes

IMPLICIT NONE
PRIVATE

PUBLIC :: BSearch, NumberOfFields, CoerceNumber, CoerceExponent
PUBLIC :: C_DDCalc_FreeWIMPs, C_DDCalc_FreeHalos, C_DDCalc_FreeDetectors
PUBLIC :: C_DDCalc_FreeAll

CONTAINS


!=======================================================================
! C++ BOOKKEEPING
!=======================================================================

SUBROUTINE C_DDCalc_FreeWIMPs() BIND(C,NAME='C_DDUtils_ddcalc_freewimps')
  IMPLICIT NONE
  INTEGER I
  IF (N_WIMPs .EQ. 0) RETURN
  DO I = 1, N_WIMPs
    DEALLOCATE(WIMPs(I)%p)
  ENDDO
  N_WIMPs = 0
END SUBROUTINE

SUBROUTINE C_DDCalc_FreeHalos() BIND(C,NAME='C_DDUtils_ddcalc_freehalos')
  IMPLICIT NONE
  INTEGER I
  IF (N_Halos .EQ. 0) RETURN
  DO I = 1, N_Halos
    DEALLOCATE(Halos(I)%p)
  ENDDO
  N_Halos = 0
END SUBROUTINE

SUBROUTINE C_DDCalc_FreeDetectors() BIND(C,NAME='C_DDUtils_ddcalc_freedetectors')
  IMPLICIT NONE
  INTEGER I
  IF (N_Detectors .EQ. 0) RETURN
  DO I = 1, N_Detectors
    DEALLOCATE(Detectors(I)%p)
  ENDDO
  N_Detectors = 0
END SUBROUTINE

SUBROUTINE C_DDCalc_FreeAll() BIND(C,NAME='C_DDUtils_ddcalc_freeall')
  IMPLICIT NONE
  CALL C_DDCalc_FreeWIMPs()
  CALL C_DDCalc_FreeHalos()
  CALL C_DDCalc_FreeDetectors()
END SUBROUTINE


!=======================================================================
! NUMBER FORMATTING
!=======================================================================

! ----------------------------------------------------------------------
! Converts the given number to a string of width w, with nl characters
! to the left of the decimal point, np digits of precision, and ne
! digits in the exponential (set ne=0 for floating point only).  The
! number will be printed in floating point format if possible,
! otherwise exponential notation is used.  Numbers are coerced into
! a representable range.
! 
! Requirements:
!   nl + np + 2 + ne <= w  [ne > 0]
! 
PURE FUNCTION NumberToString(x,w,nl,np,ne) RESULT(s)
  IMPLICIT NONE
  REAL*8, INTENT(IN) :: x
  INTEGER, INTENT(IN) :: w
  INTEGER, INTENT(IN), OPTIONAL :: nl,np,ne
  CHARACTER*(w) :: s
  LOGICAL :: no_precision,use_fixed
  INTEGER :: nl0,nr0,np0,ne0,xsgn,exponent,expmax,decloc
  REAL*8 :: x0,eps,mantissa
  CHARACTER*64 :: XFMT,MFMT,EFMT,t
  
  ! Determine number of characters
  nl0 = 2
  np0 = 0
  ne0 = 0
  IF (PRESENT(nl)) nl0 = nl
  IF (PRESENT(np)) np0 = np
  IF (PRESENT(ne)) ne0 = ne
  nl0 = MAX(nl0,2)
  nr0 = w - nl0 - 1
  ne0 = MIN(MAX(ne0,0),9)
  no_precision = .FALSE.
  IF (Np0 .LE. 0) THEN
    no_precision = .TRUE.
    IF (Ne0 .GT. 0) THEN
      Np0 = Nr0 - 1 - Ne
    ELSE
      Np0 = Nr0 + 1
    END IF
  END IF
  
  ! Check for bad cases
  IF ((nr0 .LT. 0) .OR. (Np0 .LE. 0)) THEN
    s = REPEAT('*',w)
    RETURN
  END IF
  
  ! Use modifiable variable (x is fixed)
  x0 = x
  
  IF (x0 .GE. 0) THEN
    xsgn = +1
  ELSE
    xsgn = -1
  END IF
  eps = 10d0**(-Np0)
  
  ! Exponential notation components
  IF (x0 .NE. 0d0) THEN
    exponent = FLOOR(LOG10(ABS(x0*(1+0.51d0*eps))))
  ELSE
    exponent = 0d0
  END IF
  mantissa = x * 10d0**(-exponent)
  
  ! Coerce number into representable range
  IF (Ne0 .GT. 0) THEN
    expmax = 10**Ne0
    IF (ABS(exponent) .GE. expmax) THEN
      IF (exponent .GT. 0) THEN
        mantissa = xsgn * 10 * (1-eps)
        exponent = expmax - 1
      ELSE IF (Np0 - 1 + exponent .GT. -expmax) THEN
        mantissa = mantissa * 10d0**(exponent+expmax-1)
        exponent = -expmax + 1
      ELSE
        mantissa = 0d0
        exponent = 0
        x0       = 0d0
      END IF
    END IF
  ELSE
    IF ((x0 .GT. 0) .AND. (exponent + 1 .GT. Nl0)) THEN
      IF (no_precision) eps = MAX(10d0**(-w+2),EPSILON(eps))
      mantissa = 10 * (1-eps)
      exponent = Nl0 - 1
      x0       = mantissa * 10d0**exponent
    ELSE IF ((x0 .LT. 0) .AND. (exponent + 2 .GT. Nl0)) THEN
      IF (no_precision) eps = MAX(10d0**(-w+3),EPSILON(eps))
      mantissa = -10 * (1-eps)
      exponent = Nl0 - 2
      x0       = mantissa * 10d0**exponent
    END IF
  END IF
  
  ! Determine if fixed notation should be used
  use_fixed = .TRUE.
  IF (Ne0 .GT. 0) THEN
    IF (Np0 - 1 - exponent .GT. Nr0) THEN
      use_fixed = .FALSE.
    ELSE IF ((x0 .GT. 0) .AND. (exponent + 1 .GT. Nl0)) THEN
      use_fixed = .FALSE.
    ELSE IF ((x0 .LT. 0) .AND. (exponent + 2 .GT. Nl0)) THEN
      use_fixed = .FALSE.
    END IF
  END IF
  
  ! Construct string
  ! Fixed format
  IF (use_fixed) THEN
    IF (no_precision) THEN
      WRITE(XFMT,'(A,I2,A,I2,A)') '(F',w,'.',Nr0,')'
    ELSE
      WRITE(XFMT,'(A,I2,A,I2,A)') '(F',w,'.',MAX(Np0-exponent-1,0),')'
    END IF
    WRITE(t,XFMT) x0
    t = ADJUSTL(t)
    decloc = INDEX(t,'.')
    IF (decloc .EQ. 0) THEN
      decloc = LEN_TRIM(t) + 1
      !t = TRIM(t) // '.'
    END IF
    IF (decloc .LT. nl0 + 1) THEN
      s = REPEAT(' ',nl0 + 1 - decloc) // t
    ELSE IF (decloc .GT. nl0 + 1) THEN
      ! This should not happen....
      s = t(decloc-nl0:)
    ELSE
      s = t
    END IF
  ! Exponential format
  ELSE
    ! Mantissa part
    WRITE(MFMT,'(A,I2,A,I2,A)') '(F',Np0+2,'.',Np0-1,')'
    WRITE(t,MFMT) mantissa
    t = ADJUSTL(t)
    decloc = INDEX(t,'.')
    IF (decloc .EQ. 0) THEN
      decloc = LEN_TRIM(t) + 1
      !t = TRIM(t) // '.'
    END IF
    IF (decloc .LT. nl0 + 1) THEN
      s = REPEAT(' ',nl0 + 1 - decloc) // t
    ELSE IF (decloc .GT. nl0 + 1) THEN
      ! This should not happen....
      s = t(decloc-nl0:)
    ELSE
      s = t
    END IF
    ! Exponential part
    WRITE(EFMT,'(A,I1,A)') '(I',Ne0,')'
    WRITE(t,EFMT) ABS(exponent)
    t = ADJUSTL(t)
    IF (LEN_TRIM(t) .LT. Ne0) t = REPEAT('0',Ne0-LEN_TRIM(t)) // TRIM(t)
    IF (exponent .GE. 0) THEN
      t = 'E+' // TRIM(t)
    ELSE
      t = 'E-' // TRIM(t)
    END IF
    ! Combine parts
    s = TRIM(s) // t
  END IF
  
END FUNCTION


! ----------------------------------------------------------------------
! Returns the given number coerced into a range where it can be printed
! in at most w characters, optionally with a fixed number of decimal
! digits d (not using scientific notation).  Large magnitude numbers
! are reduced and very small magnitude numbers are set to zero.
! 
ELEMENTAL FUNCTION CoerceNumber(x,w,d) RESULT(y)
  IMPLICIT NONE
  REAL*8 :: y
  REAL*8, INTENT(IN) :: x
  INTEGER, INTENT(IN) :: w
  INTEGER, INTENT(IN), OPTIONAL :: d
  REAL*8 :: xmin,xmax
  y = x
  IF (PRESENT(d)) THEN
    xmin = (0.5d0+EPSILON(1d0)) * 10d0**(-d)
    IF (ABS(x) .LT. xmin) THEN
      y = 0d0
    ELSE IF (x .GE. 0) THEN
      xmax = 10d0**(w-d-1) - 2*xmin
      IF (x .GT. xmax) y = xmax
    ELSE
      xmax = 10d0**(w-d-2) - 2*xmin
      IF (ABS(x) .GT. xmax) y = -xmax
    END IF
  ELSE IF (x .GE. 0) THEN
    xmin = (0.5d0+EPSILON(1d0)) * 10d0**(2-w)
    xmax = 10d0**w - 1d0
    IF (x .LT. xmin) THEN
      y = 0d0
    ELSE IF (x .GT. xmax) THEN
      y = xmax
    END IF
  ELSE
    xmin = (0.5d0+EPSILON(1d0)) * 10d0**(3-w)
    xmax = 10d0**(w-1) - 1d0
    IF (ABS(x) .LT. xmin) THEN
      y = 0d0
    ELSE IF (ABS(x) .GT. xmax) THEN
      y = -xmax
    END IF
  END IF
  RETURN
END FUNCTION


! ----------------------------------------------------------------------
! Returns the given number coerced into a range with an exponent of at
! most N digits.  Large magnitude numbers are reduced and very small
! magnitude numbers are set to zero.  The number of significant digits
! Ns can optionally be given (controls how many 9's appear in upper
! cutoff).
! 
ELEMENTAL FUNCTION CoerceExponent(x,N,Ns) RESULT(y)
  IMPLICIT NONE
  REAL*8 :: y
  REAL*8, INTENT(IN) :: x
  INTEGER, INTENT(IN) :: N
  INTEGER, INTENT(IN), OPTIONAL :: Ns
  REAL*8 :: xmin,xmax
  xmin = 10d0**(-10**N+1)
  IF (PRESENT(Ns)) THEN
    xmax = (1-0.1d0**Ns)*10d0**(10**N)
  ELSE
    xmax = 0.9*10d0**(10**N)
  END IF
  IF (ABS(x) .LT. xmin) THEN
    y = 0d0
  ELSE IF (ABS(x) .GE. xmax) THEN
    y = SIGN(1d0,x)*xmax
  ELSE
    y = x
  END IF
  RETURN
END FUNCTION


!=======================================================================
! ARRAY SORT/SEARCH
!=======================================================================

!-----------------------------------------------------------------------
! Sorts the given array using the heap sort method.
! Faster than quicksort when arrays are close to sorted already.
!   N           length of array
!   x           array of data to be sorted
! 
PURE SUBROUTINE HSort(N,x)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: N
  REAL*8, INTENT(INOUT) :: x(N)
  INTEGER :: I,Itmp
  REAL*8 :: xtmp
  DO I = N/2, 1, -1
    CALL HSort_siftdown(N,x,I,N)
  END DO
  DO I = N, 2, -1
    xtmp = x(1)
    x(1) = x(I)
    x(I) = xtmp
    CALL HSort_siftdown(N,x,1,I-1)
  END DO
  RETURN
    
  CONTAINS
  ! ------------------------------------------
  ! Utility function for heap sort
  PURE SUBROUTINE HSort_siftdown(N,x,Jlow,Jhigh)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: N,Jlow,Jhigh
    REAL*8, INTENT(INOUT) :: x(N)
    INTEGER :: J,Jold
    REAL*8 :: x0
    x0 = x(Jlow)
    Jold = Jlow
    J = 2*Jlow
    DO WHILE (J .LE. Jhigh)
      IF ((J .LT. Jhigh) .AND. (x(J)) .LT. x(J+1)) THEN
        J = J+1
      END IF
      IF (x0 .GE. x(J)) EXIT
      x(Jold) = x(J)
      Jold = J
      J = 2*J
    END DO
    x(Jold) = x0
    RETURN
  END SUBROUTINE HSort_siftdown
END SUBROUTINE HSort


! ----------------------------------------------------------------------
! Searches the given array of size [1:N], assumed to be sorted in
! increasing order, for the given value using a binary search algorithm.
! Returns the index I such that x(I) <= x0 < x(I+1), 0 if x0 < x(1),
! or N if x(N) <= x0.
! 
! Input arguments:
!   N               Size of x array
!   x               Array of data to be searched (must be sorted)
!   x0              Value to search for
! Optional input arguments:
!   Istart          Index to start searching from
! 
PURE FUNCTION BSearch(N,x,x0,Istart) RESULT(index)
  IMPLICIT NONE
  INTEGER :: index
  INTEGER, INTENT(IN) :: N
  REAL*8, INTENT(IN) :: x(N),x0
  INTEGER, INTENT(IN), OPTIONAL :: Istart
  INTEGER :: Ilow,Ihigh,Imid,step,dir
  
  ! Check if in bounds
  IF (x0 .LT. x(1)) THEN
    index = 0
    RETURN
  ELSE IF (x0 .GE. x(N)) THEN
    index = N
    RETURN
  END IF
  
  ! If starting index given, bracket desired index
  IF (PRESENT(Istart)) THEN
    ! Find one bound
    IF (Istart .LE. 1) THEN
      Ilow = 1
      dir = +1
    ELSE IF (Istart .LT. N) THEN
      IF (x0 .GE. x(Istart)) THEN
        Ilow = Istart
        dir = +1
      ELSE
        Ihigh = Istart
        dir = -1
      END IF
    ELSE
      Ihigh = N
      dir = -1
    END IF
    ! Search up or down in increasing step sizes to find other bound
    step = 1
    IF (dir .GT. 0) THEN
      Ihigh = MIN(Ilow + step,N)
      DO WHILE (x0 .GT. x(Ihigh))
        Ilow  = Ihigh
        Ihigh = MIN(Ilow + step,N)
        step  = 2*step
      END DO
    ELSE IF (dir .LT. 0) THEN
      Ilow = MAX(Ihigh - step,1)
      DO WHILE (x0 .LT. x(Ilow))
        Ihigh = Ilow
        Ilow  = MAX(Ilow - step,1)
        step  = 2*step
      END DO
    END IF
  ! No starting index given, search entire array
  ELSE
    Ilow  = 1
    Ihigh = N
  END IF
  
  ! Binary search
  DO WHILE (Ihigh-Ilow .GT. 1)
    Imid = (Ilow+Ihigh)/2
    IF (x0 .GE. x(Imid)) THEN
      Ilow = Imid
    ELSE
      Ihigh = Imid
    END IF
  END DO
  index = Ilow
  
END FUNCTION BSearch



!=======================================================================
! STRING UTILITY FUNCTIONS
!=======================================================================

! ----------------------------------------------------------------------
! Determines number of fields in the given string.
! 
PURE FUNCTION NumberOfFields(in) RESULT(N)
  IMPLICIT NONE
  INTEGER :: N
  CHARACTER*(*), INTENT(IN) :: in
  INTEGER :: I,len,pos
  LOGICAL :: is_blank
  
  len = LEN_TRIM(in)
  
  pos      = 0
  N        = 0
  is_blank = .TRUE.
  
  DO I = 1, len
    IF (is_blank .AND. (in(I:I) .NE. ' ')) THEN
      N = N + 1
      is_blank = .FALSE.
    ELSE
      is_blank = (in(I:I) .EQ. ' ')
    END IF
  END DO
  
END FUNCTION


! ----------------------------------------------------------------------
! Changes string to uppercase.
! 
PURE SUBROUTINE ToUppercase(string)
  IMPLICIT NONE
  CHARACTER*(*), INTENT(INOUT) :: string
  INTEGER, PARAMETER :: UA = ICHAR('A')
  INTEGER, PARAMETER :: la = ICHAR('a')
  INTEGER, PARAMETER :: UZ = ICHAR('Z')
  INTEGER, PARAMETER :: lz = ICHAR('z')
  INTEGER :: I,ival
  
  DO I = 1, LEN_TRIM(string)
    ival = ICHAR(string(I:I))
    IF ((ival .GE. la) .AND. (ival .LE. lz)) THEN
      string(I:I) = CHAR(ival+UA-la)
    END IF
  END DO
  
END SUBROUTINE


! ----------------------------------------------------------------------
! Changes string to lowercase.
! 
PURE SUBROUTINE ToLowercase(string)
  IMPLICIT NONE
  CHARACTER*(*), INTENT(INOUT) :: string
  INTEGER, PARAMETER :: UA = ICHAR('A')
  INTEGER, PARAMETER :: la = ICHAR('a')
  INTEGER, PARAMETER :: UZ = ICHAR('Z')
  INTEGER, PARAMETER :: lz = ICHAR('z')
  INTEGER :: I,ival

  DO I = 1, LEN_TRIM(string)
    ival = ICHAR(string(I:I))
    IF ((ival .GE. UA) .AND. (ival .LE. UZ)) THEN
      string(I:I) = CHAR(ival-UA+la)
    END IF
  END DO

END SUBROUTINE


! ----------------------------------------------------------------------
! Obtain the file suffix or an empty string if there is no suffix.
! 
PURE SUBROUTINE FileSuffix(file,suffix)
  IMPLICIT NONE
  CHARACTER*(*), INTENT(IN) :: file
  CHARACTER*(*), INTENT(OUT) :: suffix
  INTEGER :: len,pos
  
  len = LEN_TRIM(file)
  pos = INDEX(file,'.',back=.TRUE.)
  
  ! No suffix
  IF (pos .EQ. 0) THEN
    suffix = ''
    RETURN
  END IF
  
  ! Get suffix
  suffix = file(pos+1:)
  
END SUBROUTINE


! ----------------------------------------------------------------------
! Obtain the file suffix or an empty string if there is no suffix.
! The file suffix is converted to lower case for easier parsing.
! 
PURE SUBROUTINE FileSuffixLC(file,suffix)
  IMPLICIT NONE
  CHARACTER*(*), INTENT(IN) :: file
  CHARACTER*(*), INTENT(OUT) :: suffix
  INTEGER :: len,pos
  
  len = LEN_TRIM(file)
  pos = INDEX(file,'.',back=.TRUE.)
  
  ! No suffix
  IF (pos .EQ. 0) THEN
    suffix = ''
    RETURN
  END IF
  
  ! Get suffix, convert to lower case for easier parsing
  suffix = file(pos+1:)
  CALL ToLowercase(suffix)
  
END SUBROUTINE



!=======================================================================
! TIME UTILITY FUNCTIONS
!=======================================================================

!-----------------------------------------------------------------------
! System time in seconds.
! 
FUNCTION GetTime() RESULT(t)
  IMPLICIT NONE
  REAL*8 :: t
  ! Resolution may vary by kind; ensure highest resolution
  INTEGER(KIND=8) :: count,count_rate
  CALL SYSTEM_CLOCK(count,count_rate)
  t = count*1d0 / count_rate
END FUNCTION


!-----------------------------------------------------------------------
! Time elapsed in seconds since TimeElapsedReset() was called.
! 
FUNCTION TimeElapsed() RESULT(t)
  IMPLICIT NONE
  REAL*8 :: t
  CALL TimeElapsedRoutine(t,.FALSE.)
END FUNCTION


!-----------------------------------------------------------------------
! Reset TimeElapsed timer.  Subroutine version of above.
! 
SUBROUTINE TimeElapsedReset()
  IMPLICIT NONE
  REAL*8 :: t
  CALL TimeElapsedRoutine(t,.TRUE.)
END SUBROUTINE


!-----------------------------------------------------------------------
! Time elapsed in seconds since the parameter 'reset' was set to true.
! Used by TimeElapsed() and TimeElapsedReset() routines, which should
! be used instead of direct calls to this routine.  Note that this is
! wall time, not CPU time.
! 
!   t           Time elpased
!   reset       Logical specifying if the counter should be
!               reset to the current cpu time.
! 
SUBROUTINE TimeElapsedRoutine(t,reset)
  IMPLICIT NONE
  REAL*8, INTENT(OUT) :: t
  LOGICAL,INTENT(IN) :: reset
  ! Resolution may vary by kind; ensure highest resolution
  INTEGER(KIND=8), SAVE :: icount = 0
  INTEGER(KIND=8) :: fcount,count_rate
  IF (reset) THEN
    CALL SYSTEM_CLOCK(icount)
    t = 0d0
    RETURN
  END IF
  CALL SYSTEM_CLOCK(fcount,count_rate)
  t = (fcount-icount)*1d0 / count_rate
END SUBROUTINE


!-----------------------------------------------------------------------
! CPU time elapsed in seconds since CPUTimeElapsedReset() was called.
! For multithreaded cases, this is the sum of time spent by all threads.
! 
FUNCTION CPUTimeElapsed() RESULT(t)
  IMPLICIT NONE
  REAL*8 :: t
  CALL CPUTimeElapsedRoutine(t,.FALSE.)
END FUNCTION


!-----------------------------------------------------------------------
! Reset CPUTimeElapsed timer.
! 
SUBROUTINE CPUTimeElapsedReset()
  IMPLICIT NONE
  REAL*8 :: t
  CALL CPUTimeElapsedRoutine(t,.TRUE.)
END SUBROUTINE


!-----------------------------------------------------------------------
! CPU time elapsed in seconds since the parameter 'reset' was set to
! true.  Used by CPUTimeElapsed() and CPUTimeElapsedReset() routines,
! which should be used instead of direct calls to this routine.
! 
!   t           Time elpased
!   reset       Logical specifying if the counter should be
!               reset to the current cpu time.
! 
SUBROUTINE CPUTimeElapsedRoutine(t,reset)
  IMPLICIT NONE
  REAL*8, INTENT(OUT) :: t
  LOGICAL,INTENT(IN) :: reset
  ! Resolution may vary by kind; ensure highest resolution
  REAL(KIND=8), SAVE :: t1
  REAL(KIND=8) :: t2
  IF (reset) THEN
    CALL CPU_TIME(t1)
    t = 0d0
  END IF
  CALL CPU_TIME(t2)
  t = t2 - t1
END SUBROUTINE


END MODULE
