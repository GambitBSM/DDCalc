MODULE DDCommandLine

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! DDCommandLine
!    Utility routines for processing command line arguments, mainly
!    of the form '-x', '--flag', or '--key=value'.
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

IMPLICIT NONE

CONTAINS


!=======================================================================
! There are several types of routine in this module:
!   *) GetXXXArg(index,status)
!      Parses the command line argument at the given index and returns
!      as type XXX (logical,integer,real,complex).  A negative index
!      will count from the end of the argument list.  The optional
!      'status' parameter indicates if the argument exists and was
!      properly parsed.
!   *) GetShortArg(akey)
!      Indicates if the given akey appears in any argument of the form
!      -<akey> or -<akey1><akey2><akey3>....  In this case, akey
!      is meant to be a single character.
!   *) GetLongArg(akey)
!      Indicates if the given akey appears in any argument of the form
!      --<akey> or -<akey>=<aval>.
!   *) GetLongArgXXX(akey,aval)
!      Indicates if the given akey appears in any argument of the form
!      --<akey> or -<akey>=<aval>.  If present, <aval> will be parsed
!      as type XXX (string,logical,integer,real,complex) and placed in
!      aval.  The first instance of a particular key will always be
!      parsed; if additional arguments use the same key, they will be
!      ignored.  This routine returns .FALSE. if there is a failure to
!      parse <aval>.
!=======================================================================

! ----------------------------------------------------------------------
! Obtains the command line argument at the given index as a LOGICAL.
! Negative index counts from the end of the argument list.
! The optional status argument indicates if the given argument exists
! and was properly parsed.
! 
FUNCTION GetLogicalArg(index,status) RESULT(aval)
  IMPLICIT NONE
  LOGICAL :: aval
  INTEGER, INTENT(IN) :: index
  LOGICAL, INTENT(OUT), OPTIONAL :: status
  CHARACTER*256 :: arg
  INTEGER :: index0,Narg,ios
  
  aval = .FALSE.
  IF (PRESENT(status)) status = .FALSE.
  
  ! Check index validity, count from end if necessary
  Narg   = IARGC()
  IF (index .LT. 0) THEN
    index0 = Narg + 1 + index
  ELSE
    index0 = index
  END IF
  IF ((index0 .LT. 1) .OR. (index0 .GT. Narg)) RETURN
  
  ! Parse argument
  CALL GETARG(index0,arg)
  READ(UNIT=arg,FMT=*,IOSTAT=ios) aval
  IF (ios .EQ. 0) THEN
    IF (PRESENT(status)) status = .TRUE.
  END IF
  
END FUNCTION GetLogicalArg


! ----------------------------------------------------------------------
! Obtains the command line argument at the given index as an INTEGER.
! Negative index counts from the end of the argument list.
! The optional status argument indicates if the given argument exists
! and was properly parsed.
! Accepts real valued strings, which will be rounded to the nearest
! integer.
! 
FUNCTION GetIntegerArg(index,status) RESULT(aval)
  IMPLICIT NONE
  INTEGER :: aval
  INTEGER, INTENT(IN) :: index
  LOGICAL, INTENT(OUT), OPTIONAL :: status
  CHARACTER*256 :: arg
  INTEGER :: index0,Narg,ios
  REAL*8 :: x
  
  aval = -HUGE(aval)
  IF (PRESENT(status)) status = .FALSE.
  
  ! Check index validity, count from end if necessary
  Narg   = IARGC()
  IF (index .LT. 0) THEN
    index0 = Narg + 1 + index
  ELSE
    index0 = index
  END IF
  IF ((index0 .LT. 1) .OR. (index0 .GT. Narg)) RETURN
  
  ! Parse argument
  ! Try reading as integer first, then as real and convert
  ! to integer
  CALL GETARG(index0,arg)
  READ(UNIT=arg,FMT=*,IOSTAT=ios) aval
  IF (ios .EQ. 0) THEN
    IF (PRESENT(status)) status = .TRUE.
  ELSE
    READ(UNIT=arg,FMT=*,IOSTAT=ios) x
    IF (ios .EQ. 0) THEN
      IF (PRESENT(status)) status = .TRUE.
      aval = NINT(x)
    END IF
  END IF
  
END FUNCTION GetIntegerArg


! ----------------------------------------------------------------------
! Obtains the command line argument at the given index as a REAL*8.
! Negative index counts from the end of the argument list.
! The optional status argument indicates if the given argument exists
! and was properly parsed.
! 
FUNCTION GetRealArg(index,status) RESULT(aval)
  IMPLICIT NONE
  REAL*8 :: aval
  INTEGER, INTENT(IN) :: index
  LOGICAL, INTENT(OUT), OPTIONAL :: status
  CHARACTER*256 :: arg
  INTEGER :: index0,Narg,ios
  
  aval = -HUGE(aval)
  IF (PRESENT(status)) status = .FALSE.
  
  ! Check index validity, count from end if necessary
  Narg   = IARGC()
  IF (index .LT. 0) THEN
    index0 = Narg + 1 + index
  ELSE
    index0 = index
  END IF
  IF ((index0 .LT. 1) .OR. (index0 .GT. Narg)) RETURN
  
  ! Parse argument
  CALL GETARG(index0,arg)
  READ(UNIT=arg,FMT=*,IOSTAT=ios) aval
  IF (ios .EQ. 0) THEN
    IF (PRESENT(status)) status = .TRUE.
  END IF
  
END FUNCTION GetRealArg


! ----------------------------------------------------------------------
! Obtains the command line argument at the given index as a COMPLEX*16.
! Negative index counts from the end of the argument list.
! The optional status argument indicates if the given argument exists
! and was properly parsed.
! 
FUNCTION GetComplexArg(index,status) RESULT(aval)
  IMPLICIT NONE
  COMPLEX*16 :: aval
  INTEGER, INTENT(IN) :: index
  LOGICAL, INTENT(OUT), OPTIONAL :: status
  CHARACTER*256 :: arg
  INTEGER :: index0,Narg,ios
  
  aval = -HUGE(1d0)*(1,1)
  IF (PRESENT(status)) status = .FALSE.
  
  ! Check index validity, count from end if necessary
  Narg   = IARGC()
  IF (index .LT. 0) THEN
    index0 = Narg + 1 + index
  ELSE
    index0 = index
  END IF
  IF ((index0 .LT. 1) .OR. (index0 .GT. Narg)) RETURN
  
  ! Parse argument
  CALL GETARG(index0,arg)
  READ(UNIT=arg,FMT=*,IOSTAT=ios) aval
  IF (ios .EQ. 0) THEN
    IF (PRESENT(status)) status = .TRUE.
  END IF
  
END FUNCTION GetComplexArg


! ----------------------------------------------------------------------
! Checks if the given argument key is given as a command line argument
! of the form -akey or -<akey1><akey2><akey3>...
! Returns .TRUE. if valid key was found, .FALSE. otherwise.
! 
FUNCTION GetShortArg(akey) RESULT(status)
  IMPLICIT NONE
  LOGICAL :: status
  CHARACTER*(*), INTENT(IN) :: akey
  CHARACTER*256 :: arg
  INTEGER :: I,Narg,pos
  
  Narg   = IARGC()
  status = .FALSE.
  
  DO I=1,Narg
    CALL GETARG(I,arg)
    IF (LEN_TRIM(arg) .LT. 2)  CYCLE
    IF (arg(1:1) .NE. '-')  CYCLE
    IF (arg(2:2) .EQ. '-')  CYCLE
    pos = INDEX(arg,akey)
    IF (pos .GE. 2) THEN
      status = .TRUE.
      RETURN
    END IF
  END DO
  
END FUNCTION GetShortArg


! ----------------------------------------------------------------------
! Checks if the given argument key is given as a command line argument
! of the form --akey or --akey=<val>.
! Returns .TRUE. if valid key was found, .FALSE. otherwise.
! 
FUNCTION GetLongArg(akey) RESULT(status)
  IMPLICIT NONE
  LOGICAL :: status
  CHARACTER*(*), INTENT(IN) :: akey
  CHARACTER*256 :: arg
  INTEGER :: I,Narg,len
  
  Narg   = IARGC()
  len    = LEN_TRIM(akey)
  status = .FALSE.
  
  DO I=1,Narg
    CALL GETARG(I,arg)
    IF (arg .EQ. '--' // akey) THEN
      status = .TRUE.
      RETURN
    ELSE IF (arg(:3+len) .EQ. '--' // akey // '=') THEN
      status = .TRUE.
      RETURN
    END IF
  END DO
  
END FUNCTION GetLongArg


! ----------------------------------------------------------------------
! Checks if the given argument key is given as a command line argument
! of the form --akey[<suffix>] or --akey[<suffix>]=<val>.
! This differs from above in that it indicates if there are any command
! line arguments where the key name begins with the given string.
! Returns .TRUE. if valid key was found, .FALSE. otherwise.
! 
FUNCTION GetLongArgPrefix(prefix) RESULT(status)
  IMPLICIT NONE
  LOGICAL :: status
  CHARACTER*(*), INTENT(IN) :: prefix
  CHARACTER*256 :: arg
  INTEGER :: I,Narg,len
  
  Narg   = IARGC()
  len    = LEN_TRIM(prefix)
  status = .FALSE.
  
  DO I=1,Narg
    CALL GETARG(I,arg)
    IF (arg(:2+len) .EQ. '--' // prefix) THEN
      status = .TRUE.
      RETURN
    END IF
  END DO
  
END FUNCTION GetLongArgPrefix


! ----------------------------------------------------------------------
! Checks if the given argument key is given as a command line argument
! of the form --akey=<aval> and inserts the appropriate string into
! aval.
! Returns .TRUE. if valid key was found, .FALSE. otherwise.
! 
FUNCTION GetLongArgString(akey,aval) RESULT(status)
  IMPLICIT NONE
  LOGICAL :: status
  CHARACTER*(*), INTENT(IN) :: akey
  CHARACTER*(*), INTENT(OUT) :: aval
  CHARACTER*1024 :: arg
  INTEGER :: I,Narg,len,pos
  
  Narg   = IARGC()
  len    = LEN_TRIM(akey)
  status = .FALSE.
  
  DO I=1,Narg
    CALL GETARG(I,arg)
    IF (arg(:2+len) .EQ. '--' // akey) THEN
      pos = INDEX(arg,'=')
      IF (pos .EQ. 3+len) THEN
        status = .TRUE.
        aval = TRIM(arg(pos+1:))
        RETURN
      END IF
    END IF
  END DO
  
END FUNCTION GetLongArgString


! ----------------------------------------------------------------------
! Checks if the given argument key is given as a command line argument
! of the form --akey=<aval1>,<aval2>,... and inserts the appropriate
! CHARACTER*N values into aval (ALLOCATABLE array of length Nval).
! 
! Currently requires allocatable array of fixed length strings,
! declared as:
!   CHARACTER(LEN=N), DIMENSION(:), ALLOCATABLE :: aval
! 
! A better implementation (but requires gfortran 4.6+) uses deferred-
! length strings, where aval is declared as:
!   CHARACTER(LEN=:), DIMENSION(:), ALLOCATABLE :: aval
! 
! Returns .TRUE. if valid key was found, .FALSE. otherwise.
! 
! NOTE: This function must be defined within a module or the use of
!       interfaces will be required.
! 
FUNCTION GetLongArgStrings(akey,N,aval,Nval) RESULT(status)
  IMPLICIT NONE
  LOGICAL :: status
  CHARACTER*(*), INTENT(IN) :: akey
  INTEGER, INTENT(IN) :: N
  ! Fixed length strings
  CHARACTER(LEN=N), DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: aval
  ! Better: use deferred length strings (but needs gfortran 4.6+)
  !CHARACTER(LEN=:), DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: aval
  INTEGER, INTENT(OUT) :: Nval
  CHARACTER*1024 :: sval
  INTEGER :: I1,I2,Ia,ios
  
  IF (.NOT. GetLongArgString(akey,sval)) THEN
    status = .FALSE.
    RETURN
  END IF
  
  status = .TRUE.
  
  ! Find number of values
  I1   = 1
  Nval = 1
  DO WHILE (INDEX(sval(I1:),',') .GT. 0)
    I1 = I1 + INDEX(sval(I1:),',')
    Nval = Nval + 1
  END DO
  IF (ALLOCATED(aval)) DEALLOCATE(aval)
  ALLOCATE(CHARACTER(LEN=N) :: aval(1:Nval))
  aval = ''
  
  ! Parse values
  READ(UNIT=sval(I1:),FMT='(A)',IOSTAT=ios) aval(Nval)
  IF (ios .NE. 0) THEN
    status = .FALSE.
    Nval   = 0
    IF (ALLOCATED(aval)) DEALLOCATE(aval)
    ALLOCATE(CHARACTER(LEN=N) :: aval(0))
    RETURN
  END IF
  I2 = 1
  DO Ia = 1, Nval-1
    I1 = I2
    I2 = I1 + INDEX(sval(I1:),',')
    READ(UNIT=sval(I1:I2-2),FMT='(A)',IOSTAT=ios) aval(Ia)
    IF (ios .NE. 0) THEN
      status = .FALSE.
      Nval   = 0
      IF (ALLOCATED(aval)) DEALLOCATE(aval)
      ALLOCATE(CHARACTER(LEN=N) :: aval(0))
      RETURN
    END IF
  END DO
  
END FUNCTION GetLongArgStrings


! ----------------------------------------------------------------------
! Checks if the given argument key is given as a command line argument
! of the form --akey=<aval> and inserts the appropriate LOGICAL value
! into aval.
! Returns .TRUE. if valid key was found, .FALSE. otherwise.
! 
FUNCTION GetLongArgLogical(akey,aval) RESULT(status)
  IMPLICIT NONE
  LOGICAL :: status
  CHARACTER*(*), INTENT(IN) :: akey
  LOGICAL, INTENT(OUT) :: aval
  CHARACTER*256 :: sval
  INTEGER :: ios
  IF (.NOT. GetLongArgString(akey,sval)) THEN
    status = .FALSE.
    RETURN
  END IF
  READ(UNIT=sval,FMT=*,IOSTAT=ios) aval
  status = (ios .EQ. 0)
END FUNCTION GetLongArgLogical


! ----------------------------------------------------------------------
! Checks if the given argument key is given as a command line argument
! of the form --akey=<aval1>,<aval2>,... and inserts the appropriate
! LOGICAL values into aval (ALLOCATABLE array of length Nval).
! Returns .TRUE. if valid key was found, .FALSE. otherwise.
! 
! NOTE: This function must be defined within a module or the use of
!       interfaces will be required.
! 
FUNCTION GetLongArgLogicals(akey,aval,Nval) RESULT(status)
  IMPLICIT NONE
  LOGICAL :: status
  CHARACTER*(*), INTENT(IN) :: akey
  LOGICAL, ALLOCATABLE, INTENT(OUT) :: aval(:)
  INTEGER, INTENT(OUT) :: Nval
  ! Fixed length strings
  INTEGER, PARAMETER :: N = 64
  CHARACTER(LEN=N), ALLOCATABLE :: sval(:)
  ! Better: use deferred length strings (but needs gfortran 4.6+)
  !CHARACTER(LEN=:), ALLOCATABLE :: sval(:)
  INTEGER :: I,ios
  
  IF (.NOT. GetLongArgStrings(akey,N,sval,Nval)) THEN
    status = .FALSE.
    Nval   = 0
    IF (ALLOCATED(aval)) DEALLOCATE(aval)
    ALLOCATE(aval(0))
    RETURN
  END IF
  
  status = .TRUE.
  IF (ALLOCATED(aval)) DEALLOCATE(aval)
  ALLOCATE(aval(1:Nval))
  
  DO I=1,Nval
    READ(UNIT=sval(I),FMT=*,IOSTAT=ios) aval(I)
    IF (ios .NE. 0) THEN
      status = .FALSE.
      Nval   = 0
      IF (ALLOCATED(aval)) DEALLOCATE(aval)
      ALLOCATE(aval(0))
      RETURN
    END IF
  END DO
  
END FUNCTION GetLongArgLogicals


! ----------------------------------------------------------------------
! Checks if the given argument key is given as a command line argument
! of the form --akey=<aval> and inserts the appropriate INTEGER value
! into aval.
! Accepts real valued strings, which will be rounded to the nearest
! integer.
! Returns .TRUE. if valid key was found, .FALSE. otherwise.
! 
FUNCTION GetLongArgInteger(akey,aval) RESULT(status)
  IMPLICIT NONE
  LOGICAL :: status
  CHARACTER*(*), INTENT(IN) :: akey
  INTEGER, INTENT(OUT) :: aval
  CHARACTER*256 :: sval
  INTEGER :: ios
  REAL*8 :: x
  IF (.NOT. GetLongArgString(akey,sval)) THEN
    status = .FALSE.
    RETURN
  END IF
  ! Try reading as integer first, then as real and convert to integer
  READ(UNIT=sval,FMT=*,IOSTAT=ios) aval
  IF (ios .NE. 0) THEN
    READ(UNIT=sval,FMT=*,IOSTAT=ios) x
    IF (ios .EQ. 0) aval = NINT(x)
  END IF
  status = (ios .EQ. 0)
END FUNCTION GetLongArgInteger


! ----------------------------------------------------------------------
! Checks if the given argument key is given as a command line argument
! of the form --akey=<aval1>,<aval2>,... and inserts the appropriate
! INTEGER values into aval (ALLOCATABLE array of length Nval).
! Returns .TRUE. if valid key was found, .FALSE. otherwise.
! 
! NOTE: This function must be defined within a module or the use of
!       interfaces will be required.
! 
FUNCTION GetLongArgIntegers(akey,aval,Nval) RESULT(status)
  IMPLICIT NONE
  LOGICAL :: status
  CHARACTER*(*), INTENT(IN) :: akey
  INTEGER, ALLOCATABLE, INTENT(OUT) :: aval(:)
  INTEGER, INTENT(OUT) :: Nval
  ! Fixed length strings
  INTEGER, PARAMETER :: N = 64
  CHARACTER(LEN=N), ALLOCATABLE :: sval(:)
  ! Better: use deferred length strings (but needs gfortran 4.6+)
  !CHARACTER(LEN=:), ALLOCATABLE :: sval(:)
  INTEGER :: I,ios
  
  IF (.NOT. GetLongArgStrings(akey,N,sval,Nval)) THEN
    status = .FALSE.
    Nval   = 0
    IF (ALLOCATED(aval)) DEALLOCATE(aval)
    ALLOCATE(aval(0))
    RETURN
  END IF
  
  status = .TRUE.
  IF (ALLOCATED(aval)) DEALLOCATE(aval)
  ALLOCATE(aval(1:Nval))
  
  DO I=1,Nval
    READ(UNIT=sval(I),FMT=*,IOSTAT=ios) aval(I)
    IF (ios .NE. 0) THEN
      status = .FALSE.
      Nval   = 0
      IF (ALLOCATED(aval)) DEALLOCATE(aval)
      ALLOCATE(aval(0))
      RETURN
    END IF
  END DO
  
END FUNCTION GetLongArgIntegers


! ----------------------------------------------------------------------
! Checks if the given argument key is given as a command line argument
! of the form --akey=<aval> and inserts the appropriate REAL*8 value
! into aval.
! Returns .TRUE. if valid key was found, .FALSE. otherwise.
! 
FUNCTION GetLongArgReal(akey,aval) RESULT(status)
  IMPLICIT NONE
  LOGICAL :: status
  CHARACTER*(*), INTENT(IN) :: akey
  REAL*8, INTENT(OUT) :: aval
  CHARACTER*256 :: sval
  INTEGER :: ios
  IF (.NOT. GetLongArgString(akey,sval)) THEN
    status = .FALSE.
    RETURN
  END IF
  READ(UNIT=sval,FMT=*,IOSTAT=ios) aval
  status = (ios .EQ. 0)
END FUNCTION GetLongArgReal


! ----------------------------------------------------------------------
! Checks if the given argument key is given as a command line argument
! of the form --akey=<aval1>,<aval2>,... and inserts the appropriate
! REAL*8 values into aval (ALLOCATABLE array of length Nval).
! Returns .TRUE. if valid key was found, .FALSE. otherwise.
! 
! NOTE: This function must be defined within a module or the use of
!       interfaces will be required.
! 
FUNCTION GetLongArgReals(akey,aval,Nval) RESULT(status)
  IMPLICIT NONE
  LOGICAL :: status
  CHARACTER*(*), INTENT(IN) :: akey
  REAL*8, ALLOCATABLE, INTENT(OUT) :: aval(:)
  INTEGER, INTENT(OUT) :: Nval
  ! Fixed length strings
  INTEGER, PARAMETER :: N = 64
  CHARACTER(LEN=N), ALLOCATABLE :: sval(:)
  ! Better: use deferred length strings (but needs gfortran 4.6+)
  !CHARACTER(LEN=:), ALLOCATABLE :: sval(:)
  INTEGER :: I,ios
  
  IF (.NOT. GetLongArgStrings(akey,N,sval,Nval)) THEN
    status = .FALSE.
    Nval   = 0
    IF (ALLOCATED(aval)) DEALLOCATE(aval)
    ALLOCATE(aval(0))
    RETURN
  END IF
  
  status = .TRUE.
  IF (ALLOCATED(aval)) DEALLOCATE(aval)
  ALLOCATE(aval(1:Nval))
  
  DO I=1,Nval
    READ(UNIT=sval(I),FMT=*,IOSTAT=ios) aval(I)
    IF (ios .NE. 0) THEN
      status = .FALSE.
      Nval   = 0
      IF (ALLOCATED(aval)) DEALLOCATE(aval)
      ALLOCATE(aval(0))
      RETURN
    END IF
  END DO
  
END FUNCTION GetLongArgReals


! ----------------------------------------------------------------------
! Checks if the given argument key is given as a command line argument
! of the form --akey=<aval> and inserts the appropriate COMPLEX*16 value
! into aval.
! Returns .TRUE. if valid key was found, .FALSE. otherwise.
! 
FUNCTION GetLongArgComplex(akey,aval) RESULT(status)
  IMPLICIT NONE
  LOGICAL :: status
  CHARACTER*(*), INTENT(IN) :: akey
  COMPLEX*16, INTENT(OUT) :: aval
  CHARACTER*256 :: sval
  INTEGER :: ios
  IF (.NOT. GetLongArgString(akey,sval)) THEN
    status = .FALSE.
    RETURN
  END IF
  READ(UNIT=sval,FMT=*,IOSTAT=ios) aval
  status = (ios .EQ. 0)
END FUNCTION GetLongArgComplex


! ----------------------------------------------------------------------
! Returns the full command line used to run the program.
! This routine is nearly identical to the Fortran 2003 routine
! GET_COMMAND, but is provided here because some older compilers
! do not provide that routine.  White space between arguments will
! not be conserved: arguments will be separated by a single space,
! regardless of actual white space on the command line.
! 
SUBROUTINE GetFullCommand(cmd)
  IMPLICIT NONE
  CHARACTER*(*), INTENT(OUT) :: cmd
  INTEGER :: pos,I,Narg,alen,clen
  CHARACTER*256 :: arg
  
  ! Fortran 2003
  !CALL GET_COMMAND(cmd)
  
  cmd = ''
  pos = 1
  clen = LEN(cmd)
  
  ! Command
  !CALL GETARG(0,arg,alen)
  CALL GETARG(0,arg)
  alen = LEN_TRIM(arg)
  cmd = arg
  pos = pos + alen
  
  ! Arguments
  Narg = IARGC()
  DO I=1, Narg
    !CALL GETARG(I,arg,alen)
    CALL GETARG(I,arg)
    alen = LEN_TRIM(arg)
    cmd(pos:pos) = ' '
    pos = pos + 1
    IF (pos+alen-1 .GT. clen) THEN
      cmd(pos:) = arg
    ELSE
      cmd(pos:pos+alen-1) = arg(1:alen)
    END IF
    pos = pos + alen
    IF (pos .GE. clen) EXIT
  END DO
  
  ! Following necessary to terminate/clear string
  cmd(pos:) = ''
  
END SUBROUTINE GetFullCommand


! ----------------------------------------------------------------------
! Checks if the command line contains any requests to show the program
! usage.
! 
FUNCTION ShowUsageRequested() RESULT(flag)
  IMPLICIT NONE
  LOGICAL :: flag
  INTEGER :: I,Narg
  CHARACTER*32 :: arg
  
  flag = .FALSE.
  Narg = IARGC()
  
  DO I=1,Narg
    CALL GETARG(I,arg)
    IF (arg(1:1) .EQ. '?') THEN
      flag = .TRUE.
      RETURN
    ELSE IF ( (arg(1:2) .EQ. '-?') .OR. (arg(1:2) .EQ. '-h') ) THEN
      flag = .TRUE.
      RETURN
    ELSE IF (arg(1:6) .EQ. '--help') THEN
      flag = .TRUE.
      RETURN
    END IF
  END DO
  
END FUNCTION ShowUsageRequested


! ----------------------------------------------------------------------
! Shows usage.
! Prints out contents of given file.
! 
! Optional input argument:
!   file        File containing usage contents.  Default is the
!               calling program's name with '.help' appended, in the dir 'doc'.
!   
SUBROUTINE ShowUsage(file)
  IMPLICIT NONE
  CHARACTER*(*), INTENT(IN), OPTIONAL :: file
  CHARACTER*1024 :: file0
  INTEGER, PARAMETER :: USAGE_FID = 19
  CHARACTER*1024 :: line
  INTEGER :: ios,K
  
  IF (PRESENT(file)) THEN
    file0 = file
  ELSE
    CALL GETARG(0,file0)
    K = INDEX(file0,'/',.TRUE.)
    IF (K .NE. 0) file0 = file0(K+1:)
    file0 = TRIM(DDCALC_DIR) // '/doc/' // TRIM(file0) // '.help'
  END IF
  
  OPEN(FILE=file0,UNIT=USAGE_FID,STATUS='OLD',                          &
       FORM='FORMATTED',ACTION='READ',IOSTAT=ios)
  IF (ios .NE. 0) THEN
    WRITE(0,*) 'Unable to display usage:'
    WRITE(0,*) 'Could not open usage file ' // TRIM(file0)
    WRITE(0,*)
    RETURN
  END IF
  
  DO WHILE (.TRUE.)
    READ(UNIT=USAGE_FID,FMT='(A)',IOSTAT=ios) line
    IF (ios .LT. 0) EXIT
    WRITE(*,'(A)') TRIM(line)
  END DO
  
  CLOSE(UNIT=USAGE_FID,IOSTAT=ios)
  
END SUBROUTINE ShowUsage

END MODULE
