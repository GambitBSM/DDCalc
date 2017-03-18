MODULE DDInput

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! DDInput
!    DDCalc routines for import of data.
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

USE DDUtils

IMPLICIT NONE
PRIVATE

PUBLIC :: LoadArrays, LoadTable

! Prefixes to place at beginning of comment and data lines
CHARACTER*(*), PUBLIC, PARAMETER :: COMMENT_PREFIX = '# '
CHARACTER*(*), PUBLIC, PARAMETER :: DATA_PREFIX    = '  '
CHARACTER*(*), PUBLIC, PARAMETER :: COMMENT_LINE   = &
    '# ----------------------------------'   &
    // '------------------------------------'

CONTAINS


! ----------------------------------------------------------------------
! Returns an available (unused) file I/O unit number or -1 if no
! available unit number could be found.
! NOTE: A free unit number is a unit number not currently connected
!       at the time this function is called.  Numbers are not reserved,
!       so do not assume this number will remain free if it is not used
!       to open a file immediately.
! 
FUNCTION FreeIOUnit() RESULT(Nio)
  IMPLICIT NONE
  INTEGER :: Nio
  INTEGER, PARAMETER :: START_UNIT = 50
  INTEGER, PARAMETER :: MAX_UNIT   = 99
  LOGICAL opened
  
  ! Cycle through unit numbers until we find one that is not connected
  Nio = START_UNIT
  INQUIRE(UNIT=Nio,OPENED=opened)
  DO WHILE (opened .AND. (Nio .LT. MAX_UNIT))
    Nio = Nio + 1
    INQUIRE(UNIT=Nio,OPENED=opened)
  END DO
  
  ! Did not find unopen unit number in allowed range
  IF (opened) THEN
    WRITE(UNIT=5,FMT=*) 'ERROR: an error occurred while trying '        &
        // 'to find an available unit number.'
    Nio = -1
    RETURN
  END IF
  
END FUNCTION


! ----------------------------------------------------------------------
! Determines if the given line is a comment (non-data) line.
! A comment line is any blank line or line beginning with the comment
! character (default: '#').
!   line            the line to check
!   commentchar     (Optional) comment character.  Default is '#'.
! 
LOGICAL FUNCTION IsCommentLine(line,commentchar)
  IMPLICIT NONE
  CHARACTER*(*), INTENT(IN) :: line
  CHARACTER, INTENT(IN), OPTIONAL :: commentchar
  
  IsCommentLine = .FALSE.
  IF (LEN_TRIM(line) .EQ. 0) THEN
    IsCommentLine = .TRUE.
  ELSE IF (PRESENT(commentchar)) THEN
    IF (line(1:1) .EQ. commentchar) IsCommentLine = .TRUE.
  ELSE
    IF (line(1:1) .EQ. '#') IsCommentLine = .TRUE.
  END IF
  
END FUNCTION


! ----------------------------------------------------------------------
! Determines the number of lines in the given file.
!   fid             The identifier for the file to be read
!   status          (Optional) .TRUE. if file was successfully read
! 
FUNCTION FileLines(fid,status) RESULT(lines)
  IMPLICIT NONE
  INTEGER :: lines
  INTEGER, INTENT(IN) :: fid
  LOGICAL, INTENT(OUT), OPTIONAL :: status
  INTEGER :: ios
  
  lines = 0
  IF (PRESENT(status)) status = .FALSE.
  
  REWIND(UNIT=fid,IOSTAT=ios)
  IF (ios .LT. 0) RETURN
  
  DO WHILE (ios .GE. 0)
    lines = lines + 1
    READ(UNIT=fid,FMT=*,IOSTAT=ios)
  END DO
  
  ! Count included failed final read, so subtract 1
  lines = lines - 1
  
  REWIND(UNIT=fid,IOSTAT=ios)
  IF (PRESENT(status)) status = .TRUE.
  
END FUNCTION


! ----------------------------------------------------------------------
! Determines the number of data (non-comment) lines in the given file.
! Blank lines and lines beginning with the comment character (default:
! '#') are considered comment lines.
!   fid             The identifier for the file to be read
!   commentchar     (Optional) comment character.  Default is '#'.
!   status          (Optional) .TRUE. if file was successfully read
! 
FUNCTION FileDataLines(fid,commentchar,status) RESULT(lines)
  IMPLICIT NONE
  INTEGER :: lines
  INTEGER, INTENT(IN) :: fid
  CHARACTER, INTENT(IN), OPTIONAL :: commentchar
  LOGICAL, INTENT(OUT), OPTIONAL :: status
  INTEGER :: ios
  CHARACTER :: commentchar0
  CHARACTER*1024 :: line
  
  IF (PRESENT(commentchar)) THEN
    commentchar0 = commentchar
  ELSE
    commentchar0 = '#'
  END IF
  
  lines = 0
  IF (PRESENT(status)) status = .FALSE.
  
  REWIND(UNIT=fid,IOSTAT=ios)
  IF (ios .LT. 0) RETURN
  
  DO WHILE (ios .GE. 0)
    READ(UNIT=fid,FMT='(A)',IOSTAT=ios) line
    IF (ios .GE. 0) THEN
      IF (.NOT. IsCommentLine(line)) lines = lines + 1
    END IF
  END DO
  
  REWIND(UNIT=fid,IOSTAT=ios)
  IF (PRESENT(status)) status = .TRUE.
  
END FUNCTION


! ----------------------------------------------------------------------
! Determines the number of data columns in the given file.
! This is based upon the number of fields separated by spaces on a
! data line.
! NOTE: Currently limited to first 1024 characters in a line.
!       This can be increased below.
!   fid             The identifier for the file to be read
!   commentchar     (Optional) comment character.  Default is '#'.
!   status          (Optional) .TRUE. if file was successfully read
! 
FUNCTION FileDataColumns(fid,commentchar,status) RESULT(columns)
  IMPLICIT NONE
  INTEGER :: columns
  INTEGER, INTENT(IN) :: fid
  CHARACTER, INTENT(IN), OPTIONAL :: commentchar
  LOGICAL, INTENT(OUT), OPTIONAL :: status
  INTEGER :: ios
  CHARACTER :: commentchar0
  CHARACTER*1024 :: dataline    ! Increase for wide files
  
  IF (PRESENT(commentchar)) THEN
    commentchar0 = commentchar
  ELSE
    commentchar0 = '#'
  END IF
  
  columns = 0
  IF (PRESENT(status)) status = .FALSE.
  
  REWIND(UNIT=fid,IOSTAT=ios)
  IF (ios .LT. 0) RETURN
  
  ! Find a data line
  ios = 0
  dataline = ''
  DO WHILE ((ios .GE. 0) .AND. IsCommentLine(dataline,commentchar0))
    READ(UNIT=fid,FMT='(A)',IOSTAT=ios) dataline
  END DO
  ! No data lines?
  IF (ios .LT. 0) RETURN
  
  columns = NumberOfFields(dataline)
  
  REWIND(UNIT=fid,IOSTAT=ios)
  IF (PRESENT(status)) status = .TRUE.
  
END FUNCTION


! ----------------------------------------------------------------------
! Loads data from a formatted (text) file.
! Take an allocatable 2D array as an argument, which will be allocated
! to the correct size.  Any line beginning with a letter or symbol is
! treated as a comment line and ignored.
! Optional input arguments (at least one required):
!     file       The name of the file to be read
!     fid        The identifier for the file to be read.
!                File must already be open if given.
! Other optional input arguments:
!     commentchar  Comment character.  Default is '#'.
! Output arguments:
!     Nrow,Ncol  Number of data rows and columns; this indicates the
!                size of the allocated array
!     data       Allocatable REAL*8 array that will be allocated to
!                size data(Nrow,Ncol) and filled with data from the
!                file
! Optional output arguments:
!     success    Set to .TRUE. if the data was successfully loaded;
!                otherwise .FALSE.
! 
SUBROUTINE LoadTable(file,fid,commentchar,Nrow,Ncol,data,status)
  IMPLICIT NONE
  CHARACTER*(*), INTENT(IN), OPTIONAL :: file
  INTEGER, INTENT(IN), OPTIONAL :: fid
  CHARACTER, INTENT(IN), OPTIONAL :: commentchar
  INTEGER, INTENT(OUT), OPTIONAL :: Nrow,Ncol
  REAL*8, ALLOCATABLE, INTENT(OUT), OPTIONAL :: data(:,:)
  LOGICAL, INTENT(OUT), OPTIONAL :: status
  INTEGER :: fid0,Nrow0,Ncol0,ios,I
  LOGICAL :: status0
  CHARACTER :: commentchar0
  CHARACTER*256 :: buf
  
  ! Argument checking
  IF (.NOT. (PRESENT(file) .OR. PRESENT(fid))) THEN
    WRITE(0,*) 'ERROR: LoadTable requires a file or fid argument.'
    STOP
  END IF
  
  IF (PRESENT(commentchar)) THEN
    commentchar0 = commentchar
  ELSE
    commentchar0 = '#'
  END IF
  
  IF (PRESENT(status)) status = .FALSE.
  
  ! Open file, if necessary
  IF (PRESENT(fid)) THEN
    fid0 = fid
  ELSE
    fid0 = FreeIOUnit()
    OPEN(UNIT=fid0,FILE=file,STATUS='OLD',FORM='FORMATTED',             &
         ACTION='READ',IOSTAT=ios)
    IF (ios .NE. 0) RETURN
  END IF
  
  ! Rewind file, in case it has already been read from
  REWIND(UNIT=fid0,IOSTAT=ios)
  IF (ios .LT. 0) RETURN
  
  ! Get number of data lines
  Nrow0 = FileDataLines(fid0,commentchar0,status0)
  IF (.NOT. status0) RETURN
  IF (PRESENT(Nrow)) Nrow = Nrow0
  
  ! Get number of data columns
  Ncol0 = FileDataColumns(fid0,commentchar0,status0)
  IF (.NOT. status0) RETURN
  IF (PRESENT(Ncol)) Ncol = Ncol0
  
  ! No data array argument: nothing left to do
  IF (.NOT. PRESENT(data)) THEN
    IF (.NOT. PRESENT(fid)) THEN
      CLOSE(UNIT=fid0,IOSTAT=ios)
    END IF
    IF (PRESENT(status)) status = .TRUE.
    RETURN
  END IF
  
  ! Load data into array
  IF (ALLOCATED(data)) DEALLOCATE(data)
  ALLOCATE(data(Nrow0,Ncol0))
  I = 1
  DO WHILE (I .LE. Nrow0)
    READ(UNIT=fid0,FMT='(A)',IOSTAT=ios) buf
    ! Check if end of file or read error (ios < 0 for eof and
    ! ios > 0 for an I/O error)
    IF (ios .NE. 0) RETURN
    ! Check if comment line; if so, go back and read another line
    IF (IsCommentLine(buf)) CYCLE
    ! Parse line
    BACKSPACE(UNIT=fid0)
    READ(UNIT=fid0,FMT=*) data(I,:)
    I = I+1
  END DO
  
  ! Cleanup
  IF (.NOT. PRESENT(fid)) THEN
    CLOSE(UNIT=fid0,IOSTAT=ios)
  END IF
  
  IF (PRESENT(status)) status = .TRUE.
  
END SUBROUTINE


! ----------------------------------------------------------------------
! Loads data columns from a formatted (text) file.
! Currently allows up to 5 selectable columns to be read.
! Takes allocatable arrays as arguments, which will be allocated to
! the correct size.  Any line beginning with a letter or symbol is
! treated as a comment line and ignored.
! Optional input arguments (at least one required):
!     file       The name of the file to be read
!     fid        The identifier for the file to be read.
!                File must already be open if given.
! Other optional input arguments:
!     commentchar  Comment character.  Default is '#'.
! Output arguments:
!     N          Length of data columns; this indicates the length
!                of the allocated arrays
! Optional input arguments:
!     N1,N2,N3,N4,N5
!                File columns to read into corresponding arrays.
!                If not specified, columns 1-5 will be read.
!                A negative number counts from the end (-1 is last
!                column).
! Optional output arguments:
!     C1,C2,C3,C4,C5
!                Allocatable REAL*8 arrays that will be allocated and
!                filled with data from the columns in the file
!                corresponding to N1 - N5.
!     success    Set to .TRUE. if the data was successfully loaded;
!                otherwise .FALSE.
! 
SUBROUTINE LoadArrays(file,fid,commentchar,N,N1,C1,N2,C2,N3,C3,N4,C4,   &
                      N5,C5,status)
  IMPLICIT NONE
  CHARACTER*(*), INTENT(IN), OPTIONAL :: file
  INTEGER, INTENT(IN), OPTIONAL :: fid
  CHARACTER, INTENT(IN), OPTIONAL :: commentchar
  INTEGER, INTENT(OUT), OPTIONAL :: N
  INTEGER, INTENT(IN), OPTIONAL :: N1,N2,N3,N4,N5
  REAL*8, ALLOCATABLE, INTENT(OUT), OPTIONAL :: C1(:),C2(:),C3(:),      &
                                                C4(:),C5(:)
  LOGICAL, INTENT(OUT), OPTIONAL :: status
  INTEGER :: fid0,N0,Ncol0,ios,I,J,J1,J2,J3,J4,J5
  LOGICAL :: status0
  REAL*8, ALLOCATABLE :: data(:,:)
  CHARACTER :: commentchar0
  CHARACTER*256 :: buf
  
  ! Argument checking
  IF (.NOT. (PRESENT(file) .OR. PRESENT(fid))) THEN
    WRITE(0,*) 'ERROR: LoadArrays requires a file or fid argument.'
    STOP
  END IF

  IF (PRESENT(commentchar)) THEN
    commentchar0 = commentchar
  ELSE
    commentchar0 = '#'
  END IF
  
  IF (PRESENT(status)) status = .FALSE.
  
  ! Open file, if necessary
  IF (PRESENT(fid)) THEN
    fid0 = fid
  ELSE
    fid0 = FreeIOUnit()
    OPEN(UNIT=fid0,FILE=file,STATUS='OLD',FORM='FORMATTED',             &
         ACTION='READ',IOSTAT=ios)
    IF (ios .NE. 0) RETURN
  END IF
  
  ! Rewind file, in case it has already been read from
  REWIND(UNIT=fid0,IOSTAT=ios)
  IF (ios .LT. 0) RETURN
  
  ! Get number of data lines
  N0 = FileDataLines(fid0,commentchar0,status0)
  IF (.NOT. status0) RETURN
  IF (PRESENT(N)) N = N0
  
  ! Get number of data columns
  Ncol0 = FileDataColumns(fid0,commentchar0,status0)
  IF (.NOT. status0) RETURN
  
  ! Find how many columns should be read
  ! Following sets indices to 0, default, or appropriate
  ! optional argument
  IF (UpdateColNum(PRESENT(C1),PRESENT(N1),1,J1)) J1 = N1
  IF (PRESENT(C1)) CALL UpdateColNum2(Ncol0,J1,status0)
  IF (UpdateColNum(PRESENT(C2),PRESENT(N2),2,J2)) J2 = N2
  IF (PRESENT(C2)) CALL UpdateColNum2(Ncol0,J2,status0)
  IF (UpdateColNum(PRESENT(C3),PRESENT(N3),3,J3)) J3 = N3
  IF (PRESENT(C3)) CALL UpdateColNum2(Ncol0,J3,status0)
  IF (UpdateColNum(PRESENT(C4),PRESENT(N4),4,J4)) J4 = N4
  IF (PRESENT(C4)) CALL UpdateColNum2(Ncol0,J4,status0)
  IF (UpdateColNum(PRESENT(C5),PRESENT(N5),5,J5)) J5 = N5
  IF (PRESENT(C5)) CALL UpdateColNum2(Ncol0,J5,status0)
  IF (.NOT. status0) RETURN
  
  ! Only read in to last necessary column
  Ncol0 = MAX(J1,J2,J3,J4,J5)
  
  ! Load data into temporary array
  ALLOCATE(data(N0,Ncol0))
  I = 1
  DO WHILE (I .LE. N0)
    READ(UNIT=fid0,FMT='(A)',IOSTAT=ios) buf
    ! Check if end of file or read error (ios < 0 for eof and
    ! ios > 0 for an I/O error)
    IF (ios .NE. 0) RETURN
    ! Check if comment line; if so, go back and read another line
    IF (IsCommentLine(buf)) CYCLE
    ! Parse line
    BACKSPACE(UNIT=fid0)
    READ(UNIT=fid0,FMT=*) data(I,:)
    I = I+1
  END DO
  
  ! Allocate and fill arrays
  IF (PRESENT(C1)) CALL UpdateArray(J1,C1)
  IF (PRESENT(C2)) CALL UpdateArray(J2,C2)
  IF (PRESENT(C3)) CALL UpdateArray(J3,C3)
  IF (PRESENT(C4)) CALL UpdateArray(J4,C4)
  IF (PRESENT(C5)) CALL UpdateArray(J5,C5)
  
  ! Cleanup
  DEALLOCATE(data)
  IF (.NOT. PRESENT(fid)) THEN
    CLOSE(UNIT=fid0,IOSTAT=ios)
  END IF
  
  IF (PRESENT(status)) status = .TRUE.
  
  
  CONTAINS
    
  ! ------------------------------------------
  ! Utility to determine column number
  ! Will set column number to 0 if no column will be read or
  ! the default defCN if a column will be read.  If givenCN
  ! is .TRUE. (indicating if the optional NJ argument was given)
  ! and the column is to be read, this function will return
  ! .TRUE. indicating CN should be read from NJ instead.
  FUNCTION UpdateColNum(usecol,givenCN,defCN,CN)
    IMPLICIT NONE
    LOGICAL :: UpdateColNum
    LOGICAL, INTENT(IN) :: usecol,givenCN
    INTEGER, INTENT(IN) :: defCN
    INTEGER, INTENT(OUT) :: CN
    UpdateColNum = .FALSE.
    IF (usecol) THEN
      CN = defCN
      IF (givenCN) UpdateColNum = .TRUE.
    ELSE
      CN = 0
    END IF
  END FUNCTION UpdateColNum
  
  ! ------------------------------------------
  ! Utility to update column number
  ! Negative indices (indicating column number from end) are
  ! converted to true indices.  Validity of the column number
  ! is checked.
  SUBROUTINE UpdateColNum2(Ncol,CN,status0)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: Ncol
    INTEGER, INTENT(INOUT) :: CN
    LOGICAL, INTENT(INOUT) :: status0
    IF (Ncol .LT. 0) THEN
      CN = Ncol + 1 - CN
    END IF
    status0 = status0 .AND. (CN .GT. 0) .AND. (CN .LE. Ncol)
  END SUBROUTINE UpdateColNum2
  
  ! ------------------------------------------
  ! Allocates and fills an array arr with data from column J
  SUBROUTINE UpdateArray(J,arr)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: J
    REAL*8, ALLOCATABLE, INTENT(OUT) :: arr(:)
    IF (ALLOCATED(arr)) DEALLOCATE(arr)
    ALLOCATE(arr(1:N0))
    arr = data(:,J)
  END SUBROUTINE UpdateArray
  
END SUBROUTINE


END MODULE
