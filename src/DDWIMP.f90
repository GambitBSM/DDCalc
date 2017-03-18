MODULE DDWIMP

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! DDWIMP
!    Routines for initializing, modifying, or viewing WIMP properties,
!    mainly the mass and couplings.
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

USE DDTypes
USE DDUtils
USE DDCommandLine
USE DDCouplings

IMPLICIT NONE
PRIVATE

! WIMP mass & couplings routines
PUBLIC :: DDCalc_GetWIMP,DDCalc_SetWIMP,     &
          DDCalc_InitWIMP,DDCalc_InitWIMPCommandLine, &
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
INTERFACE DDCalc_InitWIMPCommandLine
  MODULE PROCEDURE InitWIMPCommandLine
END INTERFACE

CONTAINS


! ----------------------------------------------------------------------
! Get various WIMP quantities.
! 
! Optional output arguments:
!   m           WIMP mass [GeV].
!   GpSI        Spin-independent WIMP-proton coupling [GeV^-2].
!   GnSI        Spin-independent WIMP-neutron coupling [GeV^-2].
!   GpSD        Spin-dependent WIMP-proton coupling [GeV^-2].
!   GnSD        Spin-dependent WIMP-neutron coupling [GeV^-2].
!   GpSI0       Reference spin-independent WIMP-proton coupling [GeV^-2]
!               corresponding to \sigma_{SI,p} = 1 pb.
!   GnSI0       Reference spin-independent WIMP-neutron coupling [GeV^-2]
!               corresponding to \sigma_{SI,n} = 1 pb.
!   GpSD0       Reference spin-dependent WIMP-proton coupling [GeV^-2]
!               corresponding to \sigma_{SD,p} = 1 pb.
!   GnSD0       Reference spin-dependent WIMP-neutron coupling [GeV^-2]
!               corresponding to \sigma_{SD,n} = 1 pb.
!   fp          Spin-independent WIMP-proton coupling [GeV^-2].
!               Related by GpSI = 2 fp.
!   fn          Spin-independent WIMP-neutron coupling [GeV^-2].
!               Related by GnSI = 2 fn.
!   ap          Spin-dependent WIMP-proton coupling [unitless].
!               Related by GpSD = 2\sqrt{2} G_F ap.
!   an          Spin-dependent WIMP-neutron coupling [unitless].
!               Related by GnSD = 2\sqrt{2} G_F an.
!   sigmapSI    Spin-independent WIMP-proton cross-section [pb].
!   sigmanSI    Spin-independent WIMP-neutron cross-section [pb].
!   sigmapSD    Spin-dependent WIMP-proton cross-section [pb].
!   sigmanSD    Spin-dependent WIMP-neutron cross-section [pb].
! 
SUBROUTINE GetWIMP(WIMP, m,GpSI,GnSI,GpSD,GnSD,GpSI0,GnSI0,GpSD0,GnSD0,       &
                   fp,fn,ap,an,sigmapSI,sigmanSI,sigmapSD,sigmanSD)

  IMPLICIT NONE
  TYPE(WIMPStruct), INTENT(IN) :: WIMP
  REAL*8, INTENT(OUT), OPTIONAL :: m,GpSI,GnSI,GpSD,GnSD,GpSI0,GnSI0,GpSD0,GnSD0, &
           fp,fn,ap,an,sigmapSI,sigmanSI,sigmapSD,sigmanSD
  IF (PRESENT(m))     m     = WIMP%m
  IF (PRESENT(GpSI))  GpSI  = WIMP%GpSI
  IF (PRESENT(GnSI))  GnSI  = WIMP%GnSI
  IF (PRESENT(GpSD))  GpSD  = WIMP%GpSD
  IF (PRESENT(GnSD))  GnSD  = WIMP%GnSD
  IF (PRESENT(GpSI0)) GpSI0 = WIMP%GpSI0
  IF (PRESENT(GnSI0)) GnSI0 = WIMP%GnSI0
  IF (PRESENT(GpSD0)) GpSD0 = WIMP%GpSD0
  IF (PRESENT(GnSD0)) GnSD0 = WIMP%GnSD0
  IF (PRESENT(fp))    fp    = GToF(WIMP%GpSI)
  IF (PRESENT(fn))    fn    = GToF(WIMP%GnSI)
  IF (PRESENT(ap))    ap    = GToA(WIMP%GpSD)
  IF (PRESENT(an))    an    = GToA(WIMP%GnSD)
  IF (PRESENT(sigmapSI)) sigmapSI = GpToSigmapSI(WIMP%m,WIMP%GpSI)
  IF (PRESENT(sigmanSI)) sigmanSI = GnToSigmanSI(WIMP%m,WIMP%GnSI)
  IF (PRESENT(sigmapSD)) sigmapSD = GpToSigmapSD(WIMP%m,WIMP%GpSD)
  IF (PRESENT(sigmanSD)) sigmanSD = GnToSigmanSD(WIMP%m,WIMP%GnSD)
END SUBROUTINE


! ----------------------------------------------------------------------
! Set various WIMP quantities.
! 
! Optional input arguments:
!   m           WIMP mass [GeV].
!   GpSI        Spin-independent WIMP-proton coupling [GeV^-2].
!   GnSI        Spin-independent WIMP-neutron coupling [GeV^-2].
!   GpSD        Spin-dependent WIMP-proton coupling [GeV^-2].
!   GnSD        Spin-dependent WIMP-neutron coupling [GeV^-2].
!   fp          Spin-independent WIMP-proton coupling [GeV^-2].
!               Related by GpSI = 2 fp.
!   fn          Spin-independent WIMP-neutron coupling [GeV^-2].
!               Related by GnSI = 2 fn.
!   ap          Spin-dependent WIMP-proton coupling [unitless].
!               Related by GpSD = 2\sqrt{2} G_F ap.
!   an          Spin-dependent WIMP-neutron coupling [unitless].
!               Related by GnSD = 2\sqrt{2} G_F an.
! Optional cross-section arguments (give negative value to set
! corresponding coupling negative):
!   sigmapSI    Spin-independent WIMP-proton cross-section [pb].
!   sigmanSI    Spin-independent WIMP-neutron cross-section [pb].
!   sigmapSD    Spin-dependent WIMP-proton cross-section [pb].
!   sigmanSD    Spin-dependent WIMP-neutron cross-section [pb].
!   sigmaSI     Sets both sigmapSI and sigmanSI to the given value [pb].
!   sigmaSD     Sets both sigmapSD and sigmanSD to the given value [pb].
! 
SUBROUTINE SetWIMP(WIMP,m,GpSI,GnSI,GpSD,GnSD,fp,fn,ap,an,                   &
                   sigmapSI,sigmanSI,sigmapSD,sigmanSD,sigmaSI,sigmaSD)

  IMPLICIT NONE
  TYPE(WIMPStruct), INTENT(INOUT) :: WIMP
  REAL*8, INTENT(IN), OPTIONAL :: m,GpSI,GnSI,GpSD,GnSD,fp,fn,ap,an,    &
           sigmapSI,sigmanSI,sigmapSD,sigmanSD,sigmaSI,sigmaSD
  IF (PRESENT(m)) THEN
     WIMP%m = MAX(m,SQRT(TINY(1d0)))
     WIMP%GpSI0 = SigmapSIToGp(WIMP%m,1d0)
     WIMP%GnSI0 = SigmanSIToGn(WIMP%m,1d0)
     WIMP%GpSD0 = SigmapSDToGp(WIMP%m,1d0)
     WIMP%GnSD0 = SigmanSDToGn(WIMP%m,1d0)
  END IF
  IF (PRESENT(GpSI))  WIMP%GpSI = GpSI
  IF (PRESENT(GnSI))  WIMP%GnSI = GnSI
  IF (PRESENT(GpSD))  WIMP%GpSD = GpSD
  IF (PRESENT(GnSD))  WIMP%GnSD = GnSD
  IF (PRESENT(fp))    WIMP%GpSI = FToG(fp)
  IF (PRESENT(fn))    WIMP%GnSI = FToG(fn)
  IF (PRESENT(ap))    WIMP%GpSD = AToG(ap)
  IF (PRESENT(an))    WIMP%GnSD = AToG(an)
  IF (PRESENT(sigmapSI)) WIMP%GpSI = SigmapSIToGp(WIMP%m,sigmapSI)
  IF (PRESENT(sigmanSI)) WIMP%GnSI = SigmanSIToGn(WIMP%m,sigmanSI)
  IF (PRESENT(sigmapSD)) WIMP%GpSD = SigmapSDToGp(WIMP%m,sigmapSD)
  IF (PRESENT(sigmanSD)) WIMP%GnSD = SigmanSDToGn(WIMP%m,sigmanSD)
  IF (PRESENT(sigmaSI)) THEN
    WIMP%GpSI = SigmapSIToGp(WIMP%m,sigmaSI)
    WIMP%GnSI = SigmanSIToGn(WIMP%m,sigmaSI)
  END IF
  IF (PRESENT(sigmaSD)) THEN
    WIMP%GpSD = SigmapSDToGp(WIMP%m,sigmaSD)
    WIMP%GnSD = SigmanSDToGn(WIMP%m,sigmaSD)
  END IF

END SUBROUTINE


!-----------------------------------------------------------------------
! Initializes WIMP.
! Simply sets some default values for the WIMP parameters.
! 
FUNCTION InitWIMP() RESULT(WIMP)

  IMPLICIT NONE
  TYPE(WIMPStruct) :: WIMP

  ! Default mass of 100 GeV
  CALL SetWIMP(WIMP,m=100d0)
  
  ! Default cross-sections of 1 pb.
  WIMP%GpSI = WIMP%GpSI0
  WIMP%GnSI = WIMP%GnSI0
  WIMP%GpSD = WIMP%GpSD0
  WIMP%GnSD = WIMP%GnSD0
  
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
! Initializes WIMP from command-line parameters.
! 
! Possible options:
!   --m=<value>          ! WIMP mass [GeV]
!   --GpSI=<value>       ! Spin-independent WIMP-proton coupling [GeV^-2].
!   --GnSI=<value>       ! Spin-independent WIMP-neutron coupling [GeV^-2].
!   --GpSD=<value>       ! Spin-dependent WIMP-proton coupling [GeV^-2].
!   --GnSD=<value>       ! Spin-dependent WIMP-neutron coupling [GeV^-2].
!   --fp=<value>         ! Spin-independent WIMP-proton coupling [GeV^-2].
!                        ! Related by GpSI = 2 fp.
!   --fn=<value>         ! Spin-independent WIMP-neutron coupling [GeV^-2].
!                        ! Related by GnSI = 2 fn.
!   --ap=<value>         ! Spin-dependent WIMP-proton coupling [unitless].
!                        ! Related by GpSD = 2\sqrt{2} G_F ap.
!   --an=<value>         ! Spin-dependent WIMP-neutron coupling [unitless].
!                        ! Related by GnSD = 2\sqrt{2} G_F an.
! Cross-section options may be given as negative values to indicate the
! corresponding coupling should be negative:
!   --sigmapSI=<value>   ! Spin-independent WIMP-proton cross-section [pb].
!   --sigmanSI=<value>   ! Spin-independent WIMP-neutron cross-section [pb].
!   --sigmapSD=<value>   ! Spin-dependent WIMP-proton cross-section [pb].
!   --sigmanSD=<value>   ! Spin-dependent WIMP-neutron cross-section [pb].
!   --sigmaSI=<value>    ! Sets both sigmapSI and sigmanSI to the given value [pb].
!   --sigmaSD=<value>    ! Sets both sigmapSD and sigmanSD to the given value [pb].
! 
FUNCTION InitWIMPCommandLine(Arguments) RESULT(WIMP)

  IMPLICIT NONE
  TYPE(ArgumentStruct), INTENT(IN) :: Arguments
  TYPE(WIMPStruct) :: WIMP
  LOGICAL :: status
  REAL*8 :: x
  
  ! Process mass: default value of 100 GeV
  CALL SetWIMP(WIMP,m=100d0)
  IF (GetLongArgReal('m',x)) CALL SetWIMP(WIMP,m=x)
  IF (Arguments%Nparameters .GE. 1) THEN
    x = Arguments%values(1)
    IF (x .GT. 0d0) CALL SetWIMP(WIMP,m=x)
  END IF
  
  ! Process couplings: defaults of 1 pb.
  WIMP%GpSI = WIMP%GpSI0
  WIMP%GnSI = WIMP%GnSI0
  WIMP%GpSD = WIMP%GpSD0
  WIMP%GnSD = WIMP%GnSD0
  IF (GetLongArgReal('GpSI',x)) CALL SetWIMP(WIMP,GpSI=x)
  IF (GetLongArgReal('GnSI',x)) CALL SetWIMP(WIMP,GnSI=x)
  IF (GetLongArgReal('GpSD',x)) CALL SetWIMP(WIMP,GpSD=x)
  IF (GetLongArgReal('GnSD',x)) CALL SetWIMP(WIMP,GnSD=x)
  IF (GetLongArgReal('fp',x))   CALL SetWIMP(WIMP,fp=x)
  IF (GetLongArgReal('fn',x))   CALL SetWIMP(WIMP,fn=x)
  IF (GetLongArgReal('ap',x))   CALL SetWIMP(WIMP,ap=x)
  IF (GetLongArgReal('an',x))   CALL SetWIMP(WIMP,an=x)
  IF (GetLongArgReal('sigmapSI',x)) CALL SetWIMP(WIMP,sigmapSI=x)
  IF (GetLongArgReal('sigmanSI',x)) CALL SetWIMP(WIMP,sigmanSI=x)
  IF (GetLongArgReal('sigmapSD',x)) CALL SetWIMP(WIMP,sigmapSD=x)
  IF (GetLongArgReal('sigmanSD',x)) CALL SetWIMP(WIMP,sigmanSD=x)
  IF (GetLongArgReal('sigmaSI',x))  CALL SetWIMP(WIMP,sigmaSI=x)
  IF (GetLongArgReal('sigmaSD',x))  CALL SetWIMP(WIMP,sigmaSD=x)
  
  ! Process command-line arguments (if more than just mass)
  IF (Arguments%Nparameters .GE. 1) THEN
    status = ParseWIMPParameters(Arguments%Nparameters,Arguments%values,WIMP)
  END IF
  
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
  
  ! Check for bad cases (e.g. non-positive mass)
  status = .FALSE.
  IF (N .LT. 1) RETURN
  IF (p(1) .LE. 0d0) RETURN
  
  status = .TRUE.
  
  ! Meaning of parameters depends on number of parameters,
  ! but first parameter is always mass
  CALL SetWIMP(WIMP,m=p(1))
  SELECT CASE (N)
  CASE (1)
    ! Form: m
    ! Set WIMP couplings to 1 pb
    WIMP%GpSI = WIMP%GpSI0
    WIMP%GnSI = WIMP%GnSI0
    WIMP%GpSD = WIMP%GpSD0
    WIMP%GnSD = WIMP%GnSD0
  CASE (2)
    ! Form: m sigmaSI
    ! Set SD couplings to zero
    CALL SetWIMP(WIMP,sigmaSI=p(2),sigmaSD=0d0)
  CASE (3)
    ! Form: m sigmaSI sigmaSD
    CALL SetWIMP(WIMP,sigmaSI=p(2),sigmaSD=p(3))
  CASE (4)
    ! Form: m sigmaSI sigmapSD sigmanSD
    CALL SetWIMP(WIMP,sigmaSI=p(2),sigmapSD=p(3),sigmanSD=p(4))
  CASE (5)
    !! Form: m sigmapSI sigmanSI sigmapSD sigmanSD
    CALL SetWIMP(WIMP,sigmapSI=p(2),sigmanSI=p(3),sigmapSD=p(4),sigmanSD=p(5))
    ! Form: m GpSI GnSI GpSD GnSD
    !CALL SetWIMP(WIMP,GpSI=p(2),GnSI=p(3),GpSD=p(4),GnSD=p(5))
  CASE (6:)
    status = .FALSE.
  END SELECT
  
END FUNCTION


END MODULE
