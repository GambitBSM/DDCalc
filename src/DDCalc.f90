!#######################################################################
MODULE DDCalc
!#######################################################################

!#######################################################################
! DIRECT DETECTION RATES AND LIKELIHOODS
! Routines to calculate the dark matter direct detection event rates
! and corresponding likelihoods/exclusion levels.
! 
!   Created by Chris Savage
!   University of Utah   (2013 - 2014)
!   Nordita              (2014 - 2015)
! 
!   With contributions from:
!     Andre Scaffidi            University of Adelaide (2014)
!     Lauren Hsu                FermiLab (2015)
!     Pat Scott                 Imperial College London (2016)
!     Felix Kahlhoefer          RWTH Aachen (2018)
!     Sebastian Wild		DESY (2018)
! 
! 
! 
!------========<<<<<<<< INTERFACE DESCRIPTION >>>>>>>>========----------
!
! TYPES: WIMPStruct, HaloStruct, DetectorStruct
!===============================================
! These three types are the bedrock of DDCalc.  (Almost) every calculation
! needs to be provided with an instance of each of these to do its job.
!
!
! WIMP MODEL
!============
! WIMP model initialization:
!     TYPE(WIMPStruct) FUNCTION DDCalc_InitWIMP()
!
! WIMP parameter setting:
!     SUBROUTINE DDCalc_SetWIMP_mfa(WIMP,m,fp,fn,ap,an)
!     SUBROUTINE DDCalc_SetWIMP_mG(WIMP,m,GpSI,GnSI,GpSD,GnSD)
!   where m is the WIMP mass [GeV], f is spin-independent (SI) WIMP-
!   nucleon coupling [GeV^-2], a is the spin-dependent (SD) WIMP-nucleon
!   coupling [unitless], G are the effective 4-fermion vertex couplings,
!   related by:
!     GpSI = 2 fp        GpSD = 2\sqrt{2} G_F ap
!     GnSI = 2 fn        GnSD = 2\sqrt{2} G_F an
!
!     SUBROUTINE DDCalc_SetWIMP_msigma(WIMP,m,sigmaSIp,sigmaSIn,sigmaSDp,sigmaSDn)
!   where m is the WIMP mass [GeV], sigma is the WIMP-nucleon scattering cross-section [pb].
!   Negative cross-sections indicate the corresponding coupling should
!   be negative.  In all cases, 'p' refers to proton and 'n' to neutron.
!
! WIMP parameter retrieval:
!     SUBROUTINE DDCalc_GetWIMP_mfa(WIMP,m,fp,fn,ap,an)
!     SUBROUTINE DDCalc_GetWIMP_mG(WIMP,m,fsp,fsn,app,apn)
!     SUBROUTINE DDCalc_GetWIMP_msigma(WIMP,m,sigmaSIp,sigmaSIn,sigmaSDp,sigmaSDn)
!   Same as SetWIMP above, but retrieves current WIMP parameters.
!
! For advanced setters and getters for WIMP properties, see:
! DDCalc_SetWIMP() and DDCalc_GetWIMP() in DDWIMP.f90.
!
!
! HALO MODEL
!============
! Halo model initialization:
!     TYPE(HaloStruct) FUNCTION DDCalc_InitHalo()
!
! Standard Halo Model (SHM) settings:
!     SUBROUTINE DDCalc_SetSHM(Halo,rho,vrot,v0,vesc)
! where the arguments are the local dark matter density [GeV/cm^3],
! the local disk rotation speed [km/s], the most probable dark matter
! speed in the galactic rest frame [km/s], and the galactic escape
! speed [km/s].  This routine need only be called if non-default
! values are desired (0.4, 235, 235, and 550, respectively).
!
! For advanced setters and getters for Halo properties, see:
! DDCalc_SetHalo() and DDCalc_GetHalo() in DDHalo.f90.
!
!
! EXPERIMENTAL DETECTOR AND ANALYSIS
!====================================
! Detector/Analysis initialization:
!     TYPE(DetectorStruct) FUNCTION DDCalc_InitDetector()
! Initialise an object carrying the information about the experimental
! analysis to consider. Returns a detector structure containing the
! analysis details.  Non-default analyses can be obtained with
! specific functions:
!     TYPE(DetectorStruct) FUNCTION [Analysis]_Init()
! See DDExperiments.f90 and analyses/[Analysis].f90 for more details.
! Note that these functions are only available directly by USE-ing the
! DDExperiments module in the calling program, and will not be made
! available simply by USE-ing the DDCalc module.
!
! To run any of the following routines with a particular experiment, the
! corresponding Init routine must be called only once, and the resulting
! Detector object passed to the analysis routines to be run on it.
!
! Set dectector minimum recoil energy:
!     SUBROUTINE DDCalc_SetDetectorEmin(Detector,Emin)
! Emin in keV (initially set to 0 keV). Note the efficiency curves
! already account for detector and analysis thresholds regardless of
! this setting, so setting this to 0 keV (the default behavior when
! initialization is performed) does not imply that very low energy
! recoils actually contribute to the signal.
!
! Do rate calculations:
!     SUBROUTINE DDCalc_CalcRates(Detector,WIMP,Halo)
! Performs the rate calculations used for likelihoods/confidence
! intervals.  Results are saved internally in the Detector analysis
! object, and can be accessed using the following routines.
!
! Get results of rate calculations:
!
! 1. Number of observed events in the analysis:
!     INTEGER FUNCTION DDCalc_Events(Detector)
!
! 2. Average expected number of background events in the analysis:
!     REAL*8 FUNCTION DDCalc_Background(Detector)
!
! 3. Average expected number of signal events in the analysis:
!     REAL*8 FUNCTION DDCalc_Signal(Detector)
!
! 4. Log-likelihood:
!     REAL*8 FUNCTION DDCalc_LogLikelihood(Detector)
!    Calculates the likelihood or the p-value depending on the StatisticFlag
!    of the detector (see DDStats for details).
!
! 5. Factor by which the WIMP cross-sections must be multiplied to
!    achieve a given p-value:
!     REAL*8 FUNCTION DDCalc_ScaleToPValue(lnp)
!    Calculates the factor x by which the cross-sections must be scaled
!    (sigma -> x*sigma) to achieve the desired p-value (given as log(p)).f
!
! For experiments with binned analyses, one can also use the following functions
!
! 6. Number of bins in the analysis:
!     INTEGER FUNCTION DDCalc_Bins(Detector)
!
! 7. Number of observed events in the bin with index ibin:
!     INTEGER FUNCTION DDCalc_BinEvents(Detector, ibin)
!
! 8. Average expected number of background events in the bin with index ibin:
!     REAL*8 FUNCTION DDCalc_BinBackground(Detector, ibin)
!
! 9. Average expected number of signal events in the bin with index ibin:
!     REAL*8 FUNCTION DDCalc_BinSignal(Detector, ibin)

!
! C/C++ INTERFACE 
!=================
! For ease of use in linking to these routines from C/C++ code, a second
! (wrapper) version of each of the interface routines described above
! is defined with a 'C_' prefix.  These use C-compatible types only,
! and are given explicitly-specified symbol names, to get around
! name-mangling inconsistencies between different compilers.  These
! functions work just like the ones above, but do neither accept nor
! return WIMPStructs, HaloStructs nor DetectorStructs directly.  Instead,
! they return and accept integers corresponding to entries in an internal
! array of Fortran objects held in trust for the C/C++ calling program
! by DDCalc. The functions C_DDUtils_FreeWIMPs, C_DDUtils_FreeHalos,
! C_DDUtils_FreeDetectors and C_DDUtils_FreeAll can be used to delete
! the objects held internally.
!
!#######################################################################

USE DDExperiments
USE DDConstants
USE DDTypes
USE DDNumerical
USE DDWIMP
USE DDDetectors
USE DDRates
USE DDStats
USE DDHalo
USE DDCouplings
USE DDNREffectiveTheory

IMPLICIT NONE
PRIVATE


!############## Basic Interface ####################

! DDCalc types for interfacing with the library
PUBLIC :: WIMPStruct
PUBLIC :: HaloStruct
PUBLIC :: DetectorStruct

! Default initialisation for each type
PUBLIC :: DDCalc_InitWIMP
PUBLIC :: DDCalc_InitHalo
PUBLIC :: DDCalc_InitDetector
! C/C++ Wrappers
PUBLIC :: C_DDCalc_InitWIMP
PUBLIC :: C_DDCalc_InitHalo
PUBLIC :: C_DDCalc_InitDetector

! Initialisation of implemented expermimental analyses
! (see DDExperiments.f90 and analyses/[Analysis].f90.)
!PUBLIC :: [Analysis]_Init
! C/C++ Wrappers
!PUBLIC :: C_[Analysis]_Init

! Simple setter and getter interfaces
PUBLIC :: DDCalc_SetSHM
PUBLIC :: DDCalc_HaloFromFile
PUBLIC :: DDCalc_SetWIMP_mfa
PUBLIC :: DDCalc_GetWIMP_mfa
PUBLIC :: DDCalc_SetWIMP_mG
PUBLIC :: DDCalc_GetWIMP_mG
PUBLIC :: DDCalc_SetWIMP_msigma
PUBLIC :: DDCalc_GetWIMP_msigma
PUBLIC :: DDCalc_SetWIMP_longrange
PUBLIC :: DDCalc_GetWIMP_longrange
PUBLIC :: DDCalc_SetWIMP_NREffectiveTheory
PUBLIC :: DDCalc_SetWIMP_NREFT_CPT
PUBLIC :: DDCalc_SetNRCoefficient
PUBLIC :: DDCalc_GetNRCoefficient
PUBLIC :: DDCalc_SetDetectorEmin
! C/C++ Wrappers
PUBLIC :: C_DDCalc_SetSHM
PUBLIC :: C_DDCalc_HaloFromFile
PUBLIC :: C_DDCalc_SetWIMP_mfa
PUBLIC :: C_DDCalc_GetWIMP_mfa
PUBLIC :: C_DDCalc_SetWIMP_mG
PUBLIC :: C_DDCalc_GetWIMP_mG
PUBLIC :: C_DDCalc_SetWIMP_msigma
PUBLIC :: C_DDCalc_GetWIMP_msigma
PUBLIC :: C_DDCalc_SetWIMP_longrange
PUBLIC :: C_DDCalc_GetWIMP_longrange
PUBLIC :: C_DDCalc_SetWIMP_NREffectiveTheory
PUBLIC :: C_DDCalc_SetWIMP_NREFT_CPT
PUBLIC :: C_DDCalc_SetNRCoefficient
PUBLIC :: C_DDCalc_GetNRCoefficient
PUBLIC :: C_DDCalc_SetDetectorEmin

! Rate calculation
PUBLIC :: DDCalc_CalcRates
! C/C++ Wrapper
PUBLIC :: C_DDCalc_CalcRates

! Result inspection
PUBLIC :: DDCalc_Events
PUBLIC :: DDCalc_Background
PUBLIC :: DDCalc_Signal
PUBLIC :: DDCalc_LogLikelihood
PUBLIC :: DDCalc_ScaleToPValue
! C/C++ Wrappers
PUBLIC :: C_DDCalc_Events
PUBLIC :: C_DDCalc_Background
PUBLIC :: C_DDCalc_Signal
PUBLIC :: C_DDCalc_LogLikelihood
PUBLIC :: C_DDCalc_ScaleToPValue

! Memory cleanup for C/C++ interface
PUBLIC :: C_DDCalc_FreeWIMPs
PUBLIC :: C_DDCalc_FreeHalos
PUBLIC :: C_DDCalc_FreeDetectors
PUBLIC :: C_DDCalc_FreeAll

! More advanced setter and getter interfaces
PUBLIC :: DDCalc_GetWIMP
PUBLIC :: DDCalc_SetWIMP
PUBLIC :: DDCalc_GetHalo
PUBLIC :: DDCalc_SetHalo
PUBLIC :: DDCalc_GetDetector
PUBLIC :: DDCalc_SetDetector

CONTAINS


!#######################################################################
! ROUTINES
!#######################################################################


!=======================================================================
! SIMPLE INTERFACE ROUTINES
! Basic versions of routines that are required for using this module
! externally.  These are meant to allow for a simpler interface to
! this module; other routines are more robust and provide more
! capabilities.
!=======================================================================


!-----------------------------------------------------------------------
! Sets the dark matter halo to the Standard Halo Model (SHM) with the
! given parameters.  Need only be run if non-default parameters are
! to be used.
! 
! For more detailed halo settings, see:
!   SetHalo() [interface name: DDCalc_SetHalo]
! 
! Input arguments:
!   rho         Local dark matter density [GeV/cm^3].  Default is
!               0.4 GeV/cm^3.
!   vrot        Local galactic disk rotation speed [km/s].  Default is
!               235 km/s.
!   v0          Most probable speed [km/s] in the galactic rest frame.
!               For the conventional isothermal sphere, this should be
!               the same as vrot.  Default is 235 km/s.
!   vesc        Galactic escape speed [km/s] in the galactic rest
!               frame.  Default is 550 km/s.
! 
SUBROUTINE DDCalc_SetSHM(Halo,rho,vrot,v0,vesc)
  IMPLICIT NONE
  TYPE(HaloStruct), INTENT(INOUT) :: Halo
  REAL*8, INTENT(IN) :: rho,vrot,v0,vesc
  CALL DDCalc_SetHalo(Halo,rho=rho,vrot=vrot,v0=v0,vesc=vesc)
END SUBROUTINE

! C++ interface wrapper
SUBROUTINE C_DDCalc_SetSHM(HaloIndex,rho,vrot,v0,vesc) &
           BIND(C,NAME='C_DDCalc_ddcalc_setshm')
  USE ISO_C_BINDING, only: C_DOUBLE, C_INT
  IMPLICIT NONE
  REAL(KIND=C_DOUBLE), INTENT(IN) :: rho,vrot,v0,vesc
  INTEGER(KIND=C_INT), INTENT(IN) :: HaloIndex
  IF (.NOT. ASSOCIATED(Halos(HaloIndex)%p)) stop 'Invalid halo index given to C_DDCalc_SetSHM'
  CALL DDCalc_SetSHM(Halos(HaloIndex)%p,rho=REAL(rho,KIND=8),vrot=REAL(vrot,KIND=8), &
                     v0=REAL(v0,KIND=8),vesc=REAL(vesc,KIND=8))
END SUBROUTINE


!-----------------------------------------------------------------------
! Tells DDCalc to read halo model from file
! 
! For more detailed halo settings, see:
!   SetHalo() [interface name: DDCalc_SetHalo]
! 
! Input arguments:
!   rho         Local dark matter density [GeV/cm^3].
!   Nvmin       Number of entries in external file.
!   gcol        Column to use for velocity integral.
! 
SUBROUTINE DDCalc_HaloFromFile(Halo,rho,Nvmin,gcol)
  IMPLICIT NONE
  TYPE(HaloStruct), INTENT(INOUT) :: Halo
  REAL*8, INTENT(IN) :: rho
  INTEGER, INTENT(IN) :: Nvmin,gcol
  CALL DDCalc_SetHalo(Halo,rho=rho,g_file='halo.txt',g_column=gcol,Nvmin=Nvmin)
END SUBROUTINE

! C++ interface wrapper
SUBROUTINE C_DDCalc_HaloFromFile(HaloIndex,rho,Nvmin,gcol) &
           BIND(C,NAME='C_DDCalc_ddcalc_halofromfile')
  USE ISO_C_BINDING, only: C_DOUBLE, C_INT
  IMPLICIT NONE
  REAL(KIND=C_DOUBLE), INTENT(IN) :: rho
  INTEGER(KIND=C_INT), INTENT(IN) :: HaloIndex, Nvmin, gcol
  IF (.NOT. ASSOCIATED(Halos(HaloIndex)%p)) stop 'Invalid halo index given to C_DDCalc_HaloFromFile'
  CALL DDCalc_HaloFromFile(Halos(HaloIndex)%p,rho=REAL(rho,KIND=8),Nvmin=INT(Nvmin,KIND=C_INT),&
                       gcol=INT(gcol,KIND=C_INT))
END SUBROUTINE


!-----------------------------------------------------------------------
! Sets/gets the WIMP mass and couplings for standard spin-independent
! and spin-dependent interactions.  Couplings are specified via
! the commonly used fp/fn (spin-independent) and ap/an (spin-dependent)
! normalizations.
! 
! For more detailed WIMP settings, see:
!   SetWIMP() [interface name: DDCalc_SetWIMP]
!   GetWIMP() [interface name: DDCalc_GetWIMP]
! 
! Input/output arguments:
!   m           WIMP mass [GeV].
!   fp          Spin-independent WIMP-proton coupling [GeV^-2].
!   fn          Spin-independent WIMP-neutron coupling [GeV^-2].
!   ap          Spin-dependent WIMP-proton coupling [unitless].
!   an          Spin-dependent WIMP-neutron coupling [unitless].
! 
SUBROUTINE DDCalc_SetWIMP_mfa(WIMP,m,fp,fn,ap,an)
  IMPLICIT NONE
  REAL*8, INTENT(IN) :: m,fp,fn,ap,an
  TYPE(WIMPStruct), INTENT(INOUT) :: WIMP
  CHARACTER(LEN=24) :: DMtype = 'SISD'
  CALL DDCalc_SetWIMP(WIMP,m=m,DMtype=DMtype,params=[fp,fn,ap,an])
END SUBROUTINE

SUBROUTINE DDCalc_GetWIMP_mfa(WIMP,m,fp,fn,ap,an)
  IMPLICIT NONE
  REAL*8, INTENT(OUT) :: m,fp,fn,ap,an
  TYPE(WIMPStruct), INTENT(IN) :: WIMP
  CHARACTER(LEN=24) :: DMtype
  REAL*8, ALLOCATABLE :: params(:)
  LOGICAL :: found = .FALSE.

  ALLOCATE(params(WIMP%Nparams))

  CALL DDCalc_GetWIMP(WIMP,m=m,DMtype=DMtype,params=params)

  IF ( DMtype .EQ. 'SIonly' ) THEN
    fp = params(1)
    fn = params(2)
    ap = 0
    an = 0
    found = .TRUE.
  END IF
  IF ( DMtype .EQ. 'SDonly' ) THEN
    fp = 0
    fn = 0
    ap = params(1)
    an = params(2)
    found = .TRUE.
  END IF
  IF ( DMtype .EQ. 'SISD' ) THEN
    fp = params(1)
    fn = params(2)
    ap = params(3)
    an = params(4)
    found = .TRUE.
  END IF

  IF ( .NOT. found ) stop 'Invalid WIMP type given to DDCalc_GetWIMP_mfa' 

END SUBROUTINE

! C++ interface wrappers
SUBROUTINE C_DDCalc_SetWIMP_mfa(WIMPIndex,m,fp,fn,ap,an) &
           BIND(C,NAME='C_DDCalc_ddcalc_setwimp_mfa')
  USE ISO_C_BINDING, only: C_DOUBLE, C_INT
  IMPLICIT NONE
  REAL(KIND=C_DOUBLE), INTENT(IN) :: m,fp,fn,ap,an
  INTEGER(KIND=C_INT), INTENT(IN) :: WIMPIndex
  IF (.NOT. ASSOCIATED(WIMPs(WIMPIndex)%p)) stop 'Invalid WIMP index given to C_DDCalc_SetWIMP_mfa' 
  CALL DDCalc_SetWIMP_mfa(WIMPs(WIMPIndex)%p,m=REAL(m,KIND=8),           &
               fp=REAL(fp,KIND=8),fn=REAL(fn,KIND=8),                  &
               ap=REAL(ap,KIND=8),an=REAL(an,KIND=8))
END SUBROUTINE

SUBROUTINE C_DDCalc_GetWIMP_mfa(WIMPIndex,m,fp,fn,ap,an) &
           BIND(C,NAME='C_DDCalc_ddcalc_getwimp_mfa')
  USE ISO_C_BINDING, only: C_DOUBLE, C_INT
  IMPLICIT NONE
  REAL(KIND=C_DOUBLE), INTENT(OUT) :: m,fp,fn,ap,an
  INTEGER(KIND=C_INT), INTENT(IN) :: WIMPIndex
  REAL*8 :: m0,fp0,fn0,ap0,an0
  IF (.NOT. ASSOCIATED(WIMPs(WIMPIndex)%p)) stop 'Invalid WIMP index given to C_DDCalc_GetWIMP_mfa' 
  CALL DDCalc_GetWIMP_mfa(WIMPs(WIMPIndex)%p,m=m0,fp=fp0,fn=fn0,ap=ap0,an=an0)
  ! Automatic type conversions here
  m  = m0
  fp = fp0
  fn = fn0
  ap = ap0
  an = an0
END SUBROUTINE

!-----------------------------------------------------------------------
! Sets/gets the WIMP mass and couplings for spin-independent and spin-dependent 
! interactions with a light mediator.
!
! Input/output arguments:
!   mdm         WIMP mass [GeV].
!   gp          Spin-independent WIMP-proton coupling [unitless].
!   gn          Spin-independent WIMP-neutron coupling [unitless].
!   mmed        Mediator mass [GeV].
! 
SUBROUTINE DDCalc_SetWIMP_longrange(WIMP,mdm,gp,gn,mmed)
  IMPLICIT NONE
  REAL*8, INTENT(IN) :: mdm,gp,gn,mmed
  TYPE(WIMPStruct), INTENT(INOUT) :: WIMP
  CHARACTER(LEN=24) :: DMtype = 'SILR'
  CALL DDCalc_SetWIMP(WIMP,m=mdm,DMtype=DMtype,params=[gp,gn,mmed])
END SUBROUTINE

SUBROUTINE DDCalc_GetWIMP_longrange(WIMP,mdm,gp,gn,mmed)
  IMPLICIT NONE
  REAL*8, INTENT(OUT) :: mdm,gp,gn,mmed
  TYPE(WIMPStruct), INTENT(IN) :: WIMP
  CHARACTER(LEN=24) :: DMtype
  REAL*8, ALLOCATABLE :: params(:)
  LOGICAL :: found = .FALSE.

  IF ( WIMP%DMtype .EQ. 'SILR' ) THEN
    gp = params(1)
    gn = params(2)
    mmed = params(3)
  ELSE
    stop 'Invalid WIMP type given to DDCalc_GetWIMP_longrange' 
  END IF

END SUBROUTINE

! C++ interface wrappers
SUBROUTINE C_DDCalc_SetWIMP_longrange(WIMPIndex,mdm,gp,gn,mmed) &
           BIND(C,NAME='C_DDCalc_ddcalc_setwimp_longrange')
  USE ISO_C_BINDING, only: C_DOUBLE, C_INT
  IMPLICIT NONE
  REAL(KIND=C_DOUBLE), INTENT(IN) :: mdm,gp,gn,mmed
  INTEGER(KIND=C_INT), INTENT(IN) :: WIMPIndex
  IF (.NOT. ASSOCIATED(WIMPs(WIMPIndex)%p)) stop 'Invalid WIMP index given to C_DDCalc_SetWIMP_longrange' 
  CALL DDCalc_SetWIMP_longrange(WIMPs(WIMPIndex)%p,mdm=REAL(mdm,KIND=8),&
               gp=REAL(gp,KIND=8),gn=REAL(gn,KIND=8),                   &
               mmed=REAL(mmed,KIND=8))
END SUBROUTINE

SUBROUTINE C_DDCalc_GetWIMP_longrange(WIMPIndex,mdm,gp,gn,mmed) &
           BIND(C,NAME='C_DDCalc_ddcalc_getwimp_longrange')
  USE ISO_C_BINDING, only: C_DOUBLE, C_INT
  IMPLICIT NONE
  REAL(KIND=C_DOUBLE), INTENT(OUT) :: mdm,gp,gn,mmed
  INTEGER(KIND=C_INT), INTENT(IN) :: WIMPIndex
  REAL*8 :: mdm0,gp0,gn0,mmed0
  IF (.NOT. ASSOCIATED(WIMPs(WIMPIndex)%p)) stop 'Invalid WIMP index given to C_DDCalc_GetWIMP_longrange' 
  CALL DDCalc_GetWIMP_longrange(WIMPs(WIMPIndex)%p,mdm=mdm0,gp=gp0,gn=gn0,mmed=mmed0)
  ! Automatic type conversions here
  mdm  = mdm0
  gp = gp0
  gn = gn0
  mmed = mmed0
END SUBROUTINE

!-----------------------------------------------------------------------
! Sets/gets the WIMP mass and couplings.  Couplings are specified via
! their effective 4-fermion vertex couplings 'G'.
! 
! For more detailed WIMP settings, see:
!   SetWIMP() [interface name: DDCalc_SetWIMP]
!   GetWIMP() [interface name: DDCalc_GetWIMP]
! 
! Input/output arguments:
!   m           WIMP mass [GeV].
!   GpSI        Spin-independent WIMP-proton coupling [GeV^-2].
!               Related by GpSI = 2 fp.
!   GnSI        Spin-independent WIMP-neutron coupling [GeV^-2].
!               Related by GnSI = 2 fn.
!   GpSD        Spin-dependent WIMP-proton coupling [GeV^-2].
!               Related by GpSD = 2\sqrt{2} G_F ap.
!   GnSD        Spin-dependent WIMP-neutron coupling [GeV^-2].
!               Related by GnSD = 2\sqrt{2} G_F an.
! 
SUBROUTINE DDCalc_SetWIMP_mG(WIMP,m,GpSI,GnSI,GpSD,GnSD)
  IMPLICIT NONE
  TYPE(WIMPStruct), INTENT(INOUT) :: WIMP
  REAL*8, INTENT(IN) :: m,GpSI,GnSI,GpSD,GnSD
  CALL DDCalc_SetWIMP_mfa(WIMP,m,GtoF(GpSI),GtoF(GnSI),GtoA(GpSD),GtoA(GnSD))
END SUBROUTINE

SUBROUTINE DDCalc_GetWIMP_mG(WIMP,m,GpSI,GnSI,GpSD,GnSD)
  IMPLICIT NONE
  TYPE(WIMPStruct), INTENT(IN) :: WIMP
  REAL*8, INTENT(OUT) :: m,GpSI,GnSI,GpSD,GnSD
  REAL*8 :: fp, fn, ap, an
  CALL DDCalc_GetWIMP_mfa(WIMP,m,fp,fn,ap,an)
  GpSI = FtoG(fp)
  GnSI = FtoG(fn)
  GpSD = AtoG(ap)
  GnSD = AtoG(an)
END SUBROUTINE

! C++ interface wrappers
SUBROUTINE C_DDCalc_SetWIMP_mG(WIMPIndex,m,GpSI,GnSI,GpSD,GnSD) &
           BIND(C,NAME='C_DDCalc_ddcalc_setwimp_mg')
  USE ISO_C_BINDING, only: C_DOUBLE, C_INT
  IMPLICIT NONE
  REAL(KIND=C_DOUBLE), INTENT(IN) :: m,GpSI,GnSI,GpSD,GnSD
  INTEGER(KIND=C_INT), INTENT(IN) :: WIMPIndex
  IF (.NOT. ASSOCIATED(WIMPs(WIMPIndex)%p)) stop 'Invalid WIMP index given to C_DDCalc_SetWIMP_mG' 
  CALL DDCalc_SetWIMP_mG(WIMPs(WIMPIndex)%p,m=REAL(m,KIND=8),            &
               GpSI=REAL(GpSI,KIND=8),GnSI=REAL(GnSI,KIND=8),          &
               GpSD=REAL(GpSD,KIND=8),GnSD=REAL(GnSD,KIND=8))
END SUBROUTINE

SUBROUTINE C_DDCalc_GetWIMP_mG(WIMPIndex,m,GpSI,GnSI,GpSD,GnSD) &
           BIND(C,NAME='C_DDCalc_ddcalc_getwimp_mg')
  USE ISO_C_BINDING, only: C_DOUBLE, C_INT
  IMPLICIT NONE
  REAL(KIND=C_DOUBLE), INTENT(OUT) :: m,GpSI,GnSI,GpSD,GnSD
  INTEGER(KIND=C_INT), INTENT(IN) :: WIMPIndex
  REAL*8 :: m0,GpSI0,GnSI0,GpSD0,GnSD0
  IF (.NOT. ASSOCIATED(WIMPs(WIMPIndex)%p)) stop 'Invalid WIMP index given to C_DDCalc_GetWIMP_mG' 
  CALL DDCalc_GetWIMP_mG(WIMPs(WIMPIndex)%p,m=m0,GpSI=GpSI0,GnSI=GnSI0,GpSD=GpSD0,GnSD=GnSD0)
  ! Automatic type conversions here
  m    = m0
  GpSI = GpSI0
  GnSI = GnSI0
  GpSD = GpSD0
  GnSD = GnSD0
END SUBROUTINE


!-----------------------------------------------------------------------
! Sets/gets the WIMP mass and couplings for standard spin-independent
! and spin-dependent interactions, specified via the spin-independent and
! and spin-dependent scattering cross sections in [pb].  
! A negative sign of a cross sections corresponds to a negative value of
! the corresponding coupling parameter.
! 
! For more detailed WIMP settings, see:
!   SetWIMP() [interface name: DDCalc_SetWIMP]
!   GetWIMP() [interface name: DDCalc_GetWIMP]
! 
! Input/output arguments:
!   m           WIMP mass [GeV].
!   sigmaSIp      Spin-independent WIMP-proton cross section [pb].
!   sigmaSIn      Spin-independent WIMP-neutron cross section [pb].
!   sigmaSDp      Spin-dependent WIMP-proton cross section [pb].
!   sigmaSDn      Spin-dependent WIMP-neutron cross section [pb].
! 
SUBROUTINE DDCalc_SetWIMP_msigma(WIMP,m,sigmaSIp,sigmaSIn,sigmaSDp,sigmaSDn)
  IMPLICIT NONE
  REAL*8, INTENT(IN) :: m,sigmaSIp,sigmaSIn,sigmaSDp,sigmaSDn
  REAL*8 :: fp,fn,ap,an
  TYPE(WIMPStruct), INTENT(INOUT) :: WIMP
  CHARACTER(LEN=24) :: DMtype = 'SISD'
 
  fp = SigmapSItoFp(m,sigmaSIp)
  fn = SigmanSItoFn(m,sigmaSIn)
  ap = SigmapSDtoAp(m,sigmaSDp)
  an = SigmanSDtoAn(m,sigmaSDn)

  CALL DDCalc_SetWIMP(WIMP,m=m,DMtype=DMtype,params=[fp,fn,ap,an])
END SUBROUTINE

SUBROUTINE DDCalc_GetWIMP_msigma(WIMP,m,sigmaSIp,sigmaSIn,sigmaSDp,sigmaSDn)
  IMPLICIT NONE
  REAL*8, INTENT(OUT) :: m,sigmaSIp,sigmaSIn,sigmaSDp,sigmaSDn
  REAL*8 :: fp,fn,ap,an
  TYPE(WIMPStruct), INTENT(IN) :: WIMP
  CHARACTER(LEN=24) :: DMtype
  REAL*8, ALLOCATABLE :: params(:)
  LOGICAL :: found = .FALSE.

  ALLOCATE(params(WIMP%Nparams))

  CALL DDCalc_GetWIMP(WIMP,m=m,DMtype=DMtype,params=params)

  IF ( DMtype .EQ. 'SISD' ) THEN
    fp = params(1)
    fn = params(2)
    ap = params(3)
    an = params(4)
    found = .TRUE.
  END IF

  IF ( .NOT. found ) stop 'Invalid WIMP type given to DDCalc_GetWIMP_msigma' 

  sigmaSIp = GpToSigmapSI(m,FtoG(fp))
  sigmaSIn = GnToSigmanSI(m,FtoG(fn))
  sigmaSDp = GpToSigmapSD(m,AtoG(ap))
  sigmaSDn = GnToSigmanSD(m,AtoG(an))

END SUBROUTINE

! C++ interface wrappers
SUBROUTINE C_DDCalc_SetWIMP_msigma(WIMPIndex,m,sigmaSIp,sigmaSIn,sigmaSDp,sigmaSDn) &
           BIND(C,NAME='C_DDCalc_ddcalc_setwimp_msigma')
  USE ISO_C_BINDING, only: C_DOUBLE, C_INT
  IMPLICIT NONE
  REAL(KIND=C_DOUBLE), INTENT(IN) :: m,sigmaSIp,sigmaSIn,sigmaSDp,sigmaSDn
  INTEGER(KIND=C_INT), INTENT(IN) :: WIMPIndex
  IF (.NOT. ASSOCIATED(WIMPs(WIMPIndex)%p)) stop 'Invalid WIMP index given to C_DDCalc_SetWIMP_msigma' 
  CALL DDCalc_SetWIMP_msigma(WIMPs(WIMPIndex)%p,m=REAL(m,KIND=8),           &
               sigmaSIp=REAL(sigmaSIp,KIND=8),sigmaSIn=REAL(sigmaSIn,KIND=8),                  &
               sigmaSDp=REAL(sigmaSDp,KIND=8),sigmaSDn=REAL(sigmaSDn,KIND=8))
END SUBROUTINE

SUBROUTINE C_DDCalc_GetWIMP_msigma(WIMPIndex,m,sigmaSIp,sigmaSIn,sigmaSDp,sigmaSDn) &
           BIND(C,NAME='C_DDCalc_ddcalc_getwimp_msigma')
  USE ISO_C_BINDING, only: C_DOUBLE, C_INT
  IMPLICIT NONE
  REAL(KIND=C_DOUBLE), INTENT(OUT) :: m,sigmaSIp,sigmaSIn,sigmaSDp,sigmaSDn
  INTEGER(KIND=C_INT), INTENT(IN) :: WIMPIndex
  REAL*8 :: m0,sigmaSIp0,sigmaSIn0,sigmaSDp0,sigmaSDn0
  IF (.NOT. ASSOCIATED(WIMPs(WIMPIndex)%p)) stop 'Invalid WIMP index given to C_DDCalc_GetWIMP_msigmas' 
  CALL DDCalc_GetWIMP_msigma(WIMPs(WIMPIndex)%p,m=m0,sigmaSIp=sigmaSIp0,&
       sigmaSIn=sigmaSIn0,sigmaSDp=sigmaSDp0,sigmaSDn=sigmaSDn0)
  ! Automatic type conversions here
  m  = m0
  sigmaSIp = sigmaSIp0
  sigmaSIn = sigmaSIn0
  sigmaSDp = sigmaSDp0
  sigmaSDn = sigmaSDn0
END SUBROUTINE





!-----------------------------------------------------------------------
! Functions for handling WIMPs in the non-relativistic effective theory
!-----------------------------------------------------------------------
!
! There are two non-equivalent ways to describe the interactions of WIMPs in the non-relativistic limit.
!
! WIMP type 'NREffectiveTheory' corresponds to the operators first defined in arXiv:1203.3542 and subsequently
! extended in arXiv:1308.6288 and arXiv:1505.03117. These operators give the most general interactions of
! DM particles with spin 0, 1/2 or 1 up to first order in velocity and second order in momentum transfer.
! 
! WIMP type 'NREFT_CPT' is also based on the operators defined in arXiv:1203.3542 for DM particles with
! spin 0 and spin 1/2 but in addition contains operators with a non-trivial momentum dependence, such as 
! 1/(m_pi^2 + q^2), which arise from chiral perturbation theory. See arXiv:1708.02678 for details.
!
! While the two sets of effective operators agree for O_1 to O_12, the definition of operators with higher index differs.
! This makes it necessary to first define the WIMP type and then pass individual operator coefficients.
! 
! DDCalc_SetWIMP_NREffectiveTheory initializes a WIMP of type 'NREffectiveTheory" and sets all coefficients to zero.
! The function requires specification of the WIMP spin and mass.
!
! DDCalc_SetWIMP_NREFT_CPT initializes a WIMP of type 'NREFT_CPT" and sets all coefficients to zero.
! The function requires specification of the WIMP spin and mass.
!
! DDCalc_SetNRCoefficient sets the value of a single operator to a given value
!   in units GeV^(-2).
! Here, OpIndex is an integer specifying the operator, e.g. 6 for O_6.
! 
! For WIMP type 'NREffectiveTheory', OpIndex = -1 stands for q^2*O_1, and 
! OpIndex = -4 stands for q^2*O_4
! 
! Furthermore, the function requires an index tau, which is interpreted differently depending on the WIMP type
! For WIMP type 'NREffectiveTheory' tau=0 for isoscalar and tau=1 for isovector
! For WIMP type 'NREFT_CPT' tau=0 for proton and tau=1 for neutron.
!
! DDCalc_GetNRCoefficient gets the values for tau=0 and tau=1 of a given operator.
! 
SUBROUTINE DDCalc_SetWIMP_NREffectiveTheory(WIMP,m,spin)
  IMPLICIT NONE
  REAL*8, INTENT(IN) :: m,spin
  TYPE(WIMPStruct), INTENT(INOUT) :: WIMP
  CHARACTER(LEN=24) :: DMtype = 'NREffectiveTheory'
  REAL*8 :: params(45)

  params(:) = 0
  params(1) = spin
  CALL DDCalc_SetWIMP(WIMP,m=m,DMtype=DMtype,params=params)
END SUBROUTINE

SUBROUTINE DDCalc_SetWIMP_NREFT_CPT(WIMP,m,spin)
  IMPLICIT NONE
  REAL*8, INTENT(IN) :: m,spin
  TYPE(WIMPStruct), INTENT(INOUT) :: WIMP
  CHARACTER(LEN=24) :: DMtype = 'NREFT_CPT'
  REAL*8 :: params(51)

  params(:) = 0
  params(51) = spin
  CALL DDCalc_SetWIMP(WIMP,m=m,DMtype=DMtype,params=params)
END SUBROUTINE

SUBROUTINE DDCalc_SetNRCoefficient(WIMP, OpIndex, tau, value)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: OpIndex, tau
  REAL*8, INTENT(IN) :: value
  TYPE(WIMPStruct), INTENT(INOUT) :: WIMP
  INTEGER :: i

  IF (( WIMP%DMtype .NE. 'NREffectiveTheory' ) .AND. ( WIMP%DMtype .NE. 'NREFT_CPT' )) THEN 
   stop 'Error in DDCalc_SetNRCoefficient: WIMP is not of type NREffectiveTheory or NREFT_CPT.'
  END IF

  CALL NRET_GetParamIndex(WIMP%DMtype, OpIndex, tau, i)
  WIMP%params(i) = value

  IF ( WIMP%DMtype .EQ. 'NREffectiveTheory' ) THEN 
    CALL NRET_UpdateNRCoefficients(WIMP%params)
  END IF

END SUBROUTINE

SUBROUTINE DDCalc_GetNRCoefficient(WIMP, OpIndex, value_0, value_1)
  IMPLICIT NONE
  REAL*8, INTENT(OUT) :: value_0, value_1
  TYPE(WIMPStruct), INTENT(IN) :: WIMP
  INTEGER, INTENT(IN) :: OpIndex
  INTEGER :: index_0, index_1

  IF (( WIMP%DMtype .NE. 'NREffectiveTheory' ) .AND. ( WIMP%DMtype .NE. 'NREFT_CPT' )) THEN 
   stop 'Error in DDCalc_GetNRCoefficient: WIMP is not of type NREffectiveTheory or NREFT_CPT.'
  END IF

  CALL NRET_GetParamIndex(WIMP%DMtype, OpIndex, 0, index_0)
  CALL NRET_GetParamIndex(WIMP%DMtype, OpIndex, 1, index_1)
  value_0 = WIMP%params(index_0)
  value_1 = WIMP%params(index_1)

END SUBROUTINE

! C++ interface wrappers
SUBROUTINE C_DDCalc_SetWIMP_NREffectiveTheory(WIMPIndex,m,spin) &
           BIND(C,NAME='C_DDCalc_ddcalc_setwimp_nreffectivetheory')
  USE ISO_C_BINDING, only: C_DOUBLE, C_INT
  IMPLICIT NONE
  REAL(KIND=C_DOUBLE), INTENT(IN) :: m,spin
  INTEGER(KIND=C_INT), INTENT(IN) :: WIMPIndex
  IF (.NOT. ASSOCIATED(WIMPs(WIMPIndex)%p)) &
    stop 'Invalid WIMP index given to C_DDCalc_SetWIMP_NREffectiveTheory' 
  CALL DDCalc_SetWIMP_NREffectiveTheory(WIMPs(WIMPIndex)%p,&
               m=REAL(m,KIND=8),spin=REAL(spin,KIND=8))
END SUBROUTINE

SUBROUTINE C_DDCalc_SetWIMP_NREFT_CPT(WIMPIndex,m,spin) &
           BIND(C,NAME='C_DDCalc_ddcalc_setwimp_nreft_cpt')
  USE ISO_C_BINDING, only: C_DOUBLE, C_INT
  IMPLICIT NONE
  REAL(KIND=C_DOUBLE), INTENT(IN) :: m,spin
  INTEGER(KIND=C_INT), INTENT(IN) :: WIMPIndex
  IF (.NOT. ASSOCIATED(WIMPs(WIMPIndex)%p)) &
    stop 'Invalid WIMP index given to C_DDCalc_SetWIMP_NREFT_CPT' 
  CALL DDCalc_SetWIMP_NREFT_CPT(WIMPs(WIMPIndex)%p,&
               m=REAL(m,KIND=8),spin=REAL(spin,KIND=8))
END SUBROUTINE

SUBROUTINE C_DDCalc_SetNRCoefficient(WIMPIndex, OpIndex, tau, value) &
           BIND(C,NAME='C_DDCalc_ddcalc_setnrcoefficient')
  USE ISO_C_BINDING, only: C_DOUBLE, C_INT
  IMPLICIT NONE
  INTEGER(KIND=C_INT), INTENT(IN) :: WIMPIndex, OpIndex, tau
  REAL(KIND=C_DOUBLE), INTENT(IN) :: value
  IF (.NOT. ASSOCIATED(WIMPs(WIMPIndex)%p)) &
    stop 'Invalid WIMP index given to C_DDCalc_SetNRCoefficient' 
  CALL DDCalc_SetNRCoefficient(WIMPs(WIMPIndex)%p,OpIndex,&
               tau,value=REAL(value,KIND=8))
END SUBROUTINE


SUBROUTINE C_DDCalc_GetNRCoefficient(WIMPIndex,OpIndex, value_0, value_1) &
           BIND(C,NAME='C_DDCalc_ddcalc_getnrcoefficient')
  USE ISO_C_BINDING, only: C_DOUBLE, C_INT
  IMPLICIT NONE
  REAL(KIND=C_DOUBLE), INTENT(OUT) :: value_0, value_1
  INTEGER(KIND=C_INT), INTENT(IN) :: WIMPIndex, OpIndex
  REAL*8 :: value_0_internal, value_1_internal
  IF (.NOT. ASSOCIATED(WIMPs(WIMPIndex)%p)) stop 'Invalid WIMP index given to C_DDCalc_GetNRCoefficient' 
  CALL DDCalc_GetNRCoefficient(WIMPs(WIMPIndex)%p,OpIndex,&
    value_0=value_0_internal,value_1=value_1_internal)
  ! Automatic type conversions here
  value_0 = value_0_internal
  value_1 = value_1_internal
END SUBROUTINE




! ----------------------------------------------------------------------
! Sets the minimum recoil energy to be included in the calculations.
! Note the efficiency curves already account for detector and analysis
! thresholds regardless of this setting, so setting this to 0 keV (the
! default behavior when initialization is performed) does not imply
! that very low energy recoils actually contribute to the signal.
! 
! Required input arguments:
!     Emin        The minimum recoil energy to consider [keV]
! 
SUBROUTINE DDCalc_SetDetectorEmin(Detector,Emin)
  IMPLICIT NONE
  REAL*8, INTENT(IN) :: Emin
  TYPE(DetectorStruct), INTENT(INOUT) :: Detector
  CALL DDCalc_SetDetector(Detector,Emin=Emin)
END SUBROUTINE

! C++ interface wrapper
SUBROUTINE C_DDCalc_SetDetectorEmin(DetectorIndex,Emin) &
           BIND(C,NAME='C_DDCalc_ddcalc_setdetectoremin')
  USE ISO_C_BINDING, only: C_DOUBLE, C_INT
  IMPLICIT NONE
  REAL(KIND=C_DOUBLE), INTENT(IN) :: Emin
  INTEGER(KIND=C_INT), INTENT(IN) :: DetectorIndex
  IF (.NOT. ASSOCIATED(Detectors(DetectorIndex)%p)) stop 'Invalid Detector index given to C_DDCalc_SetDetectorEmin' 
  CALL DDCalc_SetDetectorEmin(Detectors(DetectorIndex)%p,REAL(Emin,KIND=8))
END SUBROUTINE

!-----------------------------------------------------------------------
! Initialise default detector
! 
FUNCTION DDCalc_InitDetector() RESULT(Detector)

  IMPLICIT NONE
  TYPE(DetectorStruct) :: Detector
  Detector = DDCalc_InitDefaultDetector()

END FUNCTION

!-----------------------------------------------------------------------
! C/C++ wrapper for DDCalc_InitDetector
!
INTEGER(KIND=C_INT) FUNCTION C_DDCalc_InitDetector() &
 BIND(C,NAME='C_DDExperiments_ddcalc_initdetector') 
  USE ISO_C_BINDING, only: C_INT
  IMPLICIT NONE
  N_Detectors = N_Detectors + 1
  IF (N_Detectors .GT. Max_Detectors) stop 'DDCalc: Max_Detectors exceeded.&
   Please run FreeDetectors or modify Max_Detectors in DDTypes.f90.'
  ALLOCATE(Detectors(N_Detectors)%p)
  Detectors(N_Detectors)%p = DDCalc_InitDetector()
  C_DDCalc_InitDetector = N_Detectors
END FUNCTION

END MODULE


 
