!#######################################################################
MODULE DDCalc
!#######################################################################

!#######################################################################
! DIRECT DETECTION RATES AND LIKELIHOODS
! Routines to calculate the dark matter direct detection event rates
! and corresponding likelihoods/exclusion levels.
! 
! To compile the command-line interface 'DDCalc_run':
!     make DDCalc_run
! 
! To see usage:
!     ./DDCalc_run --help
! 
!   Created by Chris Savage
!   University of Utah   (2013 - 2014)
!   Nordita              (2014 - 2015)
! 
!   With contributions from:
!     Andre Scaffidi            University of Adelaide (2014)
!     Lauren Hsu                FermiLab (2015)
!     Pat Scott                 Imperial College London (2016)
! 
! 
! 
!-----------========<<<<<<<< FILE LAYOUT >>>>>>>>========---------------
! 
!  * SIMPLE INTERFACE ROUTINES
!    Basic routines that can be used to initialise objects, set
!    parameters, do caluclations and read the results.  Meant to
!    provide a simple, basic interface for external access to the
!    package.
!  * MAIN ROUTINES
!    Main routines used for each program mode.
!  * INITIALIZATION
!    General initialization routines.
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
!     SUBROUTINE DDCalc_SetWIMP_msigma(WIMP,m,sigmapSI,sigmanSI,sigmapSD,sigmanSD)
!   where m is the WIMP mass [GeV], f is spin-independent (SI) WIMP-
!   nucleon coupling [GeV^-2], a is the spin-dependent (SD) WIMP-nucleon
!   coupling [unitless], G are the effective 4-fermion vertex couplings,
!   related by:
!     GpSI = 2 fp        GpSD = 2\sqrt{2} G_F ap
!     GnSI = 2 fn        GnSD = 2\sqrt{2} G_F an
!   and sigma is the WIMP-nucleon scattering cross-section [pb].
!   Negative cross-sections indicate the corresponding coupling should
!   be negative.  In all cases, 'p' refers to proton and 'n' to neutron.
!
! WIMP parameter retrieval:
!     SUBROUTINE DDCalc_GetWIMP_mfa(WIMP,m,fp,fn,ap,an)
!     SUBROUTINE DDCalc_GetWIMP_mG(WIMP,m,GpSI,GnSI,GpSD,GnSD)
!     SUBROUTINE DDCalc_GetWIMP_msigma(WIMP,m,sigmapSI,sigmanSI,sigmapSD,sigmanSD)
!   Same as SetWIMP above, but retrieves current WIMP parameters.
!   The only difference: cross-sections are always positive, regardless
!   of the sign of the corresponding coupling.
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
!     TYPE(DetectorStruct) FUNCTION DDCalc_InitDetector(intervals)
! Initialise an object carrying the information about the default
! experimental analysis to consider.  Here intervals is a flag
! indicating if calculations should be performed for analysis
! sub-intervals (i.e. intervals between observed events).  This is
! only necessary for maximum gap calculations and can be set to .FALSE.
! for likelihood analyses. Returns a detector structure containing the
! analysis details.  Non-default analyses can be obtained with
! specific functions:
!     TYPE(DetectorStruct) FUNCTION [Analysis]_Init(intervals)
! See DDExperiments.f90 and analyses/[Analysis].f90 for more details.
! Note that these functions are only available directly by USE-ing the
! DDExperiments module in the calling program, and will not be made
! available simply by USE-ing the DDCalc module.  (They can however be
! accessed indirectly by the various Main programs later in this file,
! but those are part of the advanced interface.)
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
!    Or the separate spin-independent and spin-dependent contributions:
!     REAL*8 FUNCTION DDCalc_SignalSI(Detector)
!     REAL*8 FUNCTION DDCalc_SignalSD(Detector)
!
! 4. Log-likelihood:
!     REAL*8 FUNCTION DDCalc_LogLikelihood(Detector)
!    Uses a Poisson distribution in the number of observed events N:
!      P(N|s+b)
!    where s is the average expected signal and b is the average expected
!    background.
!
! 5. Logarithm of the p-value:
!     REAL*8 FUNCTION DDCalc_LogPValue(Detector)
!    Uses the maximum gap method if <Detector>_Init was called with
!    argument intervals set .TRUE. and the anlysis contains the necessary
!    interval information to allow such a method; otherwise uses a Poisson
!    distribution in the number of observed events N:
!      P(N|s),
!    where s is the average expected signal (background contributions are
!    ignored).
!
! 6. Factor by which the WIMP cross-sections must be multiplied to
!    achieve a given p-value:
!     REAL*8 FUNCTION DDCalc_ScaleToPValue(lnp)
!    Calculates the factor x by which the cross-sections must be scaled
!    (sigma -> x*sigma) to achieve the desired p-value (given as log(p)).
!    See DDCalc_LogPValue() above for a description of the statistics.
!
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
USE DDCommandLine
USE DDNumerical
USE DDExHelp
USE DDWIMP
USE DDDetectors
USE DDRates
USE DDStats
USE DDHalo
USE DDOutput
USE DDTabulation

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
PUBLIC :: DDCalc_SetWIMP_mfa
PUBLIC :: DDCalc_SetWIMP_mG
PUBLIC :: DDCalc_SetWIMP_msigma
PUBLIC :: DDCalc_GetWIMP_mfa
PUBLIC :: DDCalc_GetWIMP_mG
PUBLIC :: DDCalc_GetWIMP_msigma
PUBLIC :: DDCalc_SetDetectorEmin
! C/C++ Wrappers
PUBLIC :: C_DDCalc_SetSHM
PUBLIC :: C_DDCalc_SetWIMP_mfa
PUBLIC :: C_DDCalc_GetWIMP_mfa
PUBLIC :: C_DDCalc_SetWIMP_mG
PUBLIC :: C_DDCalc_GetWIMP_mG
PUBLIC :: C_DDCalc_SetWIMP_msigma
PUBLIC :: C_DDCalc_GetWIMP_msigma
PUBLIC :: C_DDCalc_SetDetectorEmin

! Rate calculation
PUBLIC :: DDCalc_CalcRates
! C/C++ Wrapper
PUBLIC :: C_DDCalc_CalcRates

! Result inspection
PUBLIC :: DDCalc_Events
PUBLIC :: DDCalc_Background
PUBLIC :: DDCalc_Signal
PUBLIC :: DDCalc_SignalSI
PUBLIC :: DDCalc_SignalSD
PUBLIC :: DDCalc_LogLikelihood
PUBLIC :: DDCalc_LogPValue
PUBLIC :: DDCalc_ScaleToPValue
! C/C++ Wrappers
PUBLIC :: C_DDCalc_Events
PUBLIC :: C_DDCalc_Background
PUBLIC :: C_DDCalc_Signal
PUBLIC :: C_DDCalc_SignalSI
PUBLIC :: C_DDCalc_SignalSD
PUBLIC :: C_DDCalc_LogLikelihood
PUBLIC :: C_DDCalc_LogPValue
PUBLIC :: C_DDCalc_ScaleToPValue

! Memory cleanup for C/C++ interface
PUBLIC :: C_DDCalc_FreeWIMPs
PUBLIC :: C_DDCalc_FreeHalos
PUBLIC :: C_DDCalc_FreeDetectors
PUBLIC :: C_DDCalc_FreeAll

!############## Advanced Interface ####################
!               (no C/C++ wrappers)

! More advanced setter and getter interfaces
PUBLIC :: DDCalc_GetWIMP
PUBLIC :: DDCalc_SetWIMP
PUBLIC :: DDCalc_GetHalo
PUBLIC :: DDCalc_SetHalo
PUBLIC :: DDCalc_GetDetector
PUBLIC :: DDCalc_SetDetector

! Main routines (each corresponds to a different run mode)
PUBLIC :: DDCalc_Main
PUBLIC :: DDCalc_MainEventsAndLikelihoods                
PUBLIC :: DDCalc_MainEventsAndLikelihoodsInteractive                  
PUBLIC :: DDCalc_MainLogLikelihood
PUBLIC :: DDCalc_MainLogLikelihoodInteractive
PUBLIC :: DDCalc_MainLogPValue
PUBLIC :: DDCalc_MainLogPValueInteractive
PUBLIC :: DDCalc_MainSpectrum
PUBLIC :: DDCalc_MainEventsByMass
PUBLIC :: DDCalc_MainConstraintsSI
PUBLIC :: DDCalc_MainConstraintsSD
PUBLIC :: DDCalc_MainLimitsSI
PUBLIC :: DDCalc_MainLimitsSD

! Initialization routines
PRIVATE :: InitializeCommandLine
PRIVATE :: InitVerbosity
PRIVATE :: ReadArguments
PRIVATE :: ReadLogPValue

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
! Sets/gets the WIMP mass and couplings.  Couplings are specified via
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
!               Related by GpSI = 2 fp.
!   fn          Spin-independent WIMP-neutron coupling [GeV^-2].
!               Related by GnSI = 2 fn.
!   ap          Spin-dependent WIMP-proton coupling [unitless].
!               Related by GpSD = 2\sqrt{2} G_F ap.
!   an          Spin-dependent WIMP-neutron coupling [unitless].
!               Related by GnSD = 2\sqrt{2} G_F an.
! 
SUBROUTINE DDCalc_SetWIMP_mfa(WIMP,m,fp,fn,ap,an)
  IMPLICIT NONE
  REAL*8, INTENT(IN) :: m,fp,fn,ap,an
  TYPE(WIMPStruct), INTENT(INOUT) :: WIMP
  CALL DDCalc_SetWIMP(WIMP,m=m,fp=fp,fn=fn,ap=ap,an=an)
END SUBROUTINE

SUBROUTINE DDCalc_GetWIMP_mfa(WIMP,m,fp,fn,ap,an)
  IMPLICIT NONE
  REAL*8, INTENT(OUT) :: m,fp,fn,ap,an
  TYPE(WIMPStruct), INTENT(IN) :: WIMP
  CALL DDCalc_GetWIMP(WIMP,m=m,fp=fp,fn=fn,ap=ap,an=an)
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
  CALL DDCalc_SetWIMP(WIMP,m=m,GpSI=GpSI,GnSI=GnSI,GpSD=GpSD,GnSD=GnSD)
END SUBROUTINE

SUBROUTINE DDCalc_GetWIMP_mG(WIMP,m,GpSI,GnSI,GpSD,GnSD)
  IMPLICIT NONE
  TYPE(WIMPStruct), INTENT(IN) :: WIMP
  REAL*8, INTENT(OUT) :: m,GpSI,GnSI,GpSD,GnSD
  CALL DDCalc_GetWIMP(WIMP,m=m,GpSI=GpSI,GnSI=GnSI,GpSD=GpSD,GnSD=GnSD)
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
! Sets/gets the WIMP mass and couplings.  Couplings are specified via
! their corresponding WIMP-nucleon cross-sections.
! 
! For more detailed WIMP settings, see:
!   SetWIMP() [interface name: DDCalc_SetWIMP].
!   GetWIMP() [interface name: DDCalc_GetWIMP]
! 
! Input/output arguments.  As input, give negative value to cross-
! sections to set the corresponding coupling negative:
!   m           WIMP mass [GeV].
!   sigmapSI    Spin-independent WIMP-proton cross-section [pb].
!   sigmanSI    Spin-independent WIMP-neutron cross-section [pb].
!   sigmapSD    Spin-dependent WIMP-proton cross-section [pb].
!   sigmanSD    Spin-dependent WIMP-neutron cross-section [pb].
! 
SUBROUTINE DDCalc_SetWIMP_msigma(WIMP,m,sigmapSI,sigmanSI,sigmapSD,sigmanSD)
  IMPLICIT NONE
  TYPE(WIMPStruct), INTENT(INOUT) :: WIMP
  REAL*8, INTENT(IN) :: m,sigmapSI,sigmanSI,sigmapSD,sigmanSD
  CALL DDCalc_SetWIMP(WIMP,m=m,sigmapSI=sigmapSI,sigmanSI=sigmanSI, &
               sigmapSD=sigmapSD,sigmanSD=sigmanSD)
END SUBROUTINE

SUBROUTINE DDCalc_GetWIMP_msigma(WIMP,m,sigmapSI,sigmanSI,sigmapSD,sigmanSD)
  IMPLICIT NONE
  TYPE(WIMPStruct), INTENT(IN) :: WIMP
  REAL*8, INTENT(OUT) :: m,sigmapSI,sigmanSI,sigmapSD,sigmanSD
  CALL DDCalc_GetWIMP(WIMP,m=m,sigmapSI=sigmapSI,sigmanSI=sigmanSI, &
               sigmapSD=sigmapSD,sigmanSD=sigmanSD)
END SUBROUTINE


! C++ interface wrappers
SUBROUTINE C_DDCalc_SetWIMP_msigma(WIMPIndex,m,sigmapSI,sigmanSI,sigmapSD,sigmanSD) &
           BIND(C,NAME='C_DDCalc_ddcalc_setwimp_msigma')
  USE ISO_C_BINDING, only: C_DOUBLE, C_INT
  IMPLICIT NONE
  REAL(KIND=C_DOUBLE), INTENT(IN) :: m,sigmapSI,sigmanSI,sigmapSD,sigmanSD
  INTEGER(KIND=C_INT), INTENT(IN) :: WIMPIndex
  IF (.NOT. ASSOCIATED(WIMPs(WIMPIndex)%p)) stop 'Invalid WIMP index given to C_DDCalc_SetWIMP_msigma' 
  CALL DDCalc_SetWIMP_msigma(WIMPs(WIMPIndex)%p,m=REAL(m,KIND=8), &
               sigmapSI=REAL(sigmapSI,KIND=8),sigmanSI=REAL(sigmanSI,KIND=8),&
               sigmapSD=REAL(sigmapSD,KIND=8),sigmanSD=REAL(sigmanSD,KIND=8))
END SUBROUTINE

SUBROUTINE C_DDCalc_GetWIMP_msigma(WIMPIndex,m,sigmapSI,sigmanSI,sigmapSD,sigmanSD) &
           BIND(C,NAME='C_DDCalc_ddcalc_getwimp_msigma')
  USE ISO_C_BINDING, only: C_DOUBLE, C_INT
  IMPLICIT NONE
  REAL(KIND=C_DOUBLE), INTENT(OUT) :: m,sigmapSI,sigmanSI,sigmapSD,sigmanSD
  INTEGER(KIND=C_INT), INTENT(IN) :: WIMPIndex
  REAL*8 :: m0,sigmapSI0,sigmanSI0,sigmapSD0,sigmanSD0
  IF (.NOT. ASSOCIATED(WIMPs(WIMPIndex)%p)) stop 'Invalid WIMP index given to C_DDCalc_GetWIMP_msigma' 
  CALL DDCalc_GetWIMP_msigma(WIMPs(WIMPIndex)%p,m=m0,sigmapSI=sigmapSI0, &
                             sigmanSI=sigmanSI0,sigmapSD=sigmapSD0,    &
                             sigmanSD=sigmanSD0)
  ! Automatic type conversions here
  m        = m0
  sigmapSI = sigmapSI0
  sigmanSI = sigmanSI0
  sigmapSD = sigmapSD0
  sigmanSD = sigmanSD0
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


!=======================================================================
! MAIN ROUTINES
!=======================================================================

! ----------------------------------------------------------------------
! Generic main routine.
! 
SUBROUTINE DDCalc_Main()
  IMPLICIT NONE
  LOGICAL :: interactive
  TYPE (ArgumentStruct) :: Arguments
  
  ! Show usage and exit
  IF (ShowUsageRequested()) THEN
    CALL ShowUsage()
    STOP
  END IF
  
  ! Load arguments
  Arguments = ReadArguments()
  
  ! Interactive mode?
  ! True if --interactive given or if the WIMP mass is not specified
  ! through the command-line (either with --m=<val> option or non-option
  ! argument).
  interactive = GetLongArg('interactive') .OR.                          &
                ((Arguments%Nparameters .EQ. 0) .AND. .NOT. GetLongArg('m'))
  
  ! Determine program mode =============
  
  ! ------------------------------------
  ! Calculates the log-likelihood (w/ background)
  IF (GetLongArg('log-likelihood')) THEN
    IF (interactive) THEN
      CALL DDCalc_MainLogLikelihoodInteractive()
    ELSE
      CALL DDCalc_MainLogLikelihood()
    END IF
  ! ------------------------------------
  ! Calculates the log of the p-value (no background subtraction)
  ELSE IF (GetLongArg('log-pvalue')) THEN
    IF (interactive) THEN
      CALL DDCalc_MainLogPValueInteractive()
    ELSE
      CALL DDCalc_MainLogPValue()
    END IF
  ! ------------------------------------
  ! Prints the raw recoil spectrum dR/dE
  ELSE IF (GetLongArg('spectrum')) THEN
    CALL DDCalc_MainSpectrum()
  ! ------------------------------------
  ! Calculates expected events, tabulated by WIMP mass
  ELSE IF (GetLongArg('events-by-mass')) THEN
    CALL DDCalc_MainEventsByMass()
  ! ------------------------------------
  ! Calculates spin-independent likelihood contraints,
  ! tabulated by WIMP mass
  ELSE IF (GetLongArg('constraints-SI')) THEN
    CALL DDCalc_MainConstraintsSI()
  ! ------------------------------------
  ! Calculates spin-dependent likelihood contraints,
  ! tabulated by WIMP mass
  ELSE IF (GetLongArg('constraints-SD')) THEN
    CALL DDCalc_MainConstraintsSD()
  ! ------------------------------------
  ! Calculates spin-independent no-background-subtraction limits,
  ! tabulated by WIMP mass
  ELSE IF (GetLongArg('limits-SI')) THEN
    CALL DDCalc_MainLimitsSI()
  ! ------------------------------------
  ! Calculates spin-dependent no-background-subtraction limits,
  ! tabulated by WIMP mass
  ELSE IF (GetLongArg('limits-SD')) THEN
    CALL DDCalc_MainLimitsSD()
  ! ------------------------------------
  ! Default case:
  ! Calculate both expected events and likelihoods
  ELSE
    IF (interactive) THEN
      CALL DDCalc_MainEventsAndLikelihoodsInteractive()
    ELSE IF (Arguments%Nparameters .GE. 1) THEN
      CALL DDCalc_MainEventsAndLikelihoods()
    ELSE
      CALL DDCalc_MainEventsAndLikelihoodsInteractive()
    END IF
  END IF
  
END SUBROUTINE


! ----------------------------------------------------------------------
! Main routine to calculate log-likelihood.
! 
SUBROUTINE DDCalc_MainLogLikelihood()

  IMPLICIT NONE
  TYPE(DetectorStruct) :: Detector
  TYPE(WIMPStruct) :: WIMP
  TYPE(HaloStruct) :: Halo
  REAL*8 :: lnLike
  
  ! Show usage and exit
  IF (ShowUsageRequested()) THEN
    CALL ShowUsage()
    STOP
  END IF
  
  ! Initializes detector, halo and WIMP structures.
  ! Will calculate only total rates (not sub-intervals).
  Detector = InitializeCommandLine(WIMP, Halo, intervals=.FALSE.)
  
  ! Print header
  IF (VerbosityLevel .GE. 2) THEN
    CALL WriteCommandHeader()
    CALL WriteWIMPHeader(WIMP)
    CALL WriteHaloHeader(Halo)
    CALL WriteDetectorHeader(Detector)
    CALL WriteLogLikelihoodHeader()
  END IF
  
  ! Do rate calculations
  CALL DDCalc_CalcRates(Detector,WIMP,Halo)
  
  ! Get log-likelihood
  lnLike = DDCalc_LogLikelihood(Detector)
  
  ! Print results
  WRITE(*,*)  lnLike
  
END SUBROUTINE


! ----------------------------------------------------------------------
! Main routine to calculate log-likelihood (interactive mode).
! Reads one line of input containing WIMP parameters and writes out the
! corresponding log-likelihood, terminating when the input stream ends
! (EOF) or a blank line is given.
! 
SUBROUTINE DDCalc_MainLogLikelihoodInteractive()

  IMPLICIT NONE
  TYPE(DetectorStruct) :: Detector
  TYPE(WIMPStruct) :: WIMP
  TYPE(HaloStruct) :: Halo
  REAL*8 :: lnLike
  
  ! Show usage and exit
  IF (ShowUsageRequested()) THEN
    CALL ShowUsage()
    STOP
  END IF
  
  ! Initializes detector, halo and WIMP structures.
  ! Will calculate only total rates (not sub-intervals).
  Detector = InitializeCommandLine(WIMP,Halo,intervals=.FALSE.)
  
  ! Print header
  IF (VerbosityLevel .GE. 2) THEN
    CALL WriteCommandHeader()
    CALL WriteHaloHeader(Halo)
    CALL WriteDetectorHeader(Detector)
    CALL WriteLogLikelihoodHeader()
  END IF
  
  ! Print instructions
  IF (VerbosityLevel .GE. 1) THEN
    CALL WriteInteractiveHeader()
  END IF
  
  ! Cycle over input.
  ! ParseWIMPInput() reads line containing WIMP parameters from
  ! standard input and parses it, returning false if a blank line
  ! or EOF is found.  See ParseWIMPInput() for description of input
  ! line.
  DO WHILE (ParseWIMPInput(WIMP))
    ! Do rate calculations
    CALL DDCalc_CalcRates(Detector,WIMP,Halo)
    ! Get log-likelihood
    lnLike = DDCalc_LogLikelihood(Detector)
    ! Print results to standard output
    WRITE(*,*)  lnLike
  END DO
  
END SUBROUTINE


! ----------------------------------------------------------------------
! Main routine to calculate log of p-value using the maximum gap method.
! See S. Yellin, PRD 66, 032005 (2002) [physics/0203002].  Uses Poisson
! distribution if efficiencies for sub-intervals are not available (but
! no background subtraction!).
! 
SUBROUTINE DDCalc_MainLogPValue()

  IMPLICIT NONE
  TYPE(DetectorStruct) :: Detector
  TYPE(WIMPStruct) :: WIMP
  TYPE(HaloStruct) :: Halo
  REAL*8 :: lnP
  
  ! Show usage and exit
  IF (ShowUsageRequested()) THEN
    CALL ShowUsage()
    STOP
  END IF
  
  ! Initializes detector, halo and WIMP structures.
  ! Will calculate rates for sub-intervals as these are used by
  ! maximum gap method.
  Detector = InitializeCommandLine(WIMP,Halo,intervals=.TRUE.)
  
  ! Print header
  IF (VerbosityLevel .GE. 2) THEN
    CALL WriteCommandHeader()
    CALL WriteWIMPHeader(WIMP)
    CALL WriteHaloHeader(Halo)
    CALL WriteDetectorHeader(Detector)
    CALL WriteLogPValueHeader(Detector)
  END IF
  
  ! Do rate calculations
  CALL DDCalc_CalcRates(Detector,WIMP,Halo)
  
  ! Get log of p-value
  lnP = DDCalc_LogPValue(Detector)
  
  ! Print results
  WRITE(*,*)  lnP
  
END SUBROUTINE


! ----------------------------------------------------------------------
! Main routine to calculate log of p-value using the maximum gap method
! (interactive mode).  See S. Yellin, PRD 66, 032005 (2002)
! [physics/0203002].  Uses Poisson distribution if efficiencies for
! sub-intervals are not available (but no background subtraction!).
! Reads one line of input containing WIMP parameters and writes out the
! corresponding log of the p-value, terminating when the input stream
! ends (EOF) or a blank line is given.
! 
SUBROUTINE DDCalc_MainLogPValueInteractive()

  IMPLICIT NONE
  TYPE(DetectorStruct) :: Detector
  TYPE(WIMPStruct) :: WIMP
  TYPE(HaloStruct) :: Halo
  REAL*8 :: lnP
  
  ! Show usage and exit
  IF (ShowUsageRequested()) THEN
    CALL ShowUsage()
    STOP
  END IF
  
  ! Initializes detector, halo and WIMP structures.
  ! Will calculate rates for sub-intervals as these are used by
  ! maximum gap method.
  Detector = InitializeCommandLine(WIMP,Halo,intervals=.TRUE.)
  
  ! Print header
  IF (VerbosityLevel .GE. 2) THEN
    CALL WriteCommandHeader()
    CALL WriteHaloHeader(Halo)
    CALL WriteDetectorHeader(Detector)
    CALL WriteLogPValueHeader(Detector)
  END IF
  
  ! Print instructions
  IF (VerbosityLevel .GE. 1) THEN
    CALL WriteInteractiveHeader()
  END IF
  
  ! Cycle over input.
  ! ParseWIMPInput() reads line containing WIMP parameters from
  ! standard input and parses it, returning false if a blank line
  ! or EOF is found.  See ParseWIMPInput() for description of input
  ! line.
  DO WHILE (ParseWIMPInput(WIMP))
    ! Do rate calculations
    CALL DDCalc_CalcRates(Detector,WIMP,Halo)
    ! Get log of p-value
    lnP = DDCalc_LogPValue(Detector)
    ! Print results to standard output
    WRITE(*,*)  lnP
  END DO
  
END SUBROUTINE


! ----------------------------------------------------------------------
! Main routine to print expected events and likelihoods.
! 
SUBROUTINE DDCalc_MainEventsAndLikelihoods()

  IMPLICIT NONE
  TYPE(DetectorStruct) :: Detector
  TYPE(WIMPStruct) :: WIMP
  TYPE(HaloStruct) :: Halo

  ! Show usage and exit
  IF (ShowUsageRequested()) THEN
    CALL ShowUsage()
    STOP
  END IF
  
  ! Initializes detector, halo and WIMP structures.
  Detector = InitializeCommandLine(WIMP,Halo)
  
  ! Print header
  IF (VerbosityLevel .GE. 2) THEN
    CALL WriteCommandHeader()
    IF (VerbosityLevel .EQ. 2) CALL WriteWIMPHeader(WIMP)
    CALL WriteHaloHeader(Halo)
    CALL WriteDetectorHeader(Detector)
  END IF
  IF (VerbosityLevel .GE. 1) THEN
    CALL WriteEventsAndLikelihoodsHeader(Detector)
    CALL WriteEventsAndLikelihoodsColumnHeader()
  END IF
  
  ! Do rate calculations
  CALL DDCalc_CalcRates(Detector,WIMP,Halo)
  
  ! Print results
  CALL WriteEventsAndLikelihoodsData(Detector,WIMP)
  
END SUBROUTINE


! ----------------------------------------------------------------------
! Main routine to print expected events and likelihoods.
! 
SUBROUTINE DDCalc_MainEventsAndLikelihoodsInteractive()

  IMPLICIT NONE
  TYPE(DetectorStruct) :: Detector
  TYPE(WIMPStruct) :: WIMP
  TYPE(HaloStruct) :: Halo

  ! Show usage and exit
  IF (ShowUsageRequested()) THEN
    CALL ShowUsage()
    STOP
  END IF
  
  ! Initializes detector, halo and WIMP structures.
  Detector = InitializeCommandLine(WIMP,Halo)
  
  ! Print header
  IF (VerbosityLevel .GE. 2) THEN
    CALL WriteCommandHeader()
    CALL WriteHaloHeader(Halo)
    CALL WriteDetectorHeader(Detector)
    CALL WriteEventsAndLikelihoodsHeader(Detector)
  END IF
  
  ! Print instructions
  IF (VerbosityLevel .GE. 1) THEN
    CALL WriteInteractiveHeader(1)
  END IF
  IF (VerbosityLevel .GE. 1) THEN
    CALL WriteEventsAndLikelihoodsColumnHeader()
  END IF
  
  ! Cycle over input.
  ! ParseWIMPInput() reads line containing WIMP parameters from
  ! standard input and parses it, returning false if a blank line
  ! or EOF is found.  See ParseWIMPInput() for description of input
  ! line.
  DO WHILE (ParseWIMPInput(WIMP))
    ! Do rate calculations
    CALL DDCalc_CalcRates(Detector,WIMP,Halo)
    ! Print results to standard output
    CALL WriteEventsAndLikelihoodsData(Detector,WIMP)
  END DO
  
END SUBROUTINE


! ----------------------------------------------------------------------
! Main routine to calculate raw differential rates as a function of
! energy, i.e. the raw recoil energy spectrum (not including
! efficiencies or energy resolution).
! 
SUBROUTINE DDCalc_MainSpectrum()

  IMPLICIT NONE
  TYPE(DetectorStruct) :: Detector
  TYPE(WIMPStruct) :: WIMP
  TYPE(HaloStruct) :: Halo
  LOGICAL :: use_log
  INTEGER :: NE,K
  REAL*8 :: Emin,Emax
  REAL*8, ALLOCATABLE :: E(:),Eeff(:),eff(:,:)
  TYPE(TabulationStruct) :: TS
  
  ! Show usage and exit
  IF (ShowUsageRequested()) THEN
    CALL ShowUsage()
    STOP
  END IF
  
  ! Initializes detector, halo and WIMP structures.
  Detector = InitializeCommandLine(WIMP,Halo)
  
  ! Set tabulation energies; set efficiencies to 1.
  ! Note same results can be achieved by simply giving NE=-1
  ! and NEeff=-1 arguments to SetDetector (apart from any
  ! --E-tabulation command line specification).
  Emin = 0.1d0
  Emax = 1000d0
  NE   = -50
  use_log = .TRUE.
  CALL GetTabulationArgs('E-tabulation',Emin,Emax,NE,use_log)
  CALL InitTabulation(TS,Emin,Emax,NE,.TRUE.)
  NE = TS%N+2
  ALLOCATE(E(1:NE),Eeff(1:2),eff(1:2,0:0))
  E(1) = 0d0
  DO K = 2,NE
    E(K) = TabulationValue(TS,K-2)
  END DO
  Eeff = (/ 0d0, HUGE(Eeff) /)
  eff = 1d0
  CALL DDCalc_SetDetector(Detector,NE=NE,E=E,NEeff=2,Eeff=Eeff,Neff=0,eff=eff)
  
  ! For high verbosity level, we are printing reference rates, so
  ! set "actual" rates to same.
  IF (VerbosityLevel .GE. 4) CALL DDCalc_SetWIMP(WIMP,sigmaSI=1d0,sigmaSD=1d0)
  
  ! Do rate calculations
  CALL DDCalc_CalcRates(Detector,WIMP,Halo)
  
  ! Write out header.
  IF (VerbosityLevel .GE. 2) THEN
    CALL WriteCommandHeader()
    CALL WriteWIMPHeader(WIMP)
    CALL WriteHaloHeader(Halo)
    CALL WriteDetectorHeader(Detector)
    CALL WriteSpectrumHeader()
  END IF
  
  IF (VerbosityLevel .GE. 1) THEN
    CALL WriteSpectrumColumnHeader()
  END IF
  
  ! Write out table.
  CALL WriteSpectrumData(Detector)
  
END SUBROUTINE


! ----------------------------------------------------------------------
! Main routine to calculate events as a function of mass.
! Given for fixed WIMP-nucleon cross-sections.
! 
SUBROUTINE DDCalc_MainEventsByMass()

  IMPLICIT NONE
  TYPE(DetectorStruct) :: Detector
  TYPE(WIMPStruct) :: WIMP
  TYPE(HaloStruct) :: Halo
  LOGICAL :: use_log
  INTEGER :: Nm,I
  REAL*8 :: m,mmin,mmax
  REAL*8 :: sigmapSI,sigmanSI,sigmapSD,sigmanSD,GpSI,GnSI,GpSD,GnSD
  TYPE(TabulationStruct) :: TS
  
  ! Show usage and exit
  IF (ShowUsageRequested()) THEN
    CALL ShowUsage()
    STOP
  END IF
  
  ! Initializes detector, halo and WIMP structures.
  Detector = InitializeCommandLine(WIMP,Halo)
  
  ! Determine WIMP mass tabulation
  mmin    =    1d0
  mmax    = 1000d0
  Nm      =  -20
  use_log = .TRUE.
  CALL GetTabulationArgs('m-tabulation',mmin,mmax,Nm,use_log)
  CALL InitTabulation(TS,mmin,mmax,Nm,use_log)
  
  ! Get fixed cross-sections; will need to be reset at each mass.
  ! Set to negative if coupling is negative (this is the intended
  ! meaning for negative cross-sections used as input).
  CALL DDCalc_GetWIMP(WIMP,sigmapSI=sigmapSI,sigmanSI=sigmanSI, &
                      sigmapSD=sigmapSD,sigmanSD=sigmanSD,GpSI=GpSI, &
                      GnSI=GnSI,GpSD=GpSD,GnSD=GnSD)
  IF (GpSI .LT. 0d0) sigmapSI = -ABS(sigmapSI)
  IF (GnSI .LT. 0d0) sigmanSI = -ABS(sigmanSI)
  IF (GpSD .LT. 0d0) sigmapSD = -ABS(sigmapSD)
  IF (GnSD .LT. 0d0) sigmanSD = -ABS(sigmanSD)
  
  ! Write out header.
  IF (VerbosityLevel .GE. 2) THEN
    CALL WriteCommandHeader()
    CALL WriteHaloHeader(Halo)
    CALL WriteDetectorHeader(Detector)
    CALL WriteEventsByMassHeader(WIMP)
  END IF
  
  IF (VerbosityLevel .GE. 1) THEN
    CALL WriteEventsByMassColumnHeader(Detector)
  END IF
  
  ! Cycle over masses
  DO I = 0,TS%N
    ! Set WIMP mass and cross-sections
    m = TabulationValue(TS,I)
    CALL DDCalc_SetWIMP(WIMP,m=m,sigmapSI=sigmapSI,sigmanSI=sigmanSI,sigmapSD=sigmapSD,sigmanSD=sigmanSD)
    ! Do rate calculations
    CALL DDCalc_CalcRates(Detector,WIMP,Halo)
    ! Write out table data line.
    CALL WriteEventsByMassData(Detector,WIMP)
  END DO
  
END SUBROUTINE


! ----------------------------------------------------------------------
! Main routine to calculate likelihood constraints as a function of
! mass for SI couplings.
! 
SUBROUTINE DDCalc_MainConstraintsSI()

  IMPLICIT NONE
  TYPE(DetectorStruct) :: Detector
  TYPE(WIMPStruct) :: WIMP
  TYPE(HaloStruct) :: Halo
  LOGICAL :: use_log
  INTEGER :: Nm,I
  REAL*8 :: m,mmin,mmax,x,lnp,thetaG,Gp,Gn,s1,s2
  TYPE(TabulationStruct) :: TS
  
  ! Show usage and exit
  IF (ShowUsageRequested()) THEN
    CALL ShowUsage()
    STOP
  END IF
  
  ! Initializes detector, halo and WIMP structures.
  ! Do not need sub-intervals as we only use the total events.
  Detector = InitializeCommandLine(WIMP,Halo,intervals=.FALSE.)
  
  ! Determine WIMP mass tabulation
  mmin    =    1d0
  mmax    = 1000d0
  Nm      =  -20
  use_log = .TRUE.
  CALL GetTabulationArgs('m-tabulation',mmin,mmax,Nm,use_log)
  CALL InitTabulation(TS,mmin,mmax,Nm,use_log)
  
  ! Get p-value for exclusion limit
  CALL ReadLogPValue(lnp)
  
  ! Get angle of (Gp,Gn), which will be kept fixed.
  IF (GetLongArgReal('theta-SI',x)) THEN
    thetaG = x
  ELSE IF (GetLongArgReal('theta-SI-pi',x)) THEN
    thetaG = PI*x
  ELSE
    thetaG = 0.25d0*PI
  END IF
  
  ! Will initially use these couplings at every mass;
  ! constraint calculations will perform appropriate rescaling.
  Gp = 1d0*COS(thetaG)
  Gn = 1d0*SIN(thetaG)
  IF (ABS(Gp) .LT. 1d-8) Gp = 0d0
  IF (ABS(Gn) .LT. 1d-8) Gn = 0d0
  
  ! Calculate allowed signal rates
  CALL FeldmanCousinsPoissonCI(lnp,Detector%Nevents,Detector%MuBackground,s1,s2)
  
  ! Write out header.
  IF (VerbosityLevel .GE. 2) THEN
    CALL WriteCommandHeader()
    CALL WriteHaloHeader(Halo)
    CALL WriteDetectorHeader(Detector)
    CALL WriteConstraintsSIHeader(lnp,thetaG,s1,s2)
  END IF
  
  IF (VerbosityLevel .GE. 1) THEN
    CALL WriteConstraintsSIColumnHeader()
  END IF
  
  ! Cycle over masses
  DO I = 0,TS%N
    ! Set WIMP mass and cross-sections
    m = TabulationValue(TS,I)
    CALL DDCalc_SetWIMP(WIMP,m=m,GpSI=Gp,GnSI=Gn,GpSD=0d0,GnSD=0d0)
    ! Do rate calculations
    CALL DDCalc_CalcRates(Detector,WIMP,Halo)
    ! Write out table data line.
    CALL WriteConstraintsSIData(s1,s2,Detector,WIMP)
  END DO
  
END SUBROUTINE


! ----------------------------------------------------------------------
! Main routine to calculate likelihood constraints as a function of
! mass for SD couplings.
! 
SUBROUTINE DDCalc_MainConstraintsSD()

  IMPLICIT NONE
  TYPE(DetectorStruct) :: Detector
  TYPE(WIMPStruct) :: WIMP
  TYPE(HaloStruct) :: Halo
  LOGICAL :: use_log
  INTEGER :: Nm,I
  REAL*8 :: m,mmin,mmax,x,lnp,thetaG,Gp,Gn,s1,s2
  TYPE(TabulationStruct) :: TS
  
  ! Show usage and exit
  IF (ShowUsageRequested()) THEN
    CALL ShowUsage()
    STOP
  END IF
  
  ! Initializes detector, halo and WIMP structures.
  ! Do not need sub-intervals as we only use the total events.
  Detector = InitializeCommandLine(WIMP,Halo,intervals=.FALSE.)
  
  ! Determine WIMP mass tabulation
  mmin    =    1d0
  mmax    = 1000d0
  Nm      =  -20
  use_log = .TRUE.
  CALL GetTabulationArgs('m-tabulation',mmin,mmax,Nm,use_log)
  CALL InitTabulation(TS,mmin,mmax,Nm,use_log)
  
  ! Get p-value for exclusion limit
  CALL ReadLogPValue(lnp)
  
  ! Get angle of (Gp,Gn), which will be kept fixed.
  IF (GetLongArgReal('theta-SD',x)) THEN
    thetaG = x
  ELSE IF (GetLongArgReal('theta-SD-pi',x)) THEN
    thetaG = PI*x
  ELSE
    thetaG = 0.25d0*PI
  END IF
  
  ! Will initially use these couplings at every mass;
  ! constraint calculations will perform appropriate rescaling.
  Gp = 1d0*COS(thetaG)
  Gn = 1d0*SIN(thetaG)
  IF (ABS(Gp) .LT. 1d-8) Gp = 0d0
  IF (ABS(Gn) .LT. 1d-8) Gn = 0d0
  
  ! Calculate allowed signal rates
  CALL FeldmanCousinsPoissonCI(lnp,Detector%Nevents,Detector%MuBackground,s1,s2)
  
  ! Write out header.
  IF (VerbosityLevel .GE. 2) THEN
    CALL WriteCommandHeader()
    CALL WriteHaloHeader(Halo)
    CALL WriteDetectorHeader(Detector)
    CALL WriteConstraintsSDHeader(lnp,thetaG,s1,s2)
  END IF
  
  IF (VerbosityLevel .GE. 1) THEN
    CALL WriteConstraintsSDColumnHeader()
  END IF
  
  ! Cycle over masses
  DO I = 0,TS%N
    ! Set WIMP mass and cross-sections
    m = TabulationValue(TS,I)
    CALL DDCalc_SetWIMP(WIMP,m=m,GpSI=0d0,GnSI=0d0,GpSD=Gp,GnSD=Gn)
    ! Do rate calculations
    CALL DDCalc_CalcRates(Detector,WIMP,Halo)
    ! Write out table data line.
    CALL WriteConstraintsSDData(s1,s2,Detector,WIMP)
  END DO
  
END SUBROUTINE


! ----------------------------------------------------------------------
! Main routine to calculate exclusion limits as a function of mass for
! SI couplings.
! 
SUBROUTINE DDCalc_MainLimitsSI()

  IMPLICIT NONE
  TYPE(DetectorStruct) :: Detector
  TYPE(WIMPStruct) :: WIMP
  TYPE(HaloStruct) :: Halo
  LOGICAL :: use_log
  INTEGER :: Nm,I
  REAL*8 :: m,mmin,mmax,x,lnp,thetaG,Gp,Gn
  TYPE(TabulationStruct) :: TS
  
  ! Show usage and exit
  IF (ShowUsageRequested()) THEN
    CALL ShowUsage()
    STOP
  END IF
  
  ! Initializes detector, halo and WIMP structures.
  ! Will calculate rates for sub-intervals as these are used by
  ! maximum gap method.
  Detector = InitializeCommandLine(WIMP,Halo,intervals=.TRUE.)
  
  ! Determine WIMP mass tabulation
  mmin    =    1d0
  mmax    = 1000d0
  Nm      =  -20
  use_log = .TRUE.
  CALL GetTabulationArgs('m-tabulation',mmin,mmax,Nm,use_log)
  CALL InitTabulation(TS,mmin,mmax,Nm,use_log)
  
  ! Get p-value for exclusion limit
  CALL ReadLogPValue(lnp)
  
  ! Get angle of (Gp,Gn), which will be kept fixed.
  IF (GetLongArgReal('theta-SI',x)) THEN
    thetaG = x
  ELSE IF (GetLongArgReal('theta-SI-pi',x)) THEN
    thetaG = PI*x
  ELSE
    thetaG = 0.25d0*PI
  END IF
  
  ! Will initially use these couplings at every mass;
  ! limit calculations will perform appropriate rescaling.
  Gp = 1d0*COS(thetaG)
  Gn = 1d0*SIN(thetaG)
  IF (ABS(Gp) .LT. 1d-8) Gp = 0d0
  IF (ABS(Gn) .LT. 1d-8) Gn = 0d0
  
  ! Write out header.
  IF (VerbosityLevel .GE. 2) THEN
    CALL WriteCommandHeader()
    CALL WriteHaloHeader(Halo)
    CALL WriteDetectorHeader(Detector)
    CALL WriteLimitsSIHeader(lnp,thetaG,Detector)
  END IF
  
  IF (VerbosityLevel .GE. 1) THEN
    CALL WriteLimitsSIColumnHeader()
  END IF
  
  ! Cycle over masses
  DO I = 0,TS%N
    ! Set WIMP mass and cross-sections
    m = TabulationValue(TS,I)
    CALL DDCalc_SetWIMP(WIMP,m=m,GpSI=Gp,GnSI=Gn,GpSD=0d0,GnSD=0d0)
    ! Do rate calculations
    CALL DDCalc_CalcRates(Detector,WIMP,Halo)
    ! Write out table data line.
    CALL WriteLimitsSIData(lnp,Detector,WIMP)
  END DO
  
END SUBROUTINE


! ----------------------------------------------------------------------
! Main routine to calculate exclusion limits as a function of mass for
! SD couplings.
! 
SUBROUTINE DDCalc_MainLimitsSD()

  IMPLICIT NONE
  TYPE(DetectorStruct) :: Detector
  TYPE(WIMPStruct) :: WIMP
  TYPE(HaloStruct) :: Halo
  LOGICAL :: use_log
  INTEGER :: Nm,I
  REAL*8 :: m,mmin,mmax,x,lnp,thetaG,Gp,Gn
  TYPE(TabulationStruct) :: TS
  
  ! Show usage and exit
  IF (ShowUsageRequested()) THEN
    CALL ShowUsage()
    STOP
  END IF
  
  ! Initializes detector, halo and WIMP structures.
  ! Will calculate rates for sub-intervals as these are used by
  ! maximum gap method.
  Detector = InitializeCommandLine(WIMP,Halo,intervals=.TRUE.)
  
  ! Determine WIMP mass tabulation
  mmin    =    1d0
  mmax    = 1000d0
  Nm      =  -20
  use_log = .TRUE.
  CALL GetTabulationArgs('m-tabulation',mmin,mmax,Nm,use_log)
  CALL InitTabulation(TS,mmin,mmax,Nm,use_log)
  
  ! Get p-value for exclusion limit
  CALL ReadLogPValue(lnp)
  
  ! Get angle of (Gp,Gn), which will be kept fixed.
  IF (GetLongArgReal('theta-SD',x)) THEN
    thetaG = x
  ELSE IF (GetLongArgReal('theta-SD-pi',x)) THEN
    thetaG = PI*x
  ELSE
    thetaG = 0.25d0*PI
  END IF
  
  ! Will initially use these couplings at every mass;
  ! limit calculations will perform appropriate rescaling.
  Gp = 1d0*COS(thetaG)
  Gn = 1d0*SIN(thetaG)
  IF (ABS(Gp) .LT. 1d-8) Gp = 0d0
  IF (ABS(Gn) .LT. 1d-8) Gn = 0d0
  
  ! Write out header.
  IF (VerbosityLevel .GE. 2) THEN
    CALL WriteCommandHeader()
    CALL WriteHaloHeader(Halo)
    CALL WriteDetectorHeader(Detector)
    CALL WriteLimitsSDHeader(lnp,thetaG,Detector)
  END IF
  
  IF (VerbosityLevel .GE. 1) THEN
    CALL WriteLimitsSDColumnHeader()
  END IF
  
  ! Cycle over masses
  DO I = 0,TS%N
    ! Set WIMP mass and cross-sections
    m = TabulationValue(TS,I)
    CALL DDCalc_SetWIMP(WIMP,m=m,GpSI=0d0,GnSI=0d0,GpSD=Gp,GnSD=Gn)
    ! Do rate calculations
    CALL DDCalc_CalcRates(Detector,WIMP,Halo)
    ! Write out table data line.
    CALL WriteLimitsSDData(lnp,Detector,WIMP)
  END DO
  
END SUBROUTINE



!=======================================================================
! INITIALIZATION
!=======================================================================

! ----------------------------------------------------------------------
! Initialization routine that sets various parameters to default values
! or values specified on the command line.  This must be called before
! using the module (unless one of the other initialization routines is
! called).
! 
! This version of the routine uses command-line options to set up the
! various structures.  It should _only_ be used in programs that intend
! to use such command-line options for setup; otherwise, the
! Initialize() routine above should be used instead.
! 
! Optional input arguments:
!   cmdline     Specifies if the command-line should be checked for
!               values (default: false).  Should be set false if this
!               module is used by external programs.
!   intervals   Indicates if rates for any sub-intervals available from
!               the efficiencies file should be calculated as well as
!               the total rate (default: true).  These are unnecessary
!               for some calculations.
! 
FUNCTION InitializeCommandLine(WIMP, Halo, intervals) RESULT(Detector)

  IMPLICIT NONE
  TYPE(WIMPStruct), INTENT(INOUT) :: WIMP
  TYPE(HaloStruct), INTENT(INOUT) :: Halo  
  LOGICAL, INTENT(IN), OPTIONAL :: intervals
  LOGICAL :: intervals0
  TYPE(ArgumentStruct) :: Arguments
  TYPE(DetectorStruct) :: Detector
  
  intervals0 = .TRUE.
  IF (PRESENT(intervals)) intervals0 = intervals
  
  ! Extract parameters and options from command-line.
  Arguments = ReadArguments()
  
  ! Set the output verbosity level.
  CALL InitVerbosity()
  
  ! Initialize WIMP mass and couplings.
  ! Requires ReadArguments to have been called.
  WIMP = DDCalc_InitWIMPCommandLine(Arguments)
  
  ! Initialize halo.
  Halo = DDCalc_InitHaloCommandLine()
  
  ! Initialize detector isotopes, efficiencies, array sizing, etc.
  Detector = DDCalc_ChooseAnalysisCommandLine(intervals=intervals0)
  CALL DDCalc_InitDetectorDataCommandLine(Detector)
  
END FUNCTION


! ----------------------------------------------------------------------
! Reads the command line arguments and separates them into options
! and parameters.  See the definition of the ArgumentStruct in
! DDTypes for details.
! 
FUNCTION ReadArguments() result(Arguments)

  IMPLICIT NONE
  TYPE(ArgumentStruct) :: Arguments
  LOGICAL :: L
  CHARACTER*2 :: firstchars
  INTEGER :: Narg,Nopt,Nparam,I,Iopt,Iparam,ios
  REAL*8 :: x
  
  Narg   = IARGC()
  Nopt   = 0
  Nparam = 0
  
  ! NOTE: cannot use flags of form -<flag> as negative numbers
  ! will be improperly parsed
  
  ! Count argument types
  DO I=1,Narg
    CALL GETARG(I,firstchars)
    IF (firstchars .EQ. '--') THEN
      Nopt = Nopt + 1
    ELSE
      Nparam = Nparam + 1
    END IF
  END DO
  
  Arguments%Noptions = Nopt
  IF (ALLOCATED(Arguments%options)) DEALLOCATE(Arguments%options)
  ALLOCATE(Arguments%options(1:Nopt))
  
  Arguments%Nparameters = Nparam
  IF (ALLOCATED(Arguments%parameters)) DEALLOCATE(Arguments%parameters)
  ALLOCATE(Arguments%parameters(1:Nparam))
  IF (ALLOCATED(Arguments%values)) DEALLOCATE(Arguments%values)
  ALLOCATE(Arguments%values(1:Nparam))
  
  ! Divide arguments up
  Iopt   = 0
  Iparam = 0
  DO I=1,Narg
    CALL GETARG(I,firstchars)
    IF (firstchars .EQ. '--') THEN
      Iopt = Iopt + 1
      CALL GETARG(I,Arguments%options(Iopt))
    ELSE
      Iparam = Iparam + 1
      CALL GETARG(I,Arguments%parameters(Iparam))
    END IF
  END DO
  
  ! Try to parse parameters to floating point or logical -> floating
  ! point (T=1,F=0).  Set to -HUGE(1d0) if cannot be parsed.
  DO I=1,Nparam
    READ(UNIT=Arguments%parameters(I),FMT=*,IOSTAT=ios) x
    IF (ios .EQ. 0) THEN
      Arguments%values(I) = x
    ELSE
      READ(UNIT=Arguments%parameters(I),FMT=*,IOSTAT=ios) L
      IF (ios .EQ. 0) THEN
        IF (L) THEN
          Arguments%values(I) = 1d0
        ELSE
          Arguments%values(I) = 0d0
        END IF
      ELSE
        Arguments%values(I) = -HUGE(1d0)
      END IF
    END IF
  END DO
  
END FUNCTION


! ----------------------------------------------------------------------
! Initializes the verbosity level using command line arguments or
! default values.
! 
! Possible options:
!   --verbosity=<value>  ! Level of verbosity (default is 1)
!   --verbose            ! Equivalent to --verbosity=2
!   --quiet              ! Equivalent to --verbosity=0
! 
SUBROUTINE InitVerbosity()
  IMPLICIT NONE
  INTEGER :: verb
  
  IF (GetLongArgInteger('verbosity',verb)) THEN
    VerbosityLevel = verb
  ELSE IF (GetLongArg('verbose')) THEN
    VerbosityLevel = 2
  ELSE IF (GetLongArg('quiet')) THEN
    VerbosityLevel = 0
  ELSE
    VerbosityLevel = 1
  END IF
  
END SUBROUTINE


! ----------------------------------------------------------------------
! Reads the p-value or CL from command line arguments or sets them to
! a default value [p=0.1 or 90% CL].  We use p+CL = 1 and return only
! the p-value (as the logarithm of its value).
! 
! Possible options:
!   --p-value=<val>         ! The p-value
!   --log-p-value=<val>     ! The logarithm of the p-value
!   --p-value-sigma=<val>   ! The p-value corresponding to the given number of s.d.'s
!                           ! in the normal distribution
!  --confidence-level=<val> ! The confidence level (1-p)
!  --confidence-level-sigma=<val>
!                           ! Equivalent to --p-value-sigma
! 
! Output argument:
!   lnp             The logarithm of the p-value (CL = 1-p)
! 
SUBROUTINE ReadLogPValue(lnp)
  IMPLICIT NONE
  REAL*8, INTENT(OUT) :: lnp
  LOGICAL :: calc_nsigma
  REAL*8 :: x,y,p,nsigma
  
  calc_nsigma = .FALSE.
  
  ! Default
  lnp = LOG(0.1d0)
  
  ! Process arguments
  IF (GetLongArgReal('confidence-level',x)) THEN
    lnp = LOG(MAX(1-x,TINY(1d0)))
  ELSE IF (GetLongArgReal('p-value',x)) THEN
    lnp = LOG(MAX(x,TINY(1d0)))
  ELSE IF (GetLongArgReal('log-p-value',x)) THEN
    lnp = x
  ELSE IF (GetLongArgReal('p-value-sigma',x)) THEN
    calc_nsigma = .TRUE.
    nsigma = x
  ELSE IF (GetLongArgReal('confidence-level-sigma',x)) THEN
    calc_nsigma = .TRUE.
    nsigma = x
  END IF
  
  ! calculate p-value in terms of normal distribution at nsigma s.d.'s.
  ! In that case:
  !    p = erfc(nsigma/\sqrt{2})
  IF (calc_nsigma) THEN
    IF (x .LT. 25d0) THEN
      lnp = LOG(ERFC(x/SQRT2))
    ELSE
      ! For large nsigma, use asymptotic expansion of erfc.
      ! Unnecessarily calculated to near full double precision....
      y = 1d0 / x**2
      lnp = -0.5d0*x**2 - LOG(SQRTPI*x/SQRT2) + LOGp1(-y*(1-3*y*(1-5*y*(1-7*y*(1-9*y*(1-11*y))))))
    END IF
  END IF
  
END SUBROUTINE




END MODULE


 
