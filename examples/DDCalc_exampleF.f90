!#######################################################################
! DDCALC EXAMPLE PROGRAM (FORTRAN)
! This program shows how to use the DDCalc module for calculating
! various direct detection constraints.
! 
! Run:
!   ./DDCalc_exampleF [--mfa|--mG]
! where the optional flag specifies the form in which the WIMP-nucleon
! couplings will be provided (default: --mfa).
! 
!   Created by
!   Chris Savage    Nordita/Utah/SavageArchitectures
!   Pat Scott       Imperial College
!   Martin White    Adelaide
!   Felix Kahlhoefer   DESY
!   Sebastian Wild     DESY
!   ddcalc@projects.hepforge.org
! 
!#######################################################################

PROGRAM DDCalc_exampleF

  ! The module must be loaded via a 'USE DDCalc' statement
  ! in any routine that is to make use of DDCalc routines.
  ! To use anything but the default analysis also requires
  ! explicitly loading the DDExperiments module.
  USE DDCalc
  USE DDExperiments

  IMPLICIT NONE
  CHARACTER*32 :: arg
  INTEGER :: I,Narg,type
  REAL*8 :: M,xpSI,xnSI,xpSD,xnSD
  REAL*8 :: GpSI,GnSI,GpSD,GnSD,fp,fn,ap,an
  REAL*8 :: lnp

  ! These three types are the bedrock of DDCalc.  (Almost) every calculation
  ! needs to be provided with an instance of each of these to do its job.
  ! You can have as many different instances as you want though, corresponding
  ! to e.g. different detectors/analyses, WIMP models and DM halo models.
  ! (who said you couldn't do OOP in Fortran?)
  TYPE(WIMPStruct) :: WIMP
  TYPE(HaloStruct) :: Halo
  TYPE(DetectorStruct) :: MyDetector, XENON, LUX, SCDMS, SIMPLE

  ! These constants will be used to specify the type of input parameters
  INTEGER, PARAMETER :: TYPE_MG     = 1
  INTEGER, PARAMETER :: TYPE_MFA    = 2
  
  ! Parse command line options for this example program.
  ! In particlular, determine how the WIMP parameters will be specified.
  type = TYPE_MFA
  Narg = IARGC()
  DO I = 1,Narg
    CALL GETARG(I,arg)
    IF (arg .EQ. '--mG') THEN
      type = TYPE_MG
    ELSE IF (arg .EQ. '--mfa') THEN
      type = TYPE_MFA
    ELSE IF (arg .EQ. '--help') THEN
      WRITE(*,*) "Usage:"
      WRITE(*,*) "  ./DDCalc_exampleF [--mG|--mfa]"
      WRITE(*,*) "where the optional flag specifies the form in which the WIMP-"
      WRITE(*,*) "nucleon couplings will be provided (default: --mfa)."
      STOP;
    ELSE
      WRITE(*,*) "WARNING: Ignoring unknown argument '" // TRIM(arg) // "'."
    END IF
  END DO
  
  ! Write out directions for specifying input to this example program.
  ! WriteDescription is defined below.
  CALL WriteDescription(type)
    
  ! Note that we never have to initialise DDCalc as a whole, we just
  ! have to create Detectors, WIMPs and Halos, then manipulate their
  ! parameters and hand them back to DDCalc to do calculations on.

  ! Initialise a DM Halo object to default values.  See below for how to
  ! modify these values.
  Halo = DDCalc_InitHalo()
    
  ! Initialise a WIMP object to default values.  Actually, this isn't
  ! necessary here, as we set the WIMP properties from the commandline
  ! later -- but here's how you would make a default version if needed:
  WIMP = DDCalc_InitWIMP()

  ! As we are responsible adults, we also choose our own detectors below.
  ! Here is what you'd do if you wanted to just rely on the default
  ! (currently LUX 2013):
  MyDetector = DDCalc_InitDetector(.TRUE.)

  ! Explicitly create detector objects for all the experiments to be
  ! used (set up isotopes, efficiencies, array sizing, etc.)  The   
  ! argument indicates if extra sub-interval calculations should
  ! be performed.  Those calculations are required for maximum gap
  ! analyses, but are unnecessary for calculating total rates and
  ! likelihoods.  If .FALSE. is given, a no-background-subtraction
  ! p-value can still be calculated, but a Poisson is used instead
  ! of the maximum gap.  We show some maximum gap results below, so
  ! we must use .TRUE. here (the flag is ignored for experiments
  ! that do not have the event energies necessary for a maximum gap
  ! analysis).
  XENON    = XENON100_2012_Init(.TRUE.)
  LUX      = LUX_2013_Init(.TRUE.)
  SCDMS    = SuperCDMS_2014_Init(.TRUE.)
  SIMPLE   = SIMPLE_2014_Init(.TRUE.)
  
  ! Can optionally specify a minimum recoil energy to be included in
  ! the rate calculations [keV].  Note the efficiency curves already
  ! account for detector and analysis thresholds regardless of this
  ! setting, so setting this to 0 keV (the default behavior when
  ! initialization is performed) does not imply that very low energy
  ! recoils actually contribute to the signal.
  ! EXAMPLE: Uncomment to set a minimum recoil energy of 3 keV for LUX:
  !CALL DDCalc_SetEmin(LUX,3d0)
  
  ! Advanced usage:
  ! The DDCalc_SetDetector() routine overrides aspects of the detector
  ! configuration, using optional arguments that can be specified via
  ! keywords.  All the values used in the examples here are those used
  ! for the LUX case, so these specific calls are not actually necessary.
  ! 
  ! Set element(s) to use by their atomic number, along with the
  ! stoichiometry (assumed to be 1:1 for all elements if not given).
  ! Spin-dependent interactions only implemented for a limited number
  ! of isotopes.
  CALL DDCalc_SetDetector(MyDetector,Nelem=1,Zelem=(/54/),stoich=(/1/))
  ! Give explicit list of isotopes, along with their mass fractions.
  !CALL DDCalc_SetDetector(MyDetector,Niso=7,Ziso=(/.../),Aiso=(/.../),fiso=(/.../))
  ! Change parameters.
  CALL DDCalc_SetDetector(MyDetector,mass=118d0,time=85.3d0,Nevents=1,background=0.64d0)
  ! Load efficiency curves from file.  First column is recoil energy
  ! [keV], the next column with values in [0,1] is the total detection
  ! efficiency.  Optionally (intervals=.TRUE.), additional columns are
  ! taken to be the detection efficiency for each interval between
  ! events.  If these columns are not available or should be ignored,
  ! set intervals=.FALSE.  There should be Nevents+1 intervals (not
  ! including total) if used.
  CALL DDCalc_SetDetector(MyDetector,eff_file=TRIM(DDCALC_DIR)//&
                          '/data/example_efficiencies.dat',   &
                          intervals=.TRUE.)
  ! Set the minimum recoil energy [keV] to include.
  CALL DDCalc_SetDetector(MyDetector,Emin=0d0)
  ! Print to screen what our extra detector analysis is, labeled as
  ! '(special)' in the output table.
  WRITE(*,'(A)') ''
  WRITE(*,'(A)') 'The (special) case below is identical to the LUX 2013 analysis.'
  
  ! TESTING:
  ! For comparison, set MyDetector to represent the standard LUX
  ! case, but with a 3 keV minimum recoil energy.
  !CALL DDCalc_SetDetector(MyDetector,Emin=3d0)
  !WRITE(*,'(A)') ''
  !WRITE(*,'(A)') 'The (special) case below is the LUX 2013 analysis with a 3 keV'
  !WRITE(*,'(A)') 'minimum recoil energy imposed.'
    
  ! Optionally set the Standard Halo Model parameters:
  !   rho     Local dark matter density [GeV/cm^3]
  !   vrot    Local disk rotation speed [km/s]
  !   v0      Maxwell-Boltzmann most probable speed [km/s]
  !   vesc    Galactic escape speed [km/s]
  ! This example uses the default values (and is thus optional).
  !CALL DDCalc_SetSHM(0.4d0,235d0,235d0,550d0)

  ! INPUT LOOP >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! Loop over input to this example program.
  ! GetWIMPParams is defined below.
  DO WHILE (GetWIMPParams(type,M,xpSI,xnSI,xpSD,xnSD))
    
    WRITE(*,*)
    
    ! Set the WIMP parameters.
    ! There are three ways to specify the WIMP-nucleon couplings, with
    ! the WIMP mass [GeV] always the first argument:
    !   * DDCalc_SetWIMP_mfa(WIMP,m,fp,fn,ap,an)
    !     The standard couplings fp,fn [GeV^-2] & ap,an [unitless]
    !   * DDCalc_SetWIMP_mG(WIMP,m,GpSI,GnSI,GpSD,GnSD)
    !     The effective 4 fermion vertex couplings GpSI,GnSI,GpSD,GnSD
    !     [GeV^-2], related by:
    !         GpSI = 2 fp        GpSD = 2\sqrt{2} G_F ap
    !         GnSI = 2 fn        GnSD = 2\sqrt{2} G_F an
    ! In the above, 'p' is for proton, 'n' is for neutron, 'SI' is for
    ! spin-independent, and 'SD' is for spin-dependent.
    SELECT CASE(type)
    CASE(TYPE_MG)
      CALL DDCalc_SetWIMP_mG(WIMP,M,xpSI,xnSI,xpSD,xnSD)
    CASE(TYPE_MFA)
      CALL DDCalc_SetWIMP_mfa(WIMP,M,xpSI,xnSI,xpSD,xnSD)
    END SELECT
    
    ! Get the current WIMP parameters with the same signature and units
    ! as above.  The only difference is that WIMP-nucleon cross-sections
    ! are always positive (physical) values.
    CALL DDCalc_GetWIMP_mG(WIMP,M,GpSI,GnSI,GpSD,GnSD)
    CALL DDCalc_GetWIMP_mfa(WIMP,M,fp,fn,ap,an)
    WRITE(*,'(A,1(1X,1PG12.4))') 'WIMP mass [GeV]     ',M
    WRITE(*,*)
    WRITE(*,'(A28,4(1X,A11))')     'WIMP-nucleon couplings          ',  &
        ' proton-SI ',' neutron-SI',' proton-SD ',' neutron-SD'
    WRITE(*,'(A28,4(1X,1PG11.4))') '  G [GeV^-2]                    ',  &
        GpSI,GnSI,GpSD,GnSD
    WRITE(*,'(A28,4(1X,1PG11.4))') '  f & a [GeV^-2,unitless]       ',  &
        fp,fn,ap,an
    WRITE(*,*)
    
    ! Do rate calculations.
    CALL DDCalc_CalcRates(XENON,WIMP,Halo)
    CALL DDCalc_CalcRates(LUX,WIMP,Halo)
    CALL DDCalc_CalcRates(SCDMS,WIMP,Halo)
    CALL DDCalc_CalcRates(SIMPLE,WIMP,Halo)
    CALL DDCalc_CalcRates(MyDetector,WIMP,Halo)
    
    ! Header
    WRITE(*,'(A20,5(2X,A11))') '',' XENON 2012',' LUX 2013  ',          &
        'SuCDMS 2014','SIMPLE 2014',' (special) '
    !WRITE(*,'(A20,7(1X,A12))') '','-----------','-----------',          &
    !    '-----------','-----------','-----------','-----------',        &
    !    '-----------'
    
    ! Event quantities.
    ! The observed number of events (INTEGER).
    WRITE(*,'(A20,5(2X,1X,I6,4X))') 'Observed events                 ', &
        DDCalc_Events(XENON), &
        DDCalc_Events(LUX), &
        DDCalc_Events(SCDMS), &
        DDCalc_Events(SIMPLE), &
        DDCalc_Events(MyDetector)
    ! The average expected background.
    WRITE(*,'(A20,5(2X,1PG11.4))')  'Expected background             ', &
        DDCalc_Background(XENON), &
        DDCalc_Background(LUX), &
        DDCalc_Background(SCDMS), &
        DDCalc_Background(SIMPLE), &
        DDCalc_Background(MyDetector)
    ! The average expected WIMP signal.
    WRITE(*,'(A20,5(2X,1PG11.4))')  'Expected signal                 ', &
        DDCalc_Signal(XENON), &
        DDCalc_Signal(LUX), &
        DDCalc_Signal(SCDMS), &
        DDCalc_Signal(SIMPLE), &
        DDCalc_Signal(MyDetector)
    
    ! The log-likelihoods for the current WIMP; note these are _not_
    ! multiplied by -2.  The likelihood is calculated using a Poisson
    ! given the observed number of events and expected signal+background.
    WRITE(*,'(A20,5(2X,1PG11.4))')  'Log-likelihood                  ', &
        DDCalc_LogLikelihood(XENON), &
        DDCalc_LogLikelihood(LUX), &
        DDCalc_LogLikelihood(SCDMS), &
        DDCalc_LogLikelihood(SIMPLE), &
        DDCalc_LogLikelihood(MyDetector)
    
    ! The logarithm of the p-value, calculated without background
    ! subtraction, using either the maximum gap statistic or a Poisson
    ! statistic, depending on how the detector was initialized.  Note
    ! that this is actually a conservative upper _bound_ on the p-value
    ! in the event of an unknown background and is useful for excluding
    ! WIMP parameters.  However, since it is not a true p-value, it
    ! should not be interpreted as being related to any particular
    ! likelihood.
    WRITE(*,'(A20,5(2X,1PG11.4))')  'Max gap log(p-value)            ', &
        DDCalc_LogPValue(XENON), &
        DDCalc_LogPValue(LUX), &
        DDCalc_LogPValue(SCDMS), &
        DDCalc_LogPValue(SIMPLE), &
        DDCalc_LogPValue(MyDetector)
    
    ! The factor x by which the current WIMP cross- sections must be
    ! multiplied (sigma -> x*sigma, applied to all four WIMP-nucleon
    ! cross-sections) to achieve the given p-value (specified by its
    ! logarithm).  Useful for finding the no-background-subtraction
    ! exclusion limits.  For example, if setWIMP_msigma(100d0,10d0,10d0,
    ! 0d0,0d0) is called, then x*(10. pb) would be the SI cross-section
    ! at a WIMP mass of 100 GeV at which the experiment is excluded at
    ! the 90% CL (p=1-CL).
    lnp = LOG(0.1d0)
    WRITE(*,'(A20,5(2X,1PG11.4))')  'Max gap x for 90% CL            ', &
        DDCalc_ScaleToPValue(XENON), &
        DDCalc_ScaleToPValue(LUX), &
        DDCalc_ScaleToPValue(SCDMS), &
        DDCalc_ScaleToPValue(SIMPLE), &
        DDCalc_ScaleToPValue(MyDetector)
    WRITE(*,'(A60)')  '  * Factor x such that sigma->x*sigma gives desired p-value'
    
    !WRITE(*,*)
    
  END DO  ! END INPUT LOOP <<<<<<<<<<<<<<<<<<<<<
  
  
  CONTAINS

  ! --------------------------------------------
  ! Write a description of how input parameters should be specified
  SUBROUTINE WriteDescription(type)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: type 
    WRITE(*,'(A)') ''
    WRITE(*,'(A)') 'Enter WIMP parameters below.  Only the first two are necessary.'
    WRITE(*,'(A)') 'A blank line terminates input.  The parameters are:'
    WRITE(*,'(A)') ''
    SELECT CASE(type)
    CASE(TYPE_MG)
      WRITE(*,'(A)') '  M     WIMP mass [GeV]'
      WRITE(*,'(A)') '  GpSI  Spin-independent WIMP-proton effective coupling [GeV^-2]'
      WRITE(*,'(A)') '  GnSI  Spin-independent WIMP-neutron effective coupling [GeV^-2]'
      WRITE(*,'(A)') '  GpSD  Spin-dependent WIMP-proton effective coupling [GeV^-2]'
      WRITE(*,'(A)') '  GnSD  Spin-dependent WIMP-neutron effective coupling [GeV^-2]'
    CASE(TYPE_MFA)
      WRITE(*,'(A)') '  M     WIMP mass [GeV]'
      WRITE(*,'(A)') '  fp    Spin-independent WIMP-proton effective coupling [GeV^-2]'
      WRITE(*,'(A)') '  fn    Spin-independent WIMP-neutron effective coupling [GeV^-2]'
      WRITE(*,'(A)') '  ap    Spin-dependent WIMP-proton effective coupling [unitless]'
      WRITE(*,'(A)') '  an    Spin-dependent WIMP-neutron effective coupling [unitless]'
    END SELECT
    !WRITE(*,*) ''
  END SUBROUTINE


  ! --------------------------------------------
  ! Read WIMP parameters from standard input
  FUNCTION GetWIMPParams(type,M,xpSI,xnSI,xpSD,xnSD) RESULT(valid)
    IMPLICIT NONE
    LOGICAL :: valid
    INTEGER, INTENT(IN) :: type
    REAL*8, INTENT(OUT) :: M,xpSI,xnSI,xpSD,xnSD
    CHARACTER*256 :: line
    INTEGER :: I,ios
    REAL*8 :: x(5)
    
    valid = .FALSE.
    
    WRITE(*,'(A)') ''
    WRITE(*,'(A)') '------------------------------------------------------------'
    SELECT CASE(type)
    CASE(TYPE_MG)
      WRITE(*,'(A)') 'Enter values <M GpSI GnSI GpSD GnSD>:'
    CASE(TYPE_MFA)
      WRITE(*,'(A)') 'Enter values <M fp fn ap an>:'
    END SELECT
    READ(*,'(A)') line
    IF (TRIM(line) .EQ. '') RETURN
    
    I = 6
    ios = -1
    DO WHILE ((I .GT. 2) .AND. (ios .NE. 0))
      I = I - 1
      READ(line,*,IOSTAT=ios) x(1:I)
    END DO
    IF (ios .NE. 0) RETURN
    
    valid = .TRUE.
    M    = x(1)
    xpSI = x(2)
    xnSI = x(3)
    xpSD = x(4)
    xnSD = x(5)
    IF (I .LT. 3) xnSI = xpSI
    IF (I .LT. 4) xpSD = 0d0
    IF (I .LT. 5) xnSD = xpSD
    
  END FUNCTION
  
END PROGRAM



!#######################################################################
! END OF FILE
!#######################################################################

