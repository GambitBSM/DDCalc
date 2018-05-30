!#######################################################################
! DDCALC EXAMPLE PROGRAM (FORTRAN)
! This program shows how to use the DDCalc module for calculating
! various direct detection constraints.
! 
! For various types of DM-nucleon interactions, as well as for various 
! direct detection experiments, it computes the expected signal rates,
! any by comparing to the observed number of events, also the 
! log(likelihood) value for each of the examples.
! In order to convert those likelihood values into an upper bound on
! e.g. the scattering cross section, please have a look at
! DDCalc_exclusionC.cpp.
! 
! Run:
!   ./DDCalc_exampleF 
!
!   Created by
!   Chris Savage    Nordita/Utah/SavageArchitectures
!   Pat Scott       Imperial College
!   Martin White    Adelaide
!   Felix Kahlhoefer   RWTH Aachen
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
  USE DDRates


  IMPLICIT NONE
  INTEGER :: i_bin,type
  REAL*8 :: mDM, sigmap_SI, sigman_SI, sigmap_SD, sigman_SD
  REAL*8 :: DM_spin, fp, fn, ap, an

  ! These three sets of indices refer to instance of the three types that
  ! are the bedrock of DDCalc.  (Almost) every calculation needs to be
  ! provided with an instance of each of these to do its job.  Passing the
  ! index of one of them to DDCalc tells it which one in its internal cache
  ! to use for the calculation requested. You can create as many
  ! different instances as you want, corresponding to e.g. different
  ! detectors/analyses, WIMP models and DM halo models; the factory
  ! funcions create the instances in DDCalc and return you the index of
  ! the resulting object.
  TYPE(WIMPStruct) :: WIMP
  TYPE(HaloStruct) :: Halo
  TYPE(DetectorStruct) :: Detector


  ! Initialise a DM Halo object to default values.
  Halo = DDCalc_InitHalo()
    
  ! Initialise a WIMP objects to default values.
  WIMP = DDCalc_InitWIMP()

  ! Optionally set the Standard Halo Model parameters to values different from the default choices:
  !     rho     Local dark matter density [GeV/cm^3]
  !     vrot    Local disk rotation speed [km/s]
  !     v0      Maxwell-Boltzmann most probable speed [km/s]
  !     vesc    Galactic escape speed [km/s]
  CALL DDCalc_SetSHM(Halo, 0.3d0,235d0,235d0,550d0)

  ! Optional: Use a tabulated velocity integral instead of the Standard Halo Model one
  !  CALL DDCalc_SetHalo(Halo,rho=0.3d0,g_file='Halos/gvmin_VL2_shell.txt',g_column=5,Nvmin=1000)



  ! ****************************************************************************************************************
  ! Example 1: XENON1T (2017) analysis, with standard SI/SD interactions specified
  !            by WIMP-nucleon cross sections. 
  Detector = XENON1T_2017_Init()		! Initalize the XENON1T_2017 detector.
  mDM = 100.0d0                           	! DM Mass in GeV.
  sigmap_SI = 4.0d-9				! SI WIMP-proton cross section in pb.
  sigman_SI = -0.3d-9				! SI WIMP-neutron cross section in pb.
						!   The negative value indicates that the corresponding WIMP-nucleon coupling is negative.
  sigmap_SD = 2.0d-5				! SD WIMP-proton cross section in pb.
  sigman_SD = 8.0d-5				! SD WIMP-neutron cross section in pb.
  CALL DDCalc_SetWIMP_msigma(WIMP, mDM, &
	sigmap_SI, sigman_SI, sigmap_SD, sigman_SD)
  CALL DDCalc_CalcRates(Detector,WIMP,Halo)	! This performs the actual calculation of the rates.

  WRITE (*,*) '******************************************************************'
  WRITE (*,*) 'Example 1: Mixed SI and SD interactions at XENON1T (2017)'
  WRITE (*,*) '           specified by the WIMP-nucleon cross sections.'
  WRITE (*,*) 'mDM       = ', mDM, ' GeV'
  WRITE (*,*) 'sigmap_SI = ', sigmap_SI, ' pb'
  WRITE (*,*) 'sigman_SI = ', sigman_SI, ' pb'
  WRITE (*,*) 'sigmap_SD = ', sigmap_SD, ' pb'
  WRITE (*,*) 'sigman_SD = ', sigman_SD, ' pb'
  WRITE (*,*)
  WRITE (*,*) 'Expected number of signal events:      ', DDCalc_Signal(Detector)
  WRITE (*,*) 'Observed number of events:             ', DDCalc_Events(Detector)
  WRITE (*,*) 'Expected number of background events:  ', DDCalc_Background(Detector)
  WRITE (*,*) 'Log(likelihood):                       ', DDCalc_LogLikelihood(Detector)
  WRITE (*,*) '******************************************************************'
  WRITE (*,*)
  ! ****************************************************************************************************************


  ! ****************************************************************************************************************
  ! Example 2: Example 2: PICO-60 (2017) analysis, with standard SI and SD interactions
  !            specified by the couplings fp, fn, ap, an.          
  !            fp and fn are the coefficients of the DM DM N N term in the Lagrangian, 
  !	       assuming a Majorana DM particle. The definition of ap and an follow the
  !            convention in Jungman&Kamionkowski, 'SUPERSYMMETRIC DARK MATTER'.
  Detector = PICO_60_2017_Init()	! Initalize the PICO_60_2017 detector.
  mDM = 30.0d0                         	! DM Mass in GeV.
  fp = 1.0d-8				! SI WIMP-proton coupling fp, in units [GeV^(-2)].
  fn = 3.0d-8				! SI WIMP-neutron coupling fn, in units [GeV^(-2)].
  ap = 1.0d-2				! SD WIMP-proton coupling ap, unitless.
  an = 0.0d0				! SD WIMP-neutron coupling an, unitless.
  CALL DDCalc_SetWIMP_mfa(WIMP, mDM, fp, fn, ap, an)
  CALL DDCalc_CalcRates(Detector,WIMP,Halo)	! This performs the actual calculation of the rates.

  WRITE (*,*) '******************************************************************'
  WRITE (*,*) 'Example 2: Mixed SI and SD interactions at PICO-60 (2017)'
  WRITE (*,*) '           specified by the WIMP-nucleon couplings.'
  WRITE (*,*) 'mDM  = ', mDM, ' GeV'
  WRITE (*,*) 'fp   = ', fp, ' GeV^(-2)'
  WRITE (*,*) 'fn   = ', fn, ' GeV^(-2)'
  WRITE (*,*) 'ap   = ', ap
  WRITE (*,*) 'an   = ', an
  WRITE (*,*)
  WRITE (*,*) 'Expected number of signal events:      ', DDCalc_Signal(Detector)
  WRITE (*,*) 'Observed number of events:             ', DDCalc_Events(Detector)
  WRITE (*,*) 'Expected number of background events:  ', DDCalc_Background(Detector)
  WRITE (*,*) 'Log(likelihood):                       ', DDCalc_LogLikelihood(Detector)
  WRITE (*,*) '******************************************************************'
  WRITE (*,*)
  ! ****************************************************************************************************************
    
     

  ! ****************************************************************************************************************
  ! Example 3: CRESST (2017) analysis, with a couple of non-relativistic scattering operators                        
  !            set to a non-zero value.  									      
  !            In particular, the analysis of this experiment involves several bins in energy.   
  Detector = CRESST_II_Init()			! Initalize the CRESST_II detector.
  mDM = 3.0d0                           	! DM Mass in GeV.
  DM_spin = 0.5d0				! DM Spin
  CALL DDCalc_SetWIMP_NREffectiveTheory(WIMP, mDM, DM_spin)
	! This defines a WIMP within the non-relativistic effective theory of DM-nucleon interactions,
	! with a given mass in [GeV] and spin.
	! Initially, all the couplings corresponding to the various operators are set to zero

  ! The operator coefficients are set by DDCalc::SetNRCoefficient(WIMP, OpIndex, tau, value).
  !    OpIndex is an integer specifing the operator, ranging from 1 to 18 (excluding 2 and 16).
  !       Additional possible values are -1 and -4, corresponding to the momentum suppressed
  !       operators O_1 * (q^2/mp^2) and O_4 * (q^2/mp^2), respectively.
  !    tau is an integer, with 0 corresponding to the isoscalar part of the operator,
  !       and 1 to the isovector part.
  !    value specifies the operator coefficient corresponding to OpIndex and tau, in units [GeV^(-2)].
  CALL DDCalc_SetNRCoefficient(WIMP, 11, 0, 5.0d-4)
	! This sets the isoscalar part of O_11 to 5.0e-4 GeV^(-2).
  CALL DDCalc_SetNRCoefficient(WIMP, 11, 1, -2.0d-4)
	! This sets the isovector part of O_11 to -2.0e-4 GeV^(-2).
  CALL DDCalc_CalcRates(Detector,WIMP,Halo)	! This performs the actual calculation of the rates.

  WRITE (*,*) '******************************************************************'
  WRITE (*,*) 'Example 3: Non-standard scattering operators at CRESST-II'
  WRITE (*,*) '           This is also an example for an experiment with'
  WRITE (*,*) '           several bins.'
  WRITE (*,*) 'mDM       = ', mDM, ' GeV'
  WRITE (*,*) 'DM spin   = ', DM_spin
  WRITE (*,*) '               Expected signal    Observed    Expected background'
  WRITE (*,'(A, E15.5, A, I5, A, E15.5)') 'All bins    ', DDCalc_Signal(Detector), &
		'      ', &
		DDCalc_Events(Detector), '     ', DDCalc_Background(Detector)

  !WRITE (*,*) DDCalc_Bins(Detector) 
  DO i_bin = 1, DDCalc_Bins(Detector)  
    WRITE (*,'(A, I5, A, E15.5, A, I5, A, E15.5)') 'Bin', i_bin, '    ', &
	DDCalc_BinSignal(Detector, i_bin), '      ', &
	DDCalc_BinEvents(Detector, i_bin), '     ', &
        DDCalc_BinBackground(Detector, i_bin)
  END DO
  !for( int i_bin = 1; i_bin <= DDCalc::Bins(Detector); i_bin = i_bin + 1 ) {
  !   printf("Bin %d\t\t%.5e\t\t%d\t\t%.5e\n", i_bin, DDCalc::BinSignal(Detector, i_bin), DDCalc::BinEvents(Detector, i_bin), DDCalc::BinBackground(Detector, i_bin));
  !   }

  WRITE (*,*)
  WRITE (*,*) 'Log(likelihood):                       ', DDCalc_LogLikelihood(Detector)
  WRITE (*,*) '******************************************************************'
  WRITE (*,*)
  ! ****************************************************************************************************************


  
END PROGRAM



!#######################################################################
! END OF FILE
!#######################################################################

