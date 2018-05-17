PROGRAM DDCalc_exclusionF

  USE DDCALC
  USE DDExperiments

  IMPLICIT NONE
  TYPE(DetectorStruct) :: Detector
  TYPE(WIMPStruct) :: WIMP
  TYPE(HaloStruct) :: Halo
  REAL*8 :: BGlogL,mDM,sigmatest,limit,logLlimit,mDMmin,mDMmax
  INTEGER :: mDMsteps, mDMi

! Define the smallest and largest DM mass to be considered, as well as the number of intermediate values
  mDMmin = 1.d0
  mDMmax = 1.d4
  mDMsteps = 80

! Define the (arbitrary) reference cross section (in pb) for which event rates will be calculated
  sigmatest = 1.d-8

! Initialize WIMP and DM Halo object with default values
  WIMP = DDCalc_InitWIMP()
  Halo = DDCalc_InitHalo()

! Optional: Change the parameters of the Standard Halo Model to correspond to standard conventions
  CALL DDCalc_SetHalo(Halo,rho=0.3d0,vrot=220.d0,v0=220.d0)
! Optional: Use a tabulated velocity integral instead of the Standard Halo Model one
!  CALL DDCalc_SetHalo(Halo,rho=0.3d0,g_file='Halos/gvmin_VL2_shell.txt',g_column=5,Nvmin=1000)

  Detector = XENON1T_2017_Init()

! Step 1: Calculate the log likelihood for the background-only hypothesis
! This can be achieved by setting the DM couplings equal to zero (in which case the DM mass is irrelevant)
  CALL DDCalc_SetWIMP(WIMP,m=1d0,DMtype='SIonly',params=[0.d0,0.d0])
  CALL DDCalc_CalcRates(Detector, WIMP, Halo)
  BGlogL = DDCalc_LogLikelihood(Detector)

! Step 2: Calculate the critical value of the log likelihood that corresponds to the exclusion limit
! Here we assume that -2 delta log L follows a 1/2 chi^2 distribution with 1 degree of freedom
! The one-sided upper bound at 90% confidence level is then given by -2 delta log L = 1.64
  limit = 1.64
  logLlimit = BGlogL - (limit / 2.d0)

! Note: For an experiment seeing a slight excess (such that the best-fit point is different from the background-only hypothesis),
! BGlogL should be replaced by the log likelihood of the best-fit point, which must be determined from a 2d parameter scan.
! See arXiv:1407.6617 for details.

! Step 3: Calculate the spin-independent exclusion limit (assuming equal couplings to protons and neutrons)
  WRITE (*,*) 'Bound on the spin-independent cross section (in pb) as a function of the DM mass (in GeV)'
  WRITE (*,*) 'for equal couplings to protons and neutrons:'
  DO mDMi = 0,mDMsteps
    mDM = mDMmin * (mDMmax/mDMmin)**(REAL(mDMi)/mDMsteps)
    CALL DDCalc_SetWIMP_msigma(WIMP, mDM, sigmatest, sigmatest, 0d0, 0d0)
    CALL DDCalc_CalcRates(Detector, WIMP, Halo)
    WRITE (*,*) mDM,DDCalc_ScaleToPValue(Detector,logLlimit)*sigmatest
  END DO

! Step 4: Calculate spin-dependent exclusion limit (assuming couplings only to protons)
  WRITE (*,*) 'Bound on the spin-dependent cross section (in pb) as a function of the DM mass (in GeV)'
  WRITE (*,*) 'for couplings to protons only:'
  DO mDMi = 0,mDMsteps
    mDM = mDMmin * (mDMmax/mDMmin)**(REAL(mDMi)/mDMsteps)
    CALL DDCalc_SetWIMP_msigma(WIMP, mDM, 0d0, 0d0, sigmatest, 0d0)
    CALL DDCalc_CalcRates(Detector, WIMP, Halo)
    WRITE (*,*) mDM,DDCalc_ScaleToPValue(Detector,logLlimit)*sigmatest
  END DO

END PROGRAM

