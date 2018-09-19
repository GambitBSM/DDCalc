PROGRAM DDCalc_exclusionF

  USE DDCALC
  USE DDExperiments
  USE DDRates
  USE DDStats

  IMPLICIT NONE
  TYPE(DetectorStruct) :: Detector
  TYPE(WIMPStruct) :: WIMP
  TYPE(HaloStruct) :: Halo
  REAL*8 :: BGlogL,mDM,sigmatest,limit,logLlimit,mDMmin,mDMmax, sigmin, sigmax, sig, logLbest, mDMbest, sigbest, s, s1, s2
  REAL*8, ALLOCATABLE :: logL(:,:)
  INTEGER :: mDMsteps, mDMi, sigsteps, sigi

! Define the smallest and largest DM mass to be considered, as well as the number of intermediate values
  mDMmin = 1.d0
  mDMmax = 1.d4
  mDMsteps = 80

  sigmin = 1.d-11
  sigmax = 1.d-7
  sigsteps = 80

  ALLOCATE(logL(0:mDMsteps,0:sigsteps))

! Define the (arbitrary) reference cross section (in pb) for which event rates will be calculated
  sigmatest = 1.d-10

! Initialize WIMP and DM Halo object with default values
  WIMP = DDCalc_InitWIMP()
  Halo = DDCalc_InitHalo()

! Optional: Change the parameters of the Standard Halo Model to correspond to standard conventions
  CALL DDCalc_SetHalo(Halo,rho=0.3d0,vrot=220.d0,v0=220.d0)
! Optional: Use a tabulated velocity integral instead of the Standard Halo Model one
!  CALL DDCalc_SetHalo(Halo,rho=0.3d0,g_file='Halos/gvmin_VL2_shell.txt',g_column=5,Nvmin=1000)

  Detector = XENON1T_2018_Init()

! *************************************************
! *** Routines for calculating exclusion limits ***
! *************************************************

! Here we assume that the experiment sees no excess, such that the best-fit point is identical to the background-only hypothesis.
! See below for the case of experiments seeing an excess.

! Step 1: Calculate the log likelihood for the background-only hypothesis
! This can be achieved by setting the DM couplings equal to zero (in which case the DM mass is irrelevant)
  CALL DDCalc_SetWIMP(WIMP,m=1d0,DMtype='SIonly',params=[0.d0,0.d0])
  CALL DDCalc_CalcRates(Detector, WIMP, Halo)
  BGlogL = DDCalc_LogLikelihood(Detector)

! Step 2: Calculate the critical value of the log likelihood that corresponds to the exclusion limit
! Here we assume the asymptotic limit, such that -2 delta log L follows a 1/2 chi^2 distribution with 1 degree of freedom
! The one-sided upper bound at 90% confidence level is then given by -2 delta log L = 1.64
  limit = 1.64
  logLlimit = BGlogL - (limit / 2.d0)

! Step 3: Calculate the spin-independent exclusion limit (assuming equal couplings to protons and neutrons)
  WRITE (*,*) ''
  WRITE (*,*) '**************************** Spin-independent limit ****************************'
  WRITE (*,*) ''
  WRITE (*,*) 'Assuming spin-independent scattering with equal couplings to protons and neutrons.'
  WRITE (*,*) ''
  WRITE (*,*) '    log_10(m_DM/GeV)  log_10(sigma_p/pb)'
  WRITE (*,*) ''

  DO mDMi = 0,mDMsteps
    mDM = mDMmin * (mDMmax/mDMmin)**(REAL(mDMi)/mDMsteps)
    CALL DDCalc_SetWIMP_msigma(WIMP, mDM, sigmatest, sigmatest, 0d0, 0d0)
    CALL DDCalc_CalcRates(Detector, WIMP, Halo)
    WRITE (*,'(2(10X,F10.4))') log(mDM)/log(10d0),log(DDCalc_ScaleToPValue(Detector,logLlimit)*sigmatest)/log(10d0)
  END DO

! Step 4: Calculate spin-dependent exclusion limit (assuming couplings only to protons)
  WRITE (*,*) ''
  WRITE (*,*) '***************************** Spin-dependent limit *****************************'
  WRITE (*,*) ''
  WRITE (*,*) 'Assuming spin-dependent scattering with couplings to protons only.'
  WRITE (*,*) ''
  WRITE (*,*) '    log_10(m_DM/GeV)  log_10(sigma_p/pb)'
  WRITE (*,*) ''
  DO mDMi = 0,mDMsteps
    mDM = mDMmin * (mDMmax/mDMmin)**(REAL(mDMi)/mDMsteps)
    CALL DDCalc_SetWIMP_msigma(WIMP, mDM, 0d0, 0d0, sigmatest, 0d0)
    CALL DDCalc_CalcRates(Detector, WIMP, Halo)
    WRITE (*,'(2(10X,F10.4))') log(mDM)/log(10d0),log(DDCalc_ScaleToPValue(Detector,logLlimit)*sigmatest)/log(10d0)
  END DO

! *************************************************
! *** Routines for calculating best-fit regions ***
! *************************************************

! Now we consider a more general case of an experiment seeing an excess,
! i.e. we are interested in constructing confidence intervals / regions.
! See arXiv:1407.6617 for details.

! Step 5: Feldman-Cousins
! If we consider only the total number of expected and observed events, we can employ the Feldman-Cousins method
! to construct a lower and an upper bound on the spin-independent scattering cross section as a function of DM mass

  CALL FeldmanCousinsPoissonCI(-2.3d0,DDCalc_Events(Detector),DDCalc_Background(Detector),s1,s2)
! The first argument is log_10(p), where p = 1 - CL = 0.1 for a 90% confidence limit.
! See src/DDStats.f90 for details.

  WRITE (*,*) ''
  WRITE (*,*) '**************************** Feldmann-Cousins bound ****************************'
  WRITE (*,*) ''
  WRITE (*,*) 'Assuming spin-independent scattering with equal couplings to protons and neutrons.'
  WRITE (*,*) ''
  IF (s1 .GT. 0) THEN
    WRITE (*,*) '      log_10(m_DM/GeV)  log_10(sigma_min/pb)  log_10(sigma_max/pb)'
    WRITE (*,*) ''
  ELSE
    WRITE (*,*) 'Lower bound corresponds to zero cross section. Quoting upper bound only.'
    WRITE (*,*) ''
    WRITE (*,*) '      log_10(m_DM/GeV)  log_10(sigma_max/pb)'
    WRITE (*,*) ''
  END IF

  DO mDMi = 0,mDMsteps
    mDM = mDMmin * (mDMmax/mDMmin)**(REAL(mDMi)/mDMsteps)
    CALL DDCalc_SetWIMP_msigma(WIMP, mDM, sigmatest, sigmatest, 0d0, 0d0)
    CALL DDCalc_CalcRates(Detector, WIMP, Halo)
    s = DDCalc_Signal(Detector)
    IF (s .gt. 0) THEN
      IF (s1 .GT. 0) THEN
        WRITE (*,'(3(12X,F10.4))') log(mDM)/log(10d0),log(s1/s*sigmatest)/log(10d0),log(s2/s*sigmatest)/log(10d0)
      ELSE
        WRITE (*,'(2(12X,F10.4))') log(mDM)/log(10d0),log(s2/s*sigmatest)/log(10d0)
      END IF
    END IF
  END DO


! Step 6: 2d likelihood scan
! To include spectral information, we need to perform a full scan over the 2d parameter space
! We can then determine the combination of mass and cross section that give the best fit to the data and construct
! confidence regions around this point.
  logLbest = BGlogL
  DO mDMi = 0,mDMsteps
    DO sigi = 0,sigsteps
      mDM = mDMmin * (mDMmax/mDMmin)**(REAL(mDMi)/mDMsteps)
      sig = sigmin * (sigmax/sigmin)**(REAL(sigi)/sigsteps)
      CALL DDCalc_SetWIMP_msigma(WIMP, mDM, sig, sig, 0d0, 0d0)
      CALL DDCalc_CalcRates(Detector, WIMP, Halo)
      logL(mDMi,sigi) = DDCalc_LogLikelihood(Detector)
      IF (logL(mDMi,sigi) .GT. logLbest) THEN
        logLbest = logL(mDMi,sigi)
        mDMbest = mDM
        sigbest = sig
      END IF
    END DO
  END DO

  WRITE (*,*) ''
  WRITE (*,*) '****************************** 2d likelihood scan ******************************'
  WRITE (*,*) ''
  WRITE (*,*) 'Assuming spin-independent scattering with equal couplings to protons and neutrons.'
  WRITE (*,*) ''
  WRITE (*,'(A,F7.1,A,G10.3,A,F7.3)') ' Best fit point is m_DM = ',mDM,' GeV, sigma_p = ',sigbest,&
  ' pb with log L = ',logLbest
  WRITE (*,*) ''
  WRITE (*,*) '    log_10(m_DM/GeV)  log_10(sigma_p/pb)               log L      -2 Delta log L'
  WRITE (*,*) ''
  DO mDMi = 0,mDMsteps
    DO sigi = 0,sigsteps
      mDM = mDMmin * (mDMmax/mDMmin)**(REAL(mDMi)/mDMsteps)
      sig = sigmin * (sigmax/sigmin)**(REAL(sigi)/sigsteps)
      WRITE (*,'(4(10X,F10.4))'), log(mDM)/log(10d0),log(sig)/log(10d0),logL(mDMi,sigi),-2.d0*(logL(mDMi,sigi)-logLbest)
    END DO
  END DO

  CALL FeldmanCousinsPoissonCI(-2.3d0,0,0.8d0,s1,s2)

  WRITE (*,*) s1, s2

END PROGRAM

