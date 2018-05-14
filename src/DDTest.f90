PROGRAM DDTest

  USE DDCALC
  USE DDExperiments
  USE DDConstants
  USE DDTypes
  USE DDNumerical
  USE DDWIMP
  USE DDDetectors
  USE DDRates
  USE DDStats
  USE DDHalo
  USE DDNuclear
  USE DDNREffectiveTheory
  USE DDCouplings

  IMPLICIT NONE
  TYPE(DetectorStruct) :: Detector
  TYPE(WIMPStruct) :: WIMP
  TYPE(HaloStruct) :: Halo
  REAL*8 :: BGlogL, fn, fp
  REAL*8 :: params_tmp(45)
  INTEGER :: bins, ibin
  INTEGER :: events
  REAL*8 :: signal,background

  WRITE (*,*) 'This is DDTest.' 
  WIMP = DDCalc_InitWIMP()

  Halo = DDCalc_InitHalo()
  CALL DDCalc_SetHalo(Halo,rho=0.3d0,vrot=235.0d0,v0=235.d0,vesc=550.0d0)


  !!! Example 1: standard SI interaction, specified by params = [fp, fn]
  !CALL DDCalc_SetWIMP(WIMP,m=30d0,DMtype='SIonly',params=[0.7d-8, -0.3d-8])

  !!! Example 2: standard SD interaction, specified by params = [ap, an]
  !CALL DDCalc_SetWIMP(WIMP,m=30d0,DMtype='SDonly',params=[1.0d-1, 5.0d-2])

  !!! Example 3: mixed SI and SD interaction, specified by params = [fp, fn, ap, an]
  !CALL DDCalc_SetWIMP(WIMP,m=30d0,DMtype='SISD',params=[0.7d-8, -0.3d-8, 1.0d-1, 5.0d-2])

  !!! Example 4: Higgsportal interaction (scalar + pseudoscalar term), specified by fsp,fsn,app,apn
  !CALL DDCalc_SetWIMP_Higgsportal(WIMP,30d0, 0.7d-8, -0.3d-8, 1.5d-6, -3d-6)

  !!! Example 5: general non-relativistic effective theory, specified by setting some coefficients of
  !!!            the operators to non-zero values
  !params_tmp = NRET_CreateCoeffList() ! this creates an empty list of coefficients
  !CALL NRET_SetDMSpin(params_tmp, 0.5d0) ! this sets the DM spin to 1/2
  !CALL NRET_SetNRCoefficient(params_tmp, 'Op6', 0, 1.0d-3)  ! this sets the isoscalar operator 6 to a given value
  !CALL NRET_SetNRCoefficient(params_tmp, 'Op6', 1, -0.5d-3) ! this sets the isovector operator 6 to a given value
!  CALL DDCalc_SetWIMP(WIMP,m=30.0d0,DMtype='NREffectiveTheory',params=params_tmp)




  fp = SigmapSItoFp(100.0d0,1.0d-9)
  fn = SigmanSItoFn(100.0d0,0.5d-9)
  WRITE (*,*) 'fp, fn = ', fp, ', ', fn
  CALL DDCalc_SetWIMP(WIMP,m=100.0d0,DMtype='SIonly',params=[fp,fn])

  Detector = Xenon1T_2017_Init()
  CALL DDCalc_CalcRates(Detector, WIMP, Halo)


  WRITE (*,*) 'InitSuccess =', Detector%InitSuccess
  WRITE (*,*) 'StatisticFlag =', Detector%StatisticFlag
  WRITE (*,*) 'exposure =', Detector%exposure
  WRITE (*,*) 'Nevents =', Detector%Nevents
  WRITE (*,*) 'Backgr =', Detector%Backgr
  WRITE (*,*) 'MuSignal =', Detector%MuSignal
  WRITE (*,*) 'LogLikelihood =', DDCalc_LogLikelihood(Detector)

  bins = DDCalc_Bins(Detector)

  WRITE (*,*) 'Nbins =', bins

  DO ibin = 0,bins
    events = DDCalc_BinEvents(Detector, ibin)
    background = DDCalc_BinBackground(Detector, ibin)
    signal = DDCalc_BinSignal(Detector, ibin)
    WRITE (*,*) 'Observed: ',events,' Background: ',background,' Signal: ',signal
  END DO


END PROGRAM

