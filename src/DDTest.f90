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

  IMPLICIT NONE
  TYPE(DetectorStruct) :: Detector
  TYPE(WIMPStruct) :: WIMP
  TYPE(HaloStruct) :: Halo
  REAL*8 :: BGlogL
  WRITE (*,*) 'This is DDTest.' 

  WIMP = DDCalc_InitWIMP()
  !CALL DDCalc_SetWIMP(WIMP,m=20d0,DMtype='SDonly',params=[4d-2,0d0])
  !CALL DDCalc_SetWIMP(WIMP,m=20d0,DMtype='SISD',params=[1d-8, -0.5d-8, 1d-1, -0.3d-1])
  !CALL DDCalc_SetWIMP_Higgsportal(WIMP, 20d0, 1d-8, -0.5d-8, 5d-6, -0.6d-6)
  CALL DDCalc_SetWIMP(WIMP,m=1d0,DMtype='SIonly',params=[0.d0,0.d0])
  Halo = DDCalc_InitHalo()
  CALL DDCalc_SetHalo(Halo,rho=0.3d0,vrot=220.d0,v0=220.d0)
  Detector = CRESST_II_Init()
  CALL DDCalc_CalcRates(Detector, WIMP, Halo)

  WRITE (*,*) 'InitSuccess =', Detector%InitSuccess
  WRITE (*,*) 'StatisticFlag =', Detector%StatisticFlag
  WRITE (*,*) 'exposure =', Detector%exposure
  WRITE (*,*) 'Nevents =', Detector%Nevents
  WRITE (*,*) 'Backgr =', Detector%Backgr

  BGlogL = DDCalc_LogLikelihood(Detector)
  WRITE (*,*) 'Background log likelihood =',BGlogL

  CALL DDCalc_SetWIMP(WIMP,m=1d0,DMtype='SIonly',params=[8.55d-6,8.55d-6])
  CALL DDCalc_CalcRates(Detector, WIMP, Halo)
  !CALL DDCalc_SetDetector(Detector, Nevents_tot = -1)
  WRITE (*,*) 'MuSignal =',Detector%MuSignal
  WRITE (*,*) 'log likelihood =',DDCalc_LogLikelihood(Detector)
  WRITE (*,*) '-2 Delta log L =',-2*(DDCalc_LogLikelihood(Detector)-BGlogL)



END PROGRAM

