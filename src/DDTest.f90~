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

  IMPLICIT NONE
  TYPE(DetectorStruct) :: Detector
  TYPE(WIMPStruct) :: WIMP
  TYPE(HaloStruct) :: Halo
  REAL*8 :: BGlogL
  REAL*8 :: params_tmp(45)

  WRITE (*,*) 'This is DDTest.' 

  params_tmp = NRET_CreateCoeffList()
  CALL NRET_SetDMSpin(params_tmp, 0.5d0)
  CALL NRET_SetNRCoefficient(params_tmp, 'Op18', 0, 3.0d0) 
  CALL NRET_SetNRCoefficient(params_tmp, 'Op18', 1, -2.0d0) 
  !WRITE (*,*) 'params_test =', params_tmp(38:)


  WIMP = DDCalc_InitWIMP()
  CALL DDCalc_SetWIMP(WIMP,m=10.0d0,DMtype='NREffectiveTheory',params=params_tmp)
  Halo = DDCalc_InitHalo()
  CALL DDCalc_SetHalo(Halo,rho=0.3d0,vrot=220.d0,v0=220.d0)
  Detector = CRESST_II_Init()
  CALL DDCalc_CalcRates(Detector, WIMP, Halo)



  WRITE (*,*) 'InitSuccess =', Detector%InitSuccess
  WRITE (*,*) 'StatisticFlag =', Detector%StatisticFlag
  WRITE (*,*) 'exposure =', Detector%exposure
  WRITE (*,*) 'Nevents =', Detector%Nevents
  WRITE (*,*) 'Backgr =', Detector%Backgr
  WRITE (*,*) 'MuSignal =', Detector%MuSignal

  BGlogL = DDCalc_LogLikelihood(Detector)
  WRITE (*,*) 'Background log likelihood =',BGlogL



END PROGRAM

