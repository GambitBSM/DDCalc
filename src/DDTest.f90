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
  CALL NRET_SetNRCoefficient(params_tmp, 'Op1', 0, 1.0d-6) 
  CALL NRET_SetNRCoefficient(params_tmp, 'Op1', 1, -0.5d-6) 
  !WRITE (*,*) 'params_test =', params_tmp(38:)


  WIMP = DDCalc_InitWIMP()
  CALL DDCalc_SetWIMP(WIMP,m=30.0d0,DMtype='NREffectiveTheory',params=params_tmp)
  !CALL DDCalc_SetWIMP(WIMP,m=10d0,DMtype='SDonly',params=[1.0d-6,-0.5d-6])
  !CALL DDCalc_SetWIMP(WIMP,m=30d0,DMtype='SIonly',params=[0.7d-11,-0.3d-11])
  !CALL DDCalc_SetWIMP(WIMP,m=10d0,DMtype='SISD',params=[0.7d-11,-0.3d-11,1.0d-6,-0.5d-6])
  Halo = DDCalc_InitHalo()
  CALL DDCalc_SetHalo(Halo,rho=0.3d0,vrot=220.d0,v0=220.d0)
  Detector = Xenon1T_2017_Init()
  CALL DDCalc_CalcRates(Detector, WIMP, Halo)


  WRITE (*,*) 'InitSuccess =', Detector%InitSuccess
  WRITE (*,*) 'StatisticFlag =', Detector%StatisticFlag
  WRITE (*,*) 'exposure =', Detector%exposure
  WRITE (*,*) 'Nevents =', Detector%Nevents
  WRITE (*,*) 'Backgr =', Detector%Backgr
  WRITE (*,*) 'MuSignal =', Detector%MuSignal




END PROGRAM

