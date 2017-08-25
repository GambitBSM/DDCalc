



PROGRAM DDTest

  USE DDCALC
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
  USE DDNuclear

  IMPLICIT NONE
  TYPE(DetectorStruct) :: Detector
  TYPE(WIMPStruct) :: WIMP
  TYPE(HaloStruct) :: Halo

  WRITE (*,*) 'This is DDTest.'
 
  WIMP = DDCalc_InitWIMP()
  !CALL DDCalc_SetWIMP(WIMP,m=20d0,DMtype='SIonly',params=[1d-9,0d0])
  !CALL DDCalc_SetWIMP(WIMP,m=20d0,DMtype='SDonly',params=[4d-2,0d0])
  CALL DDCalc_SetWIMP(WIMP,m=20d0,DMtype='SISD',params=[1d-8, -0.5d-8, 1d-1, -0.3d-1])
  !CALL DDCalc_SetWIMP_Higgsportal(WIMP, 20d0, 1d-6, -0.5d-8, 5d-6, -0.6d-6)
  Halo = DDCalc_InitHalo()

  Detector = PICO_60_F_Init(.true.)
  CALL DDCalc_CalcRates(Detector, WIMP, Halo)
  WRITE (*,*) Detector%MuSignal  

  Detector = PICO_60_I_Init(.true.)
  CALL DDCalc_CalcRates(Detector, WIMP, Halo)
  WRITE (*,*) Detector%MuSignal\

  Detector = PICO_60_Init(.true.)
  CALL DDCalc_CalcRates(Detector, WIMP, Halo)
  WRITE (*,*) Detector%MuSignal\

  Detector = Xenon1T_2017_Init(.true.)
  CALL DDCalc_CalcRates(Detector, WIMP, Halo)
  WRITE (*,*) Detector%MuSignal\

  Detector = PICO_60_2017_Init(.true.)
  CALL DDCalc_CalcRates(Detector, WIMP, Halo)
  WRITE (*,*) Detector%MuSignal\

END PROGRAM
