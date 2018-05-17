PROGRAM DDCalc_test

  USE DDCALC
  USE DDExperiments

  IMPLICIT NONE
  TYPE(DetectorStruct) :: Detector
  TYPE(WIMPStruct) :: WIMP
  TYPE(HaloStruct) :: Halo

! Initialize WIMP, DM Halo and detector objects with default values
  WIMP = DDCalc_InitWIMP()
  Halo = DDCalc_InitHalo()
  Detector = DDCalc_InitDetector()

  WRITE (*,*) 'DDCalc successfully initialized'

END PROGRAM

