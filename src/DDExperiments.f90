MODULE DDExperiments

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!  DDExperiments
!    Lists of experimental analyses included in DDCalc.  To add any new
!    experimental analysis:
!    1.  Put it in its own individual module source file, put that file
!        in src/analyses, and then USE that module from here.
!    2.  If you want your new analysis to also be available in the
!        commandline version of DDCalc, make sure to add it to the list
!        in AvailableAnalyses below.
!    3.  If you want your new analysis to also be available when calling
!        DDCalc from C/C++, make sure to add it to
!        include/DDExperiments.hpp too.
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

USE XENON100_2012
USE LUX_2013
USE SIMPLE_2014
USE SuperCDMS_2014
USE Darwin_Ar
USE Darwin_Xe
USE LUX_2016
USE PandaX_2016
USE PandaX_2017
USE Xenon1T_2017
USE LUX_2015
USE PICO_2L
USE PICO_60_F
USE PICO_60_I
USE PICO_60
USE PICO_60_2017
USE CRESST_II
USE DummyExp

IMPLICIT NONE

! Change these if you want to pick a different analysis as the default.
INTERFACE DDCalc_InitDefaultDetector
  MODULE PROCEDURE Xenon1T_2017_Init
END INTERFACE

CONTAINS


!-----------------------------------------------------------------------
! List of the analyses available to the command-line interface.  Add to
! this if you want to make new analyses accessible from the command line.
! 
FUNCTION AvailableAnalyses() RESULT(Detector)

  USE DDTypes
  USE DDCommandLine
  USE DDDetectors

  TYPE(DetectorStruct) :: Detector

  ! Add your new analysis here.
  IF (GetLongArg('LUX-2016')) THEN
    Detector = LUX_2016_Init()
  ELSE IF      (GetLongArg('XENON100-2012'))  THEN
    Detector = XENON100_2012_Init()
  ELSE IF (GetLongArg('LUX-2013'))       THEN
    Detector = LUX_2013_Init()
  ELSE IF (GetLongArg('PandaX-2016'))    THEN
    Detector = PandaX_2016_Init()
  ELSE IF (GetLongArg('PandaX-2017'))    THEN
    Detector = PandaX_2017_Init()
  ELSE IF (GetLongArg('Xenon1T-2017'))   THEN
    Detector = Xenon1T_2017_Init()
  ELSE IF (GetLongArg('LUX-2015'))       THEN
    Detector = LUX_2015_Init()
  ELSE IF (GetLongArg('PICO-2L'))        THEN
    Detector = PICO_2L_Init()
  ELSE IF (GetLongArg('PICO-60_F'))      THEN
    Detector = PICO_60_F_Init()
  ELSE IF (GetLongArg('PICO-60_I'))      THEN
    Detector = PICO_60_I_Init()
  ELSE IF (GetLongArg('PICO-60'))        THEN
    Detector = PICO_60_Init()
  ELSE IF (GetLongArg('PICO-60_2017'))   THEN
    Detector = PICO_60_2017_Init()
  ELSE IF (GetLongArg('SuperCDMS-2014')) THEN !(low-energy analysis)
    Detector = SuperCDMS_2014_Init() 
  ELSE IF (GetLongArg('SIMPLE-2014'))    THEN
    Detector = SIMPLE_2014_Init()
  ELSE IF (GetLongArg('DARWIN-Xe'))      THEN
    Detector = DARWIN_Xe_Init()
  ELSE IF (GetLongArg('DARWIN-Ar'))      THEN
    Detector = DARWIN_Ar_Init()
  ELSE IF (GetLongArg('CRESST_II'))      THEN
    Detector = CRESST_II_Init()
  ELSE IF (GetLongArg('DummyExp'))      THEN
    Detector = DummyExp_Init()
  ELSE                                   !Default
    Detector = DDCalc_InitDefaultDetector()
  END IF

END FUNCTION


END MODULE
