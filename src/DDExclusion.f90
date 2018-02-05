PROGRAM DDLikelihood

  USE DDCALC
  USE DDExperiments
  USE DDConstants
  USE DDCouplings
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
  REAL*8 :: BGlogL,mDM,sigmaSItest,limit,logLlimit,mDMmin,mDMmax,fp,fn,lnp,Gp
  INTEGER :: mDMsteps, sigmaSIsteps, mDMi, sigmaSIi

  mDMmin = 1.d0
  mDMmax = 1.d4

  sigmaSItest = 1.d-8
!  sigmaSItest = 1.d-6

  mDMsteps = 80

  WIMP = DDCalc_InitWIMP()
  CALL DDCalc_SetWIMP(WIMP,m=1d0,DMtype='SIonly',params=[0.d0,0.d0])
  Halo = DDCalc_InitHalo()
  CALL DDCalc_SetHalo(Halo,rho=0.3d0,vrot=220.d0,v0=220.d0)
  Detector = DarkSide_Init()
  CALL DDCalc_CalcRates(Detector, WIMP, Halo)

  BGlogL = DDCalc_LogLikelihood(Detector)
  limit = 1.68
  logLlimit = BGlogL - (limit / 2.d0)

  DO mDMi = 0,mDMsteps
    mDM = mDMmin * (mDMmax/mDMmin)**(REAL(mDMi)/mDMsteps)
    fp = SigmapSItoFp(mDM,sigmaSItest)
    fn = SigmanSItoFn(mDM,sigmaSItest)
!    fp = SigmapSDtoAp(mDM,sigmaSItest)
!    fn = SigmanSDtoAn(mDM,sigmaSItest)
    CALL DDCalc_SetWIMP(WIMP,m=mDM,DMtype='SIonly',params=[fp,fn])
!    CALL DDCalc_SetWIMP(WIMP,m=mDM,DMtype='SDonly',params=[0d0,fn])
    CALL DDCalc_CalcRates(Detector, WIMP, Halo)
    WRITE (*,*) mDM,DDCalc_ScaleToPValue(Detector,logLlimit)*sigmaSItest
  END DO

END PROGRAM

