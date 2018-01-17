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
  REAL*8 :: BGlogL,mDM,sigmaSI,mDMmin,mDMmax,sigmaSImin,sigmaSImax,fp,fn,lnp
  INTEGER :: mDMsteps, sigmaSIsteps, mDMi, sigmaSIi

  mDMmin = 1.d0
  mDMmax = 1.d2
  sigmaSImin = 1.d-6
  sigmaSImax = 1.d-1

  mDMsteps = 40
  sigmaSIsteps = 70

  WIMP = DDCalc_InitWIMP()
  CALL DDCalc_SetWIMP(WIMP,m=1d0,DMtype='SIonly',params=[0.d0,0.d0])
  Halo = DDCalc_InitHalo()
  CALL DDCalc_SetHalo(Halo,rho=0.3d0,vrot=220.d0,v0=220.d0)
  Detector = CRESST_II_Init()
  CALL DDCalc_CalcRates(Detector, WIMP, Halo)

  BGlogL = DDCalc_LogLikelihood(Detector)
  WRITE (*,*) 'Background log likelihood =',BGlogL
  WRITE (*,*) 'mDM[GeV] sigmaSI(cm2) logL(signal) -2DeltalogL log(pvalue)'

  DO mDMi = 0,mDMsteps
    mDM = mDMmin * (mDMmax/mDMmin)**(REAL(mDMi)/mDMsteps)
    DO sigmaSIi = 0,sigmaSIsteps
      sigmaSI = sigmaSImin * (sigmaSImax/sigmaSImin)**(REAL(sigmaSIi)/sigmaSIsteps)
      fp = SigmapSItoFp(mDM,sigmaSI)
      fn = SigmanSItoFn(mDM,sigmaSI)
      CALL DDCalc_SetWIMP(WIMP,m=mDM,DMtype='SIonly',params=[fp,fn])
      CALL DDCalc_CalcRates(Detector, WIMP, Halo)
      WRITE (*,*) mDM,sigmaSI,DDCalc_LogLikelihood(Detector),-2*(DDCalc_LogLikelihood(Detector)-BGlogL)
    END DO
  END DO


END PROGRAM

