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
  TYPE(DetectorStruct) :: Detector1, Detector2, Detector3
  TYPE(WIMPStruct) :: WIMP
  TYPE(HaloStruct) :: Halo
  REAL*8 :: BGlogL1,BGlogL2,BGlogL3,mDM,sigmaSI,mDMmin,mDMmax,sigmaSImin,sigmaSImax,fp,fn,lnp
  INTEGER :: mDMsteps, sigmaSIsteps, mDMi, sigmaSIi

  mDMmin = 7.d-1
  mDMmax = 1.3d1
  sigmaSImin = 1.d-6
  sigmaSImax = 1.d2

  mDMsteps = 80
  sigmaSIsteps = 80

  WIMP = DDCalc_InitWIMP()
  CALL DDCalc_SetWIMP(WIMP,m=1d0,DMtype='SIonly',params=[0.d0,0.d0])
  Halo = DDCalc_InitHalo()
  CALL DDCalc_SetHalo(Halo,rho=0.3d0,vrot=220.d0,v0=220.d0)
  Detector1 = CRESST_II_Init()
  Detector2 = CDMSlite_Init()
  Detector3 = LUX_2016_Init()
  CALL DDCalc_CalcRates(Detector1, WIMP, Halo)
  CALL DDCalc_CalcRates(Detector2, WIMP, Halo)
  CALL DDCalc_CalcRates(Detector3, WIMP, Halo)

  BGlogL1 = DDCalc_LogLikelihood(Detector1)
  BGlogL2 = DDCalc_LogLikelihood(Detector2)
  BGlogL3 = DDCalc_LogLikelihood(Detector3)

!  WRITE (*,*) 'Background log likelihood =',BGlogL
!  WRITE (*,*) 'mDM[GeV] sigmaSI(cm2) logL(signal) -2DeltalogL log(pvalue)'

  DO mDMi = 0,mDMsteps
    mDM = mDMmin * (mDMmax/mDMmin)**(REAL(mDMi)/mDMsteps)
    DO sigmaSIi = 0,sigmaSIsteps
      sigmaSI = sigmaSImin * (sigmaSImax/sigmaSImin)**(REAL(sigmaSIi)/sigmaSIsteps)
      fp = SigmapSItoFp(mDM,sigmaSI)
      fn = SigmanSItoFn(mDM,sigmaSI)
      CALL DDCalc_SetWIMP(WIMP,m=mDM,DMtype='SIonly',params=[fp,fn])
      CALL DDCalc_CalcRates(Detector1, WIMP, Halo)
      CALL DDCalc_CalcRates(Detector2, WIMP, Halo)
      CALL DDCalc_CalcRates(Detector3, WIMP, Halo)
      WRITE (*,*) mDM,sigmaSI,-2*(DDCalc_LogLikelihood(Detector1)+DDCalc_LogLikelihood(Detector2) &
                                  +DDCalc_LogLikelihood(Detector3)-BGlogL1-BGlogL2-BGlogL3)
    END DO
  END DO


END PROGRAM

