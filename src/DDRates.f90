MODULE DDRates

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! DDRates
!    Routines to calculate and view rates.  The CalcRates() routine
!    here must be called after changing the WIMP mass and/or couplings.
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX


USE DDConstants
USE DDTypes
USE DDHalo
USE DDCouplings

IMPLICIT NONE
PRIVATE

! Rates
PUBLIC :: DDCalc_GetRates,DDCalc_CalcRates,C_DDCalc_CalcRates
INTERFACE DDCalc_GetRates
  MODULE PROCEDURE GetRates
END INTERFACE
INTERFACE DDCalc_CalcRates
  MODULE PROCEDURE CalcRates
END INTERFACE

! Event routines (utility)
PUBLIC :: DDCalc_Events,DDCalc_Background
PUBLIC :: DDCalc_Signal
PUBLIC :: C_DDCalc_Events,C_DDCalc_Background
PUBLIC :: C_DDCalc_Signal
INTERFACE DDCalc_Events
  MODULE PROCEDURE GetEvents
END INTERFACE
INTERFACE DDCalc_Background
  MODULE PROCEDURE GetBackground
END INTERFACE
INTERFACE DDCalc_Signal
  MODULE PROCEDURE GetSignal
END INTERFACE

CONTAINS


! ----------------------------------------------------------------------
! Get various rate quantities.  These require CalcRates() to have been
! run first.
! 
! Optional input argument:
!   D           The DetectorStruct structure to extract rate
!               quantities from.
! Optional output arguments:
!   Nevents     Number of observed events
!   background  Average expected background events
!   signal      Average expected signal events
!   signal_si   Average expected spin-independent signal events
!                  --> In the new DDCalc setup, this is set to zero as a dummy value.
!   signal_sd   Average expected spin-dependent signal events
!                  --> In the new DDCalc setup, this is set to zero as a dummy value.
!   rate        Signal rate [cpd/kg]
! Optional sub-interval/bin arguments.  Arrays will be allocated to
! size [1:Nbins]:
!   Nbins       Number of bins/intervals
!   binsignal   Allocatable array to be filled with average expected
!               signal in each bin.  Allocated to size [1:Nbins].
!               The sum of all bins is equal to 'signal'.
!   binrate      Allocatable array to be filled with recoil rate in
!               each bin.  Allocated to size [1:Nbins].
!               The sum of all bins is equal to 'rate'.
!   
SUBROUTINE GetRates(D,Nevents,background,signal,signal_si,signal_sd,    &
                    rate,Nbins,binsignal,binrate)
  IMPLICIT NONE
  TYPE(DetectorStruct), INTENT(IN) :: D
  INTEGER, INTENT(OUT), OPTIONAL :: Nevents,Nbins
  REAL*8, INTENT(OUT), OPTIONAL :: background,signal,signal_si,signal_sd,&
          rate
  REAL*8, ALLOCATABLE, INTENT(OUT), OPTIONAL :: binsignal(:),binrate(:)
  
  IF ( .NOT. D%InitSuccess ) THEN
    WRITE(*,*) 'ERROR: Cannot get information from a detector that has not been correctly initialized.'
    STOP      
  END IF

  ! Observed events and expected background events
  IF (PRESENT(Nevents))    Nevents    = D%Nevents(0)
  IF (PRESENT(background)) background = D%Backgr(0)
  
  ! Signal events
  IF (PRESENT(signal))    signal    = D%MuSignal(0)
  IF (PRESENT(signal_si)) signal_si = 0d0 ! FIXME: delete this dummy assignment
  IF (PRESENT(signal_sd)) signal_sd = 0d0 ! FIXME: delete this dummy assignment
  
  ! Signal rates
  IF (PRESENT(rate))    rate    = D%R(0)
  
  ! Bins
  IF (PRESENT(Nbins)) Nbins = D%Nbins
  
  ! Signal events by bin
  IF (PRESENT(binsignal)) THEN
    ALLOCATE(binsignal(D%Nbins))
    binsignal = D%MuSignal(1:D%Nbins)
  END IF
  
  ! Signal rates by bin
  IF (PRESENT(binrate)) THEN
    ALLOCATE(binrate(D%Nbins))
    binrate = D%MuSignal(1:D%Nbins)
  END IF
  
END SUBROUTINE




! -------------------------------------------------------------------------------
! Auxiliary functions for calculating dRdE for different interaction types
!


! standard SI scattering
FUNCTION dRdE_SI(m, rho, eta, fiso, fp, fn, Wsi1, Wsi0, WsiM1)   &
    RESULT (rate)
  IMPLICIT NONE

  REAL*8 :: rate
  REAL*8, INTENT(IN) :: m, rho, eta, fiso, fp, fn
  REAL*8, INTENT(IN) :: Wsi1, Wsi0, WsiM1
  REAL*8, PARAMETER :: TO_CPD_KG_KEV = 1.695e14 !   s / (cm^3 km GeV^4)  -->  cpd/kg/keV

    ! In the DDCalc/DarkBit convention, fp and fn are the coefficients
    ! of the DM DM N N term in the Lagrangian, assuming a Majorana DM
    ! particle. In particular, one has sigma_p = 4 mu^2 fp^2 / pi.  
  ! dRdE = fiso*rho*eta(vmin)/(2*mDM) * 
  !   ( 4 fp^2 Wsi(+1) + 4 fp fn Wsi(0) + 4 fn^2 Wsi(-1) )
  ! The weighted structure functions Wsi are related to the usual Helm form factor F(E) via
  !   Wsi(+1,E) = (1/pi) Z^2 F^2(E)        ! SI proton
  !   Wsi( 0,E) = (1/pi) 2*Z*(A-Z) F^2(E)  ! SI crossterm
  !   Wsi(-1,E) = (1/pi) (A-Z)^2 F^2(E)    ! SI neutron
  rate = TO_CPD_KG_KEV * fiso * rho * eta/(2*m) * (    & 
                4*fp**2 * Wsi1     &
              + 4*fp*fn * Wsi0     &
              + 4*fn**2 * WsiM1)
  
END FUNCTION


! standard SD scattering
FUNCTION dRdE_SD(m, rho, eta, fiso, ap, an, Wsd1, Wsd0, WsdM1)   &
    RESULT (rate)
  IMPLICIT NONE

  REAL*8 :: rate
  REAL*8, INTENT(IN) :: m, rho, eta, fiso, ap, an
  REAL*8, INTENT(IN) :: Wsd1, Wsd0, WsdM1
  REAL*8, PARAMETER :: TO_CPD_KG_KEV = 1.695e14 !   s / (cm^3 km GeV^4)  -->  cpd/kg/keV

   ! This follows the Jungman/Kamionkowski convention of defining ap and an.
  ! dRdE = 4*fiso*rho*eta(vmin)*Gfermi**2 / mDM * 
  !   ( ap^2*WSd(+1) + ap*an*WSd(0) * an^2*Wsd(-1) )
  ! The weighted structure functions Wsd are related to the Spp,Spn,Snn structure functions via
  !   Wsd(+1,E) = 4/(2J+1) Spp(E)          ! SD proton
  !   Wsd( 0,E) = 4/(2J+1) Spn(E)          ! SD crossterm
  !   Wsd(-1,E) = 4/(2J+1) Snn(E)          ! SD neutron
  rate = TO_CPD_KG_KEV * 4.0d0 * fiso * rho * eta *   &
                FERMI_COUPLING_CONSTANT**2/m * (     &
                ap**2 * Wsd1 + ap*an * Wsd0 +        &
                an**2 * WsdM1 )

END FUNCTION


!
! End of auxiliary functions 
! -----------------------------------------------------------------------------



! ----------------------------------------------------------------------
! Calculates various rate quantities using the passed WIMP and halo.
! This is the main routine intended for putting together the WIMP
! rate calculations and must be called any time the WIMP mass
! and/or couplings are modified.
! 
SUBROUTINE CalcRates(D, WIMP, Halo)
  IMPLICIT NONE
  TYPE(DetectorStruct), INTENT(INOUT) :: D
  TYPE(WIMPStruct), INTENT(IN) :: WIMP
  TYPE(HaloStruct), INTENT(IN) :: Halo
  INTEGER :: Kiso, KE, Keff, Neff

  IF ( .NOT. D%InitSuccess ) THEN
    WRITE(*,*) 'ERROR: Cannot calculate rates for a detector that has not been correctly initialized.'
    STOP      
  END IF

  ! Update mean inverse speed, set dRdEiso to zero
  D%vmin = EToVmin(D%NE,D%E,WIMP%m,D%Niso,D%Miso)
  D%eta  = MeanInverseSpeed(D%NE,D%Niso,D%vmin,Halo)
  D%dRdEiso = 0.0d0



  ! ...................... CALCULATE DIFFERENTIAL RATES ............................
  IF (WIMP%DMtype .EQ. 'SIonly') THEN
    ! fp = WIMP%params(1), fn = WIMP%params(2)
    ! In the DDCalc/DarkBit convention, fp and fn are the coefficients
    ! of the DM DM N N term in the Lagrangian, assuming a Majorana DM
    ! particle. In particular, one has sigma_p = 4 mu^2 fp^2 / pi.  
     DO KE = 1,D%NE
       DO Kiso = 1,D%Niso
         D%dRdEiso(KE,Kiso) = &
             dRdE_SI(WIMP%m, Halo%rho, D%eta(KE,Kiso), &
             D%fiso(Kiso), WIMP%params(1), WIMP%params(2), &
             D%Wsi(+1,KE,Kiso), D%Wsi(0,KE,Kiso), D%Wsi(-1,KE,Kiso))
       END DO
     END DO

  ELSE IF (WIMP%DMtype .EQ. 'SDonly') THEN
   ! ap = WIMP%params(1), an = WIMP%params(2)
   ! This follows the Jungman/Kamionkowski convention of defining ap and an.
     DO KE = 1,D%NE
       DO Kiso = 1,D%Niso
         D%dRdEiso(KE,Kiso) = &
             dRdE_SD(WIMP%m, Halo%rho, D%eta(KE,Kiso), &
             D%fiso(Kiso), WIMP%params(1), WIMP%params(2), &
             D%Wsd(+1,KE,Kiso), D%Wsd(0,KE,Kiso), D%Wsd(-1,KE,Kiso))
       END DO
     END DO


  ELSE IF (WIMP%DMtype .EQ. 'SISD') THEN
   ! combination of SIonly and SD only, see above.
   ! fp = WIMP%params(1), fn = WIMP%params(2), ap = WIMP%params(3), an = WIMP%params(4)
     DO KE = 1,D%NE
       DO Kiso = 1,D%Niso
         D%dRdEiso(KE,Kiso) = &
           dRdE_SI(WIMP%m, Halo%rho, D%eta(KE,Kiso), &
             D%fiso(Kiso), WIMP%params(1), WIMP%params(2), &
             D%Wsi(+1,KE,Kiso), D%Wsi(0,KE,Kiso), D%Wsi(-1,KE,Kiso)) + &
           dRdE_SD(WIMP%m, Halo%rho, D%eta(KE,Kiso), &
             D%fiso(Kiso), WIMP%params(3), WIMP%params(4), &
             D%Wsd(+1,KE,Kiso), D%Wsd(0,KE,Kiso), D%Wsd(-1,KE,Kiso))
       END DO
     END DO

  ELSE IF (WIMP%DMtype .EQ. 'HiggsPortal') THEN
   ! fsN * (DM DM N N) + apN * (i DM G5 GM N N), assuming a Majorana DM particle
   ! The expression for the rate originating from apN is the same as the standard SI interaction,
   !   aside from an additional factor q^2/(4*mDM^2).
   ! fsp = WIMP%params(1), fsn =  WIMP%params(2), app = WIMP%params(3), apn = WIMP%params(4)

     DO KE = 1,D%NE
       DO Kiso = 1,D%Niso
         D%dRdEiso(KE,Kiso) = &
            dRdE_SI(WIMP%m, Halo%rho, D%eta(KE,Kiso), & ! scalar-scalar term
             D%fiso(Kiso), WIMP%params(1), WIMP%params(2), &
             D%Wsi(+1,KE,Kiso), D%Wsi(0,KE,Kiso), D%Wsi(-1,KE,Kiso)) + &
            dRdE_SI(WIMP%m, Halo%rho, D%eta(KE,Kiso), & ! psuedoscalar-scalar term
             D%fiso(Kiso), WIMP%params(3), WIMP%params(4), &
             D%Wsi(+1,KE,Kiso), D%Wsi(0,KE,Kiso), D%Wsi(-1,KE,Kiso)) * &
             (2*D%Miso(Kiso)*D%E(KE)*1d-6)/(4d0 * WIMP%m**2)
       END DO
     END DO

  ELSE
     WRITE (*,*) '... this should produce some error ...'
  END IF


  ! ...................... CALCULATE INTEGRATED RATE ............................
  ! Number of intervals/bins to do calculations for.
  ! If intervals=.FALSE. then efficiency index is over [0:0]
  ! (total only).
  IF (D%intervals) THEN
    Neff = D%Nbins
  ELSE
    Neff = 0
  END IF

  ! Integrate (efficiency-weighted) to find total rates.
  ! Uses a simple trapezoidal integration.
  D%R = 0d0
  ! Cycle over E bins and efficiency curves.
  DO KE = 1,D%NE-1
    DO Kiso = 1,D%Niso
      DO Keff = 0,Neff
        D%R(Keff) = D%R(Keff)                                 &
            + 0.5d0 * (D%E(KE+1) - D%E(KE))                   &
              * (D%eff(Kiso,KE,Keff)*D%dRdEiso(KE,Kiso)           &
                 + D%eff(Kiso,KE+1,Keff)*D%dRdEiso(KE+1,Kiso))
      END DO
    END DO
  END DO

  ! calculate number of events for all Keff
  D%MuSignal(:) = D%exposure * D%R

END SUBROUTINE




!-----------------------------------------------------------------------
! C/C++ wrapper for DDCalc_CalcRates.
!
SUBROUTINE C_DDCalc_CalcRates(DetectorIndex,WIMPIndex,HaloIndex) &
 BIND(C,NAME='C_DDRates_ddcalc_calcrates')
  USE ISO_C_BINDING, only: C_INT
  IMPLICIT NONE
  INTEGER(KIND=C_INT), INTENT(IN) :: DetectorIndex,WIMPIndex,HaloIndex
  IF (.NOT. ASSOCIATED(Detectors(DetectorIndex)%p)) stop 'Invalid detector index given to C_DDCalc_CalcRates'
  IF (.NOT. ASSOCIATED(WIMPs(WIMPIndex)%p)) stop 'Invalid WIMP index given to C_DDCalc_CalcRates'
  IF (.NOT. ASSOCIATED(Halos(HaloIndex)%p)) stop 'Invalid halo index given to C_DDCalc_CalcRates'
  CALL CalcRates(Detectors(DetectorIndex)%p,WIMPs(WIMPIndex)%p,Halos(HaloIndex)%p)
END SUBROUTINE


! ----------------------------------------------------------------------
! Returns the observed number of events.
! 
! Required input argument:
!   D           A DetectorStruct containing detector info.
! 
FUNCTION GetEvents(D) RESULT(N)
  IMPLICIT NONE
  INTEGER :: N
  TYPE(DetectorStruct), INTENT(IN) :: D
  CALL GetRates(D,Nevents=N)
END FUNCTION

!-----------------------------------------------------------------------
! C/C++ wrapper for DDCalc_Events
!
INTEGER(KIND=C_INT) FUNCTION C_DDCalc_Events(DetectorIndex) &
 BIND(C,NAME='C_DDRates_ddcalc_events')
  USE ISO_C_BINDING, only: C_INT
  IMPLICIT NONE
  INTEGER(KIND=C_INT), INTENT(IN) :: DetectorIndex
  IF (.NOT. ASSOCIATED(Detectors(DetectorIndex)%p)) stop 'Invalid detector index given to C_DDCalc_Events'
  C_DDCalc_Events = GetEvents(Detectors(DetectorIndex)%p)
END FUNCTION


! ----------------------------------------------------------------------
! Returns the average expected number of background events.
! 
! Required input argument:
!   D           A DetectorStruct containing detector info.
! 
FUNCTION GetBackground(D) RESULT(b)
  IMPLICIT NONE
  REAL*8 :: b
  TYPE(DetectorStruct), INTENT(IN) :: D
  CALL GetRates(D,background=b)
END FUNCTION

!-----------------------------------------------------------------------
! C/C++ wrapper for DDCalc_Background
!
REAL(KIND=C_DOUBLE) FUNCTION C_DDCalc_Background(DetectorIndex) &
 BIND(C,NAME='C_DDRates_ddcalc_background')
  USE ISO_C_BINDING, only: C_DOUBLE, C_INT
  IMPLICIT NONE
  INTEGER(KIND=C_INT), INTENT(IN) :: DetectorIndex
  IF (.NOT. ASSOCIATED(Detectors(DetectorIndex)%p)) stop 'Invalid detector index given to C_DDCalc_Background'
  C_DDCalc_Background = REAL(GetBackground(Detectors(DetectorIndex)%p),KIND=C_DOUBLE)
END FUNCTION


! ----------------------------------------------------------------------
! Returns the average expected number of signal events for the
! current WIMP.
! 
! Required input argument:
!   D           A DetectorStruct containing rates (CalcRates(D) must
!               have already been called).
! 
FUNCTION GetSignal(D) RESULT(s)
  IMPLICIT NONE
  REAL*8 :: s
  TYPE(DetectorStruct), INTENT(IN) :: D
  CALL GetRates(D,signal=s)
END FUNCTION

!-----------------------------------------------------------------------
! C/C++ wrapper for DDCalc_Signal
!
REAL(KIND=C_DOUBLE) FUNCTION C_DDCalc_Signal(DetectorIndex) &
 BIND(C,NAME='C_DDRates_ddcalc_signal')
  USE ISO_C_BINDING, only: C_DOUBLE, C_INT
  IMPLICIT NONE
  INTEGER(KIND=C_INT), INTENT(IN) :: DetectorIndex
  IF (.NOT. ASSOCIATED(Detectors(DetectorIndex)%p)) stop 'Invalid detector index given to C_DDCalc_Signal'
  C_DDCalc_Signal = REAL(GetSignal(Detectors(DetectorIndex)%p),KIND=C_DOUBLE)
END FUNCTION

END MODULE
