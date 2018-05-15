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
USE DDNREffectiveTheory

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
PUBLIC :: DDCalc_Signal,DDCalc_Bins
PUBLIC :: DDCalc_BinEvents,DDCalc_BinBackground
PUBLIC :: DDCalc_BinSignal

PUBLIC :: C_DDCalc_Events,C_DDCalc_Background
PUBLIC :: C_DDCalc_Signal, C_DDCalc_Bins
PUBLIC :: C_DDCalc_BinEvents,C_DDCalc_BinBackground
PUBLIC :: C_DDCalc_BinSignal

INTERFACE DDCalc_Events
  MODULE PROCEDURE GetEvents
END INTERFACE
INTERFACE DDCalc_Background
  MODULE PROCEDURE GetBackground
END INTERFACE
INTERFACE DDCalc_Signal
  MODULE PROCEDURE GetSignal
END INTERFACE
INTERFACE DDCalc_Bins
  MODULE PROCEDURE GetBins
END INTERFACE
INTERFACE DDCalc_BinEvents
  MODULE PROCEDURE GetBinEvents
END INTERFACE
INTERFACE DDCalc_BinBackground
  MODULE PROCEDURE GetBinBackground
END INTERFACE
INTERFACE DDCalc_BinSignal
  MODULE PROCEDURE GetBinSignal
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
!   events      Number of observed events
!   background  Average expected background events
!   signal      Average expected signal events
!   rate        Signal rate [cpd/kg]
!   Nbins       Number of bins
! Optional bin arguments.  Arrays will be allocated to size [0:Nbins]:
!   binevents   Allocatable array to be filled with observed events
!               in each bin.  Allocated to size [0:Nbins].
!               The entry with index 0 contains the sum of all bins, which is equal to 'events'.
!   binbackground  Allocatable array to be filled with average expected
!               background in each bin.  Allocated to size [1:Nbins].
!               The entry with index 0 contains the sum of all bins, which is equal to 'background'.
!   binsignal   Allocatable array to be filled with average expected
!               signal in each bin.  Allocated to size [1:Nbins].
!               The entry with index 0 contains the sum of all bins, which is equal to 'signal'.
!   binrate      Allocatable array to be filled with recoil rate in
!               each bin.  Allocated to size [1:Nbins].
!               The entry with index 0 contains the sum of all bins, which is equal to 'rate'.
!   
SUBROUTINE GetRates(D,events,background,signal,   &
                    rate,Nbins,binevents,binbackground, &
                    binsignal,binrate)
  IMPLICIT NONE
  TYPE(DetectorStruct), INTENT(IN) :: D
  INTEGER, INTENT(OUT), OPTIONAL :: events,Nbins
  REAL*8, INTENT(OUT), OPTIONAL :: background,signal,rate

  INTEGER, ALLOCATABLE, INTENT(OUT), OPTIONAL :: binevents(:)
  REAL*8, ALLOCATABLE, INTENT(OUT), OPTIONAL :: binbackground(:),binsignal(:),binrate(:)
  
  IF ( .NOT. D%InitSuccess ) THEN
    WRITE(*,*) 'ERROR: Cannot get information from a detector that has not been correctly initialized.'
    STOP      
  END IF

  ! Observed events and expected background events
  IF (PRESENT(events))    events    = D%Nevents(0)
  IF (PRESENT(background)) background = D%Backgr(0)
  
  ! Signal events
  IF (PRESENT(signal))    signal    = D%MuSignal(0)
  
  ! Signal rates
  IF (PRESENT(rate))    rate    = D%R(0)
  
  ! Bins
  IF (PRESENT(Nbins)) Nbins = D%Nbins
  
  ! Observed events and expected background events by bin
  IF (PRESENT(binevents))  THEN
    ALLOCATE(binevents(0:D%Nbins))
    binevents = D%Nevents(0:D%Nbins)
  END IF

  IF (PRESENT(binbackground))  THEN
    ALLOCATE(binbackground(0:D%Nbins))
    binbackground = D%Backgr(0:D%Nbins)
  END IF

  ! Signal events by bin
  IF (PRESENT(binsignal)) THEN
    ALLOCATE(binsignal(0:D%Nbins))
    binsignal = D%MuSignal(0:D%Nbins)
  END IF
  
  ! Signal rates by bin
  IF (PRESENT(binrate)) THEN
    ALLOCATE(binrate(0:D%Nbins))
    binrate = D%R(0:D%Nbins)
  END IF
  
END SUBROUTINE


! -------------------------------------------------------------------------------
! Auxiliary functions for calculating dRdE for different interaction types
!


! standard SI scattering
FUNCTION dRdE_SI(m, rho, eta, fiso, fp, fn, WTilde1_00, WTilde1_01, WTilde1_11)   &
    RESULT (rate)
  IMPLICIT NONE

  REAL*8 :: rate
  REAL*8, INTENT(IN) :: m, rho, eta, fiso, fp, fn
  REAL*8, INTENT(IN) :: WTilde1_00, WTilde1_01, WTilde1_11
  REAL*8, PARAMETER :: TO_CPD_KG_KEV = 1.695e14 !   s / (cm^3 km GeV^4)  -->  cpd/kg/keV

    ! In the DDCalc/DarkBit convention, fp and fn are the coefficients
    ! of the DM DM N N term in the Lagrangian, assuming a Majorana DM
    ! particle. In particular, one has sigma_p = 4 mu^2 fp^2 / pi.  
  ! dRdE = fiso*rho*eta(vmin)/(2*pi*mDM) * 
  !   ( 4 fp^2 WTilde1_pp + 4 fp fn WTilde1_pn + 4 fn^2 WTilde1_nn )
  ! The weighted structure functions WTilde1_XX are related to the usual Helm form factor F(E) via
  !   WTilde1_pp = WTilde1_00 + WTilde1_11 + 2*WTilde1_01 
  !   WTilde1_pn = 2*(WTilde1_00 - WTilde1_11) 
  !   WTilde1_nn = WTilde1_00 + WTilde1_11 - 2*WTilde1_01 
  !
  !   WTilde1_00 = 1/4 * A^2 F^2(E)        ! SI proton
  !   WTilde1_01 = 1/4 * (2Z-A)^2 F^2(E)  ! SI crossterm
  !   WTilde1_11 = 1/4 * A(2Z-A) F^2(E)    ! SI neutron
  rate = TO_CPD_KG_KEV * fiso * rho * eta/(2*PI*m) * (    & 
                4*fp**2 * (WTilde1_00 + WTilde1_11 + 2*WTilde1_01)     &
              + 4*fp*fn * 2*(WTilde1_00 - WTilde1_11)     &
              + 4*fn**2 * (WTilde1_00 + WTilde1_11 - 2*WTilde1_01))
  
END FUNCTION


! standard SD scattering
FUNCTION dRdE_SD(m, rho, eta, fiso, ap, an, WTilde56_00, WTilde56_01, WTilde56_11)   &
    RESULT (rate)
  IMPLICIT NONE

  REAL*8 :: rate
  REAL*8, INTENT(IN) :: m, rho, eta, fiso, ap, an
  REAL*8, INTENT(IN) :: WTilde56_00, WTilde56_01, WTilde56_11
  REAL*8, PARAMETER :: TO_CPD_KG_KEV = 1.695e14 !   s / (cm^3 km GeV^4)  -->  cpd/kg/keV

   ! This follows the Jungman/Kamionkowski convention of defining ap and an, assuming a Majorana DM particle.
   ! dRdE = fiso*rho*eta(vmin)/(2*pi*mDM) * 8 GFermi^2 * (
   !          ap^2 (WTilde56_00 + WTilde56_11 + 2*WTilde56_01) +  
   !          ap*an (2*WTilde56_00 - 2*WTilde56_11) + 
   !          an^2 (WTilde56_00 + WTilde56_11 - 2*WTilde56_01) 
   !       )
   ! Here, WTilde56_XY = WTilde5_XY + WTilde6_XY
  rate = TO_CPD_KG_KEV * fiso * rho * eta/(2*PI*m) *   &
                8 * FERMI_COUPLING_CONSTANT**2 * (     &
                ap**2 * (WTilde56_00 + WTilde56_11 + 2*WTilde56_01) + &
                ap*an * (2*WTilde56_00 - 2*WTilde56_11) + &
                an**2 * (WTilde56_00 + WTilde56_11 - 2*WTilde56_01))

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
  INTEGER :: Kiso, KE, Keff, Neff, alpha
  REAL*8 :: S1S2(8)

  IF ( .NOT. D%InitSuccess ) THEN
    WRITE(*,*) 'ERROR: Cannot calculate rates for a detector that has not been correctly initialized.'
    STOP      
  END IF

  ! Update mean inverse speed and mean speed, set dRdEiso to zero
  D%vmin = EToVmin(D%NE,D%E,WIMP%m,D%Niso,D%Miso)
  D%g_vmin  = MeanInverseSpeed(D%NE,D%Niso,D%vmin,Halo)
  D%h_vmin  = MeanSpeed(D%NE,D%Niso,D%vmin,Halo)
  D%dRdEiso = 0.0d0

  ! ...................... CALCULATE DIFFERENTIAL RATES ............................
  IF (WIMP%DMtype .EQ. 'SIonly') THEN
    ! fp = WIMP%params(1), fn = WIMP%params(2)
    ! In the DDCalc/DarkBit convention, fp and fn are the coefficients
    ! of the DM DM N N term in the Lagrangian, assuming a Majorana DM
    ! particle. In particular, one has sigma_p = 4 mu^2 fp^2 / pi.
    ! Notice: WTilde1_00 corresponds to D%WTilde(1,1,KE,Kiso) 
    ! Notice: WTilde1_01 corresponds to D%WTilde(1,2,KE,Kiso) 
    ! Notice: WTilde1_11 corresponds to D%WTilde(1,4,KE,Kiso) 

    ! If the parameter PreferNewFF is set to .TRUE., use the Haxton form factor for SI scattering
    IF (PreferNewFF) THEN

     DO KE = 1,D%NE
       DO Kiso = 1,D%Niso
         D%dRdEiso(KE,Kiso) = &
             dRdE_SI(WIMP%m, Halo%rho, D%g_vmin(KE,Kiso), &
             D%fiso(Kiso), WIMP%params(1), WIMP%params(2), &
             D%WTilde(1,1,KE,Kiso), D%WTilde(1,2,KE,Kiso), D%WTilde(1,4,KE,Kiso))
       END DO
     END DO


   ! If the parameter PreferNewFF is set to .FALSE., use the Helm form factor for SI scattering
   ELSE

     DO KE = 1,D%NE
       DO Kiso = 1,D%Niso
         D%dRdEiso(KE,Kiso) = &
             dRdE_SI(WIMP%m, Halo%rho, D%g_vmin(KE,Kiso), &
             D%fiso(Kiso), WIMP%params(1), WIMP%params(2), &
             D%WTilde(0,1,KE,Kiso), D%WTilde(0,2,KE,Kiso), D%WTilde(0,4,KE,Kiso))
       END DO
     END DO

   END IF


  ELSE IF (WIMP%DMtype .EQ. 'SDonly') THEN
   ! ap = WIMP%params(1), an = WIMP%params(2)
   ! This follows the Jungman/Kamionkowski convention of defining ap and an.

   ! If the parameter PreferNewFF is set to .TRUE., use the Haxton form factor for SD scattering
   IF(PreferNewFF) THEN

     DO KE = 1,D%NE
       DO Kiso = 1,D%Niso
         D%dRdEiso(KE,Kiso) = &
             dRdE_SD(WIMP%m, Halo%rho, D%g_vmin(KE,Kiso), &
             D%fiso(Kiso), WIMP%params(1), WIMP%params(2), &
              D%WTilde(5,1,KE,Kiso)+D%WTilde(6,1,KE,Kiso), &
              D%WTilde(5,2,KE,Kiso)+D%WTilde(6,2,KE,Kiso), &
              D%WTilde(5,4,KE,Kiso)+D%WTilde(6,4,KE,Kiso))
       END DO
     END DO


   ! If the parameter PreferNewFF is set to .FALSE., use the form factor from Klos et. al. for SD scattering
   ELSE

     DO KE = 1,D%NE
       DO Kiso = 1,D%Niso
         D%dRdEiso(KE,Kiso) = &
             dRdE_SD(WIMP%m, Halo%rho, D%g_vmin(KE,Kiso), &
             D%fiso(Kiso), WIMP%params(1), WIMP%params(2), &
              D%WTilde(9,1,KE,Kiso), D%WTilde(9,2,KE,Kiso), &
              D%WTilde(9,4,KE,Kiso))
       END DO
     END DO

   END IF

  ELSE IF (WIMP%DMtype .EQ. 'SISD') THEN
   ! combination of SIonly and SD only, see above.
   ! fp = WIMP%params(1), fn = WIMP%params(2), ap = WIMP%params(3), an = WIMP%params(4)


   ! If the parameter PreferNewFF is set to .TRUE., use the Haxton form factor for SI and SD scattering
   IF(PreferNewFF) THEN

     DO KE = 1,D%NE
       DO Kiso = 1,D%Niso
         D%dRdEiso(KE,Kiso) = &
             dRdE_SI(WIMP%m, Halo%rho, D%g_vmin(KE,Kiso), &
               D%fiso(Kiso), WIMP%params(1), WIMP%params(2), &
               D%WTilde(1,1,KE,Kiso), D%WTilde(1,2,KE,Kiso), D%WTilde(1,4,KE,Kiso)) + & 
             dRdE_SD(WIMP%m, Halo%rho, D%g_vmin(KE,Kiso), &
               D%fiso(Kiso), WIMP%params(3), WIMP%params(4), &
               D%WTilde(5,1,KE,Kiso)+D%WTilde(6,1,KE,Kiso), &
               D%WTilde(5,2,KE,Kiso)+D%WTilde(6,2,KE,Kiso), &
               D%WTilde(5,4,KE,Kiso)+D%WTilde(6,4,KE,Kiso))
       END DO
     END DO


   ! If the parameter PreferNewFF is set to .FALSE., use the Helm form factor for SI scattering,
   !   and the Klos form factor for SD scattering
   ELSE

     DO KE = 1,D%NE
       DO Kiso = 1,D%Niso
         D%dRdEiso(KE,Kiso) = &
             dRdE_SI(WIMP%m, Halo%rho, D%g_vmin(KE,Kiso), &
               D%fiso(Kiso), WIMP%params(1), WIMP%params(2), &
               D%WTilde(0,1,KE,Kiso), D%WTilde(0,2,KE,Kiso), D%WTilde(0,4,KE,Kiso)) + & 
             dRdE_SD(WIMP%m, Halo%rho, D%g_vmin(KE,Kiso), &
               D%fiso(Kiso), WIMP%params(3), WIMP%params(4), &
              D%WTilde(9,1,KE,Kiso), D%WTilde(9,2,KE,Kiso), &
              D%WTilde(9,4,KE,Kiso))
       END DO
     END DO

   END IF


  ELSE IF (WIMP%DMtype .EQ. 'NREffectiveTheory') THEN

   ! WIMP%params has to be a 45-element list, interpreted as coefficients in units GeV^(-2) of the non-relativistic operators
   ! (DM spin, O1_0, O1_1, O1q2_0, O1q2_1, O3_0, O3_1, O4_0, O4_1, O4q2_0, O4q2_1, O5_0, O5_1, O6_0, O6_1, ..., O15_0, O15_1, O17_0, O17_1, O18_0, O18_1, alpha_1, ..., alpha_8)
   ! See DDTypes for more information.
   !
   ! Notice that the form factors from Haxton are used for this WIMPType. For using instead the Helm form factor for SI scattering,
   !   or the Klos form factor for SD scattering, use the types SIonly, SDonly or SISD, together with the switch PreferNewFF = .FALSE.
     
     IF (WIMP%params(1).LE.0) THEN
        WRITE (*,*) 'Error in using WIMP type NREffectiveTheory: the dark matter spin is not set correctly.'
        STOP
     END IF 

     DO KE = 1,D%NE
       DO Kiso = 1,D%Niso 
         D%dRdEiso(KE,Kiso) = 0.0
         DO alpha = 1,8
            IF (abs(WIMP%params(37 + alpha))>0) THEN ! only calculate and add those terms in alpha for which the corresponding param entry is non-zero
              S1S2 = NRET_SFunctions(D, WIMP%m, WIMP%params, alpha, KE, Kiso)
              D%dRdEiso(KE,Kiso) = D%dRdEiso(KE,Kiso) + &
                  1.6961e14 * dot_product(S1S2(:4), D%WTilde(alpha,:,KE,Kiso)) * D%g_vmin(KE,Kiso)
              D%dRdEiso(KE,Kiso) = D%dRdEiso(KE,Kiso) + &
                  1.8871e3 * dot_product(S1S2(5:), D%WTilde(alpha,:,KE,Kiso)) * D%h_vmin(KE,Kiso)
            END IF
         END DO
         D%dRdEiso(KE,Kiso) = D%fiso(Kiso)*Halo%rho*D%dRdEiso(KE,Kiso)/(2*PI*WIMP%m)
       END DO
     END DO

  ELSE
     WRITE (*,*) 'Error in CalcRates: invalid WIMP%DMType'
     STOP
  END IF


  ! ...................... CALCULATE INTEGRATED RATE ............................
  ! Number of bins to do calculations for.
  IF (D%StatisticFlag .GT. 0) THEN
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
  CALL GetRates(D,events=N)
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
  C_DDCalc_Events = INT(GetEvents(Detectors(DetectorIndex)%p),KIND=C_INT)
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


! ----------------------------------------------------------------------
! Returns the number of bins.
! 
! Required input argument:
!   D           A DetectorStruct containing detector info.
! 
FUNCTION GetBins(D) RESULT(N)
  IMPLICIT NONE
  INTEGER :: N
  TYPE(DetectorStruct), INTENT(IN) :: D
  CALL GetRates(D,Nbins=N)
END FUNCTION

!-----------------------------------------------------------------------
! C/C++ wrapper for DDCalc_Events
!
INTEGER(KIND=C_INT) FUNCTION C_DDCalc_Bins(DetectorIndex) &
 BIND(C,NAME='C_DDRates_ddcalc_bins')
  USE ISO_C_BINDING, only: C_INT
  IMPLICIT NONE
  INTEGER(KIND=C_INT), INTENT(IN) :: DetectorIndex
  IF (.NOT. ASSOCIATED(Detectors(DetectorIndex)%p)) stop 'Invalid detector index given to C_DDCalc_Events'
  C_DDCalc_Bins = INT(GetBins(Detectors(DetectorIndex)%p),KIND=C_INT)
END FUNCTION



! ----------------------------------------------------------------------
! Returns the observed number of events in a given bin.
! 
! Required input argument:
!   D           A DetectorStruct containing detector info.
!   ibin        Index of the bin to be returned.
! 
FUNCTION GetBinEvents(D, ibin) RESULT(N)
  IMPLICIT NONE
  INTEGER :: N
  TYPE(DetectorStruct), INTENT(IN) :: D
  INTEGER, INTENT(IN) :: ibin
  INTEGER :: Nbins
  INTEGER, ALLOCATABLE :: bins(:)
  Nbins = GetBins(D)
  CALL GetRates(D,binevents=bins)
  IF ((ibin >= 0) .AND. (ibin <= Nbins)) THEN
    N = bins(ibin)
  ELSE
    stop 'Bin index out of range!'
  END IF
END FUNCTION

!-----------------------------------------------------------------------
! C/C++ wrapper for DDCalc_BinEvents
!
INTEGER(KIND=C_INT) FUNCTION C_DDCalc_BinEvents(DetectorIndex, BinIndex) &
 BIND(C,NAME='C_DDRates_ddcalc_binevents')
  USE ISO_C_BINDING, only: C_INT
  IMPLICIT NONE
  INTEGER(KIND=C_INT), INTENT(IN) :: DetectorIndex
  INTEGER(KIND=C_INT), INTENT(IN) :: BinIndex
  IF (.NOT. ASSOCIATED(Detectors(DetectorIndex)%p)) stop 'Invalid detector index given to C_DDCalc_Events'
  C_DDCalc_BinEvents = INT(GetBinEvents(Detectors(DetectorIndex)%p,BinIndex),KIND=C_INT)
END FUNCTION

! ----------------------------------------------------------------------
! Returns the expected background in a given bin.
! 
! Required input argument:
!   D           A DetectorStruct containing detector info.
!   ibin        Index of the bin to be returned.
! 
FUNCTION GetBinBackground(D, ibin) RESULT(N)
  IMPLICIT NONE
  REAL*8 :: N
  TYPE(DetectorStruct), INTENT(IN) :: D
  INTEGER, INTENT(IN) :: ibin
  INTEGER :: Nbins
  REAL*8, ALLOCATABLE :: bins(:)
  Nbins = GetBins(D)
  CALL GetRates(D,binbackground=bins)
  IF ((ibin >= 0) .AND. (ibin <= Nbins)) THEN
    N = bins(ibin)
  ELSE
    stop 'Bin index out of range!'
  END IF
END FUNCTION

!-----------------------------------------------------------------------
! C/C++ wrapper for DDCalc_BinBackground
!
REAL(KIND=C_DOUBLE) FUNCTION C_DDCalc_BinBackground(DetectorIndex, BinIndex) &
 BIND(C,NAME='C_DDRates_ddcalc_binbackground')
  USE ISO_C_BINDING, only: C_DOUBLE, C_INT
  IMPLICIT NONE
  INTEGER(KIND=C_INT), INTENT(IN) :: DetectorIndex
  INTEGER(KIND=C_INT), INTENT(IN) :: BinIndex
  IF (.NOT. ASSOCIATED(Detectors(DetectorIndex)%p)) stop 'Invalid detector index given to C_DDCalc_Events'
  C_DDCalc_BinBackground = REAL(GetBinBackground(Detectors(DetectorIndex)%p, BinIndex),KIND=C_DOUBLE)
END FUNCTION

! ----------------------------------------------------------------------
! Returns the expected signal in a given bin.
! 
! Required input argument:
!   D           A DetectorStruct containing detector info.
!   ibin        Index of the bin to be returned.
! 
FUNCTION GetBinSignal(D, ibin) RESULT(N)
  IMPLICIT NONE
  REAL*8 :: N
  TYPE(DetectorStruct), INTENT(IN) :: D
  INTEGER, INTENT(IN) :: ibin
  INTEGER :: Nbins
  REAL*8, ALLOCATABLE :: bins(:)
  Nbins = GetBins(D)
  CALL GetRates(D,binsignal=bins)
  IF ((ibin >= 0) .AND. (ibin <= Nbins)) THEN
    N = bins(ibin)
  ELSE
    stop 'Bin index out of range!'
  END IF
END FUNCTION

!-----------------------------------------------------------------------
! C/C++ wrapper for DDCalc_BinSignal
!
REAL(KIND=C_DOUBLE) FUNCTION C_DDCalc_BinSignal(DetectorIndex, BinIndex) &
 BIND(C,NAME='C_DDRates_ddcalc_binsignal')
  USE ISO_C_BINDING, only: C_DOUBLE, C_INT
  IMPLICIT NONE
  INTEGER(KIND=C_INT), INTENT(IN) :: DetectorIndex
  INTEGER(KIND=C_INT), INTENT(IN) :: BinIndex
  IF (.NOT. ASSOCIATED(Detectors(DetectorIndex)%p)) stop 'Invalid detector index given to C_DDCalc_Events'
  C_DDCalc_BinSignal = REAL(GetBinSignal(Detectors(DetectorIndex)%p, BinIndex),KIND=C_DOUBLE)
END FUNCTION

END MODULE
