MODULE DDRates

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! DDRates
!    Routines to calculate and view rates.  The CalcRates() routine
!    here must be called after changing the WIMP mass and/or couplings.
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

USE DDTypes
USE DDHalo

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
PUBLIC :: DDCalc_Signal,DDCalc_SignalSI,DDCalc_SignalSD
PUBLIC :: C_DDCalc_Events,C_DDCalc_Background
PUBLIC :: C_DDCalc_Signal,C_DDCalc_SignalSI,C_DDCalc_SignalSD
INTERFACE DDCalc_Events
  MODULE PROCEDURE GetEvents
END INTERFACE
INTERFACE DDCalc_Background
  MODULE PROCEDURE GetBackground
END INTERFACE
INTERFACE DDCalc_Signal
  MODULE PROCEDURE GetSignal
END INTERFACE
INTERFACE DDCalc_SignalSI
  MODULE PROCEDURE GetSignalSI
END INTERFACE
INTERFACE DDCalc_SignalSD
  MODULE PROCEDURE GetSignalSD
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
!   signal_sd   Average expected spin-dependent signal events
!   rate        Signal rate [cpd/kg]
!   rate_si     Spin-independent signal rate [cpd/kg]
!   rate_sd     Spin-dependent signal rate [cpd/kg]
! Optional sub-interval/bin arguments.  Arrays will be allocated to
! size [1:Nbins]:
!   Nbins       Number of bins/intervals
!   binsignal   Allocatable array to be filled with average expected
!               signal in each bin.  Allocated to size [1:Nbins].
!               The sum of all bins is equal to 'signal'.
!   binsignal_si,binsignal_sd
!               Same as 'binsignal' but for only spin-independent or
!               spin-dependent events, respectively.
!   binrate      Allocatable array to be filled with recoil rate in
!               each bin.  Allocated to size [1:Nbins].
!               The sum of all bins is equal to 'rate'.
!   binrate_si,binrate_sd
!               Same as 'binrate' but for only spin-independent or
!               spin-dependent events, respectively.
!   
SUBROUTINE GetRates(D,Nevents,background,signal,signal_si,signal_sd,    &
                    rate,rate_si,rate_sd,Nbins,binsignal,binsignal_si,  &
                    binsignal_sd,binrate,binrate_si,binrate_sd)
  IMPLICIT NONE
  TYPE(DetectorStruct), INTENT(IN) :: D
  INTEGER, INTENT(OUT), OPTIONAL :: Nevents,Nbins
  REAL*8, INTENT(OUT), OPTIONAL :: background,signal,signal_si,signal_sd,&
          rate,rate_si,rate_sd
  REAL*8, ALLOCATABLE, INTENT(OUT), OPTIONAL :: binsignal(:),           &
          binsignal_si(:),binsignal_sd(:),binrate(:),binrate_si(:),     &
          binrate_sd(:)
  INTEGER :: Neff
    
  Neff = D%Neff
  
  ! Observed events and expected background events
  IF (PRESENT(Nevents))    Nevents    = D%Nevents
  IF (PRESENT(background)) background = D%MuBackground
  
  ! Signal events
  IF (PRESENT(signal))    signal    = D%MuSignal(0)
  IF (PRESENT(signal_si)) signal_si = D%MuSignalSI(0)
  IF (PRESENT(signal_sd)) signal_sd = D%MuSignalSD(0)
  
  ! Signal rates
  IF (PRESENT(rate))    rate    = D%R(0)
  IF (PRESENT(rate_si)) rate_si = D%Rsi(0)
  IF (PRESENT(rate_sd)) rate_sd = D%Rsd(0)
  
  ! Bins
  IF (PRESENT(Nbins)) Nbins = Neff
  
  ! Signal events by bin
  IF (PRESENT(binsignal)) THEN
    ALLOCATE(binsignal(Neff))
    binsignal = D%MuSignal(1:Neff)
  END IF
  IF (PRESENT(binsignal_si)) THEN
    ALLOCATE(binsignal_si(Neff))
    binsignal_si = D%MuSignalSI(1:Neff)
  END IF
  IF (PRESENT(binsignal_sd)) THEN
    ALLOCATE(binsignal_sd(Neff))
    binsignal_sd = D%MuSignalSD(1:Neff)
  END IF
  
  ! Signal rates by bin
  IF (PRESENT(binrate)) THEN
    ALLOCATE(binrate(Neff))
    binrate = D%MuSignal(1:Neff)
  END IF
  IF (PRESENT(binrate_si)) THEN
    ALLOCATE(binrate_si(Neff))
    binrate_si = D%MuSignalSI(1:Neff)
  END IF
  IF (PRESENT(binrate_sd)) THEN
    ALLOCATE(binrate_sd(Neff))
    binrate_sd = D%MuSignalSD(1:Neff)
  END IF
  
END SUBROUTINE


! ----------------------------------------------------------------------
! Calculates various rate quantities using the passed WIMP and halo.
! This is the main routine intended for putting together the WIMP
! rate calculations and must be called any time the WIMP mass
! and/or couplings are modified.
! 
! Optional input/output argument:
!   D           The DetectorStruct structure containing detector
!               data, also where calculated rate data will be placed.
!               If not given, a default structure (internally stored)
!               will be used.
! 
SUBROUTINE CalcRates(D,WIMP,Halo)
  IMPLICIT NONE
  TYPE(DetectorStruct), INTENT(INOUT) :: D
  TYPE(WIMPStruct), INTENT(IN) :: WIMP
  TYPE(HaloStruct), INTENT(IN) :: Halo
  INTEGER :: Kiso,KE,Keff,Neff0
  REAL*8 :: alphasi(-1:1),alphasd(-1:1)
  ! Constant used to convert units:
  !   s / (cm^3 km GeV^4)  -->  cpd/kg/keV
  ! Includes factor of hbar^2 c^4.
  REAL*8, PARAMETER :: TO_CPD_KG_KEV = 1.695e14
  
  !---------------------------------------------
  ! Differential rate is given by:
  !   dR/dE(E) = 1/(2m\mu^2) \sigma(q) \rho \eta(vmin)
  ! where q(E) = \sqrt{2ME}, vmin(E) = \sqrt{ME/2\mu^2}, and:
  !   \sigma(E) = \mu^2 (hbar c)^2 [W(1,E)*Gp^2 + W(0,E)*Gp*Gn + W(-1,E)*Gn^2]
  ! This gives:
  !   dR/dE(E) = 1/2m \rho \eta(E) [W(1)*Gp^2 + W(0)*Gp*Gn + W(-1)*Gn^2]
  ! The total differential rate is a mass-fraction weighted sum of
  ! dR/dE over different isotopes and is summed over spin-independent
  ! (SI) and spin-dependent (SD) contributions.
  ! 
  ! The weighted form factors W are defined as:
  !   Wsi(+1,E) = (1/pi) Z^2 F^2(E)        ! SI proton
  !   Wsi( 0,E) = (1/pi) 2*Z*(A-Z) F^2(E)  ! SI crossterm
  !   Wsi(-1,E) = (1/pi) (A-Z)^2 F^2(E)    ! SI neutron
  !   Wsd(+1,E) = 4/(2J+1) Spp(E)          ! SD proton
  !   Wsd( 0,E) = 4/(2J+1) Spn(E)          ! SD crossterm
  !   Wsd(-1,E) = 4/(2J+1) Snn(E)          ! SD neutron
  ! The above definitions give for the conventional SI and SD
  ! cross-sections:
  !   \sigma(E) = \mu^2 (hbar c)^2 [W(+1,E)*Gp^2 + W(0,E)*Gp*Gn + W(-1,E)*Gn^2]
  ! where Gp and Gn are the effective proton and neutron couplings
  ! in units of [GeV^-2] and \mu is the reduced mass.  In terms of
  ! more commonly used notation:
  !   SI (scalar):        G = 2f
  !   SD (axial-vector):  G = 2\sqrt{2} G_F a
  ! where G, f, and a have 'p' and 'n' subscripts.
  ! 
  ! These spin-independent and spin-dependent weighted form factors are
  ! tabulated over E for each of the isotopes, so the Wsi and Wsd arrays
  ! are of size of size [-1:1,1:NE,1:Niso] where NE is the number of
  ! tabulation points (energies) and Niso is the number of isotopes.
  ! As these factors are independent of the WIMP, they are calculated
  ! only once, during initialization (nothing needs to be done here).
  ! SD FORM FACTORS ARE ONLY IMPLEMENTED FOR XENON.  ALL OTHER SD FORM
  ! FACTORS ARE SET TO ZERO.
  ! 
  ! The mean inverse speed
  !   eta(E) = \int_{v>vmin} d^3v (1/v) f(v),
  ! with
  !   vmin(E) = \sqrt{M E/(2\mu^2)}
  ! the minimum WIMP speed capable of causing a nucleus to recoil with
  ! energy E, depends upon the WIMP mass, so it will need to be
  ! recalculated for every WIMP.  Both vmin and eta are tabulated over
  ! E for each of the isotopes, so the vmin and eta arrays are of size
  ! [1:NE,1:Niso].
  ! 
  ! When factoring in detector efficiencies and energy resolution, the
  ! total rate of observed events is given by:
  !   R = \int_0^\infty dE \phi(E) dR/dE(E)
  ! where \phi(E) is the efficiency/response factor defined as the
  ! fraction of events at a given energy E that will be both observed
  ! and fall within some particular S1 analysis range.  The actual
  ! observed spectrum is dR/dS1(S1), NOT the quantity \phi(E) dR/dE(E);
  ! the former cannot be directly calculated using \phi(E).  However,
  ! the definition of \phi(E) allows the integral of each spectrum to
  ! be related as we can also write the rate as:
  !   R = \int_{S1min}^{S1max} dS1 dR/dS1(S1)
  ! The \phi(E) are calculated for _specific_ ranges of S1.  As these
  ! efficiencies are independent of the WIMP, they do not need to be
  ! recalculated, which is good, since this code cannot calculate
  ! them.  Instead, the the efficiency curve(s) are loaded from a file
  ! (or uses the internal default) during the initialization stage.
  ! 
  ! The efficiency file should contain a tabulated efficiency for the
  ! full desired S1 range, but can optionally include additional
  ! tabulated efficiencies for sub-intervals of that full range.
  ! For example, LUX 2013 analyzed the S1 range [2,30] and observed one
  ! event at S1=3.1 below their nuclear recoil mean.  For the maximum
  ! gap method of determining exclusion constraints in the presence of
  ! unknown backgrounds, the expected rate in intervals between events
  ! is used, so the rate in the S1 ranges [2,3.1] and [3.1,20] must be
  ! calculated.  The default efficiencies are the \phi(E) efficiencies
  ! for these three ranges: the full range and the two sub-intervals.
  ! The tabulated efficiency array (called 'eff') is an array of size
  ! [1:NE,0:Neff] where Neff is the number of sub-intervals (0 if
  ! separate interval efficiencies are not to be used or were not
  ! available).  The second index is that of the sub-interval, with 0
  ! used for the total range.
  ! 
  ! To avoid having to deal with a variety of interpolation issues,
  ! everything is tabulated at the same energies E, given in the array
  ! E of size [1:NE].  The tabulation is taken from that of the
  ! provided efficiency file (or internal default) since efficiencies
  ! cannot be recalculated with this program (everything else can).
  ! 
  ! In the calculations below, we calculate the differential rates
  ! dR/dE(E) for each of the isotopes and each of the SI & SD couplings.
  ! Then we perform an efficiency-weighted integral over E and sum
  ! contributions together to get the total rates and expected events.
  ! The separate coupling spectra are saved to allow for possible
  ! further inspection.
  ! 
  ! To be clear, the following are required to determine the full rate:
  !   * Sum over SI & SD contributions (with a Gp^2, Gn^2, and Gp*Gn
  !     term for each)
  !   * Mass fraction weighted sum over isotope contributions
  !   * Integration over energy
  !   * All of the above possibly for multiple analysis intervals
  ! Though not necesssary, various intermediate quantities are kept
  ! to allow for possible further inspection (aside from individual
  ! isotope contributions, which are not kept).
  ! Note there is no required order for performing the sums and
  ! integrations; we take advantage of this to perform the SI/SD
  ! component sum _last_ so that the rate is very easily recalculated
  ! for different couplings.  Notably, the two cross-sections are
  ! treated as:
  !    \sigma(E) = (Gp/Gp0)^2 \sigmapp0(E) + (Gn/Gn0)^2 \sigmann0(E)
  !                + (Gp*Gn/Gp0*Gn0) \sigmapn0(E)
  ! where
  !    \sigmapp0(E) = \mu^2 (hbar c)^2 W(+1,E) * Gp0^2
  !    \sigmann0(E) = \mu^2 (hbar c)^2 W(-1,E) * Gn0^2
  !    \sigmapn0(E) = \mu^2 (hbar c)^2 W( 0,E) * Gp0*Gn0
  ! and Gp0 & Gn0 are reference couplings, which we take to be those
  ! that yield WIMP-nucleon cross-sections of 1 pb.  In that case,
  ! we can also write:
  !    \sigma(E) = (sigmap/[pb])*\sigmapp0(E) + (sigman/[pb])*\sigmann0(E)
  !                +/- \sqrt{(sigmap/[pb])(sigman/[pb])} \sigmapn0(E)
  ! with the sign of the cross-term equal to the sign of Gp*Gn.  The
  ! quantities below with names ending with 'si0' and 'sd0' have an
  ! index [-1:1] identical to that described for the W terms above,
  ! corresponds to the that quantity calculated for sigma(E) =
  ! sigmapp0(E) [+1], sigmann0(E) [-1], or sigmapn0(E) [0].
  !---------------------------------------------
    
  ! Update mean inverse speed.
  D%vmin = EToVmin(D%NE,D%E,WIMP%m,D%Niso,D%Miso)
  D%eta  = MeanInverseSpeed(D%NE,D%Niso,D%vmin,Halo)
  
  ! Below, we will calculate reference rates for reference couplings
  ! GpSI0, GnSI0, GpSD0, & GnSD0 defined such that WIMP-nucleon
  ! cross-sections are 1 pb.  To convert between reference rates and
  ! actual rates, we multiply by alpha = G^2/G0^2.  Indices of alpha
  ! match that of the W arrays.
  alphasi(+1) = (WIMP%GpSI/WIMP%GpSI0)**2                      ! SI proton
  alphasi( 0) = (WIMP%GpSI/WIMP%GpSI0)*(WIMP%GnSI/WIMP%GnSI0)  ! SI crossterm
  alphasi(-1) = (WIMP%GnSI/WIMP%GnSI0)**2                      ! SI neutron
  alphasd(+1) = (WIMP%GpSD/WIMP%GpSD0)**2                      ! SD proton
  alphasd( 0) = (WIMP%GpSD/WIMP%GpSD0)*(WIMP%GnSD/WIMP%GnSD0)  ! SD crossterm
  alphasd(-1) = (WIMP%GnSD/WIMP%GnSD0)**2                      ! SD neutron
  
  ! Reference differential rates.
  ! Keep separate contributions from Gp^2, Gn^2, and Gp*Gn terms.
  D%dRdEsi0 = 0d0
  D%dRdEsd0 = 0d0
  ! Mass-fraction weighted sum over isotopes.
  DO Kiso = 1,D%Niso
    ! SI proton term
    D%dRdEsi0(+1,:) = D%dRdEsi0(+1,:)                                 &
        +  D%fiso(Kiso) / (2 * WIMP%m)                                 &
           * Halo%rho * D%eta(:,Kiso)                                  &
           * D%Wsi(+1,:,Kiso) * WIMP%GpSI0**2                          &
           * TO_CPD_KG_KEV
    ! SI cross term
    D%dRdEsi0( 0,:) = D%dRdEsi0( 0,:)                                 &
        +  D%fiso(Kiso) / (2 * WIMP%m)                                 &
           * Halo%rho * D%eta(:,Kiso)                                  &
           * D%Wsi( 0,:,Kiso) * WIMP%GpSI0*WIMP%GnSI0                  &
           * TO_CPD_KG_KEV
    ! SI neutron term
    D%dRdEsi0(-1,:) = D%dRdEsi0(-1,:)                                 &
        +  D%fiso(Kiso) / (2 * WIMP%m)                                 &
           * Halo%rho * D%eta(:,Kiso)                                  &
           * D%Wsi(-1,:,Kiso) * WIMP%GnSI0**2                          &
           * TO_CPD_KG_KEV
    ! SD proton term
    D%dRdEsd0(+1,:) = D%dRdEsd0(+1,:)                                 &
        +  D%fiso(Kiso) / (2 * WIMP%m)                                 &
           * Halo%rho * D%eta(:,Kiso)                                  &
           * D%Wsd(+1,:,Kiso) * WIMP%GpSD0**2                          &
           * TO_CPD_KG_KEV
    ! SD cross term
    D%dRdEsd0( 0,:) = D%dRdEsd0( 0,:)                                 &
        +  D%fiso(Kiso) / (2 * WIMP%m)                                 &
           * Halo%rho * D%eta(:,Kiso)                                  &
           * D%Wsd( 0,:,Kiso) * WIMP%GpSD0*WIMP%GnSD0                  &
           * TO_CPD_KG_KEV
    ! SD neutron term
    D%dRdEsd0(-1,:) = D%dRdEsd0(-1,:)                                 &
        +  D%fiso(Kiso) / (2 * WIMP%m)                                 &
           * Halo%rho * D%eta(:,Kiso)                                  &
           * D%Wsd(-1,:,Kiso) * WIMP%GnSD0**2                          &
           * TO_CPD_KG_KEV
  END DO
  
  ! Number of intervals/bins to do calculations for.
  ! If intervals=.FALSE. then efficiency index is over [0:0]
  ! (total only).
  IF (D%intervals) THEN
    Neff0 = D%Neff
  ELSE
    Neff0 = 0
  END IF
  
  ! Integrate (efficiency-weighted) to find total rates.
  ! Uses a simple trapezoidal integration.
  D%Rsi0 = 0d0
  D%Rsd0 = 0d0
  ! Cycle over E bins and efficiency curves.
  DO KE = 1,D%NE-1
    DO Keff = 0,Neff0
      D%Rsi0(:,Keff) = D%Rsi0(:,Keff)                                 &
          + 0.5d0 * (D%E(KE+1) - D%E(KE))                             &
            * (D%eff0(KE,Keff)*D%dRdEsi0(:,KE)                        &
               + D%eff0(KE+1,Keff)*D%dRdEsi0(:,KE+1))
      D%Rsd0(:,Keff) = D%Rsd0(:,Keff)                                 &
          + 0.5d0 * (D%E(KE+1) - D%E(KE))                             &
            * (D%eff0(KE,Keff)*D%dRdEsd0(:,KE)                        &
               + D%eff0(KE+1,Keff)*D%dRdEsd0(:,KE+1))
    END DO
  END DO
  
  ! Rates for actual couplings.  Here do separate SI and SD rates
  ! as well as total.
  D%dRdEsi(:) =  alphasi(+1) * D%dRdEsi0(+1,:)                        &
                + alphasi( 0) * D%dRdEsi0( 0,:)                        &
                + alphasi(-1) * D%dRdEsi0(-1,:)
  D%dRdEsd(:) =  alphasd(+1) * D%dRdEsd0(+1,:)                        &
                + alphasd( 0) * D%dRdEsd0( 0,:)                        &
                + alphasd(-1) * D%dRdEsd0(-1,:)
  D%dRdE(:)   = D%dRdEsi(:) + D%dRdEsd(:)
  D%Rsi(:) =  alphasi(+1) * D%Rsi0(+1,:)                              &
             + alphasi( 0) * D%Rsi0( 0,:)                              &
             + alphasi(-1) * D%Rsi0(-1,:)
  D%Rsd(:) =  alphasd(+1) * D%Rsd0(+1,:)                              &
             + alphasd( 0) * D%Rsd0( 0,:)                              &
             + alphasd(-1) * D%Rsd0(-1,:)
  D%R(:)   = D%Rsi(:) + D%Rsd(:)
  
  ! Average expected event (components at reference couplings)
  D%MuSignalSI0(:,:) = D%exposure * D%Rsi0(:,:)
  D%MuSignalSD0(:,:) = D%exposure * D%Rsd0(:,:)
  
  ! Average expected events.
  ! Could also achieve this by e.g. exposure * Rsi(:).
  D%MuSignalSI(:) =  alphasi(+1) * D%MuSignalSI0(+1,:)                 &
                    + alphasi( 0) * D%MuSignalSI0( 0,:)                &
                    + alphasi(-1) * D%MuSignalSI0(-1,:)
  D%MuSignalSD(:) =  alphasd(+1) * D%MuSignalSD0(+1,:)                 &
                    + alphasd( 0) * D%MuSignalSD0( 0,:)                &
                    + alphasd(-1) * D%MuSignalSD0(-1,:)
  D%MuSignal(:)   = D%MuSignalSI(:) + D%MuSignalSD(:)
  
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


! ----------------------------------------------------------------------
! Returns the average expected number of spin-independent signal events
! for the current WIMP.
! 
! Required input argument:
!   D           A DetectorStruct containing rates (CalcRates(D) must
!               have already been called).
! 
FUNCTION GetSignalSI(D) RESULT(s)
  IMPLICIT NONE
  REAL*8 :: s
  TYPE(DetectorStruct), INTENT(IN) :: D
  CALL GetRates(D,signal_si=s)
END FUNCTION

!-----------------------------------------------------------------------
! C/C++ wrapper for DDCalc_SignalSI
!
REAL(KIND=C_DOUBLE) FUNCTION C_DDCalc_SignalSI(DetectorIndex) &
 BIND(C,NAME='C_DDRates_ddcalc_signalsi')
  USE ISO_C_BINDING, only: C_DOUBLE, C_INT
  IMPLICIT NONE
  INTEGER(KIND=C_INT), INTENT(IN) :: DetectorIndex
  IF (.NOT. ASSOCIATED(Detectors(DetectorIndex)%p)) stop 'Invalid detector index given to C_DDCalc_SignalSI'
  C_DDCalc_SignalSI = REAL(GetSignalSI(Detectors(DetectorIndex)%p),KIND=C_DOUBLE)
END FUNCTION


! ----------------------------------------------------------------------
! Returns the average expected number of spin-dependent signal events
! for the current WIMP.
! 
! Required input argument:
!   D           A DetectorStruct containing rates (CalcRates(D) must
!               have already been called).
! 
FUNCTION GetSignalSD(D) RESULT(s)
  IMPLICIT NONE
  REAL*8 :: s
  TYPE(DetectorStruct), INTENT(IN) :: D
  CALL GetRates(D,signal_sd=s)
END FUNCTION

!-----------------------------------------------------------------------
! C/C++ wrapper for DDCalc_SignalSD
!
REAL(KIND=C_DOUBLE) FUNCTION C_DDCalc_SignalSD(DetectorIndex) &
 BIND(C,NAME='C_DDRates_ddcalc_signalsd')
  USE ISO_C_BINDING, only: C_DOUBLE, C_INT
  IMPLICIT NONE
  INTEGER(KIND=C_INT), INTENT(IN) :: DetectorIndex
  IF (.NOT. ASSOCIATED(Detectors(DetectorIndex)%p)) stop 'Invalid detector index given to C_DDCalc_SignalSD'
  C_DDCalc_SignalSD = REAL(GetSignalSD(Detectors(DetectorIndex)%p),KIND=C_DOUBLE)
END FUNCTION


END MODULE
