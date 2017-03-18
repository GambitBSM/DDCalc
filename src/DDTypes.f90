MODULE DDTypes

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! DDTypes
!    Global types for DDCalc
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

IMPLICIT NONE
PRIVATE


! Commandline argument structure
TYPE, PUBLIC :: ArgumentStruct
  ! Currently, only the parameters are used (options are parsed directly).
  INTEGER :: Noptions    = -1
  INTEGER :: Nparameters = -1
  !Arguments that begin with '--', e.g. --<flag> or --<flag>=<value>
  CHARACTER*256, ALLOCATABLE :: options(:)
  !Arguments that are not options
  CHARACTER*256, ALLOCATABLE :: parameters(:)
  !Conversion of parameters to floating point (if possible)
  REAL*8, ALLOCATABLE :: values(:)
END TYPE 


! Structure to contain tabulated rates as a function of energy.
TYPE, PUBLIC :: DetectorRateStruct
  
  ! Integrated rates (reference couplings) -----
  ! Reference integrated rates calculated at \sigma = 1 pb for each
  ! coupling, with the efficiency-weighted integral done separately
  ! for each available efficiency curve.  Arrays are of size
  ! [-1:1,0:Neff].  The first index is for the proton & neutron
  ! components at sigma = 1 pb in each case and the second index is for
  ! the S1 bin/interval (0 for full range). [cpd/kg]
  REAL*8, ALLOCATABLE :: Rsi0(:,:),Rsd0(:,:)
  
  ! Integrated rates (actual couplings) --------
  ! Efficiency-corrected rates at given couplings.  Arrays are of size
  ! [0:Neff] with the index being that of the S1 bin/interval
  ! efficiency curve used in the integral (0 for full range).  Given
  ! separately for SI and SD components as well as the SI+SD total.
  ! [cpd/kg]
  REAL*8, ALLOCATABLE :: Rsi(:),Rsd(:),R(:)
  
  ! Events -------------------------------------
  ! Expected number of signal events at reference couplings.
  ! Arrays of size [-1:1,0:Neff].
  REAL*8, ALLOCATABLE :: MuSignalSI0(:,:),MuSignalSD0(:,:)
  ! Expected number of signal events.  Arrays of size [0:Neff].
  REAL*8, ALLOCATABLE :: MuSignalSI(:),MuSignalSD(:),MuSignal(:)
  
  ! Average expected background events
  REAL*8 :: MuBackground = 0d0
  
  ! Observed number of events
  INTEGER :: Nevents = -1
  
END TYPE


! WIMP Structure
TYPE, PUBLIC :: WIMPStruct
  REAL*8 :: m   ! WIMP mass [GeV]
  ! Effective couplings to the proton and neutron in units of [GeV^-2].
  ! In terms of more commonly used notation:
  !   SI (scalar):        G = 2f
  !   SD (axial-vector):  G = 2\sqrt{2} G_F a
  ! where G, f, and a have 'p' and 'n' subscripts.
  REAL*8 :: GpSI,GnSI,GpSD,GnSD  
  ! Couplings that yield \sigma = 1 pb in each case.
  REAL*8 :: GpSI0,GnSI0,GpSD0,GnSD0
END TYPE

! Structure containing halo parameters
TYPE, PUBLIC :: HaloStruct
  ! Galactic motions ---------------------------
  ! Local galactic disk rotation speed [km/s]
  REAL*8 :: vrot = 235d0
  ! Local standard of rest velocity vector [km/s], defined relative to
  ! galactic rest frame.
  REAL*8 :: vlsr(3) = (/ 0d0, 235d0, 0d0 /)
  ! Sun's peculiar velocity vector [km/s], defined relative to local
  ! standard of rest.
  REAL*8 :: vpec(3) = (/ 11d0, 12d0, 7d0 /)
  ! Sun's velocity vector [km/s], defined relative to galactic rest
  ! frame.
  REAL*8 :: vsun(3) = (/ 0d0, 235d0, 0d0 /) + (/ 11d0, 12d0, 7d0 /)
  ! Sun's speed (or observer's speed) [km/s], defined relative to
  ! galactic rest frame.
  REAL*8 :: vobs = SQRT(11d0**2 + (235d0+12d0)**2 + 7d0**2)
  
  ! Local DM density ---------------------------
  ! Local dark matter density [GeV/cm^3]
  REAL*8 :: rho = 0.4d0
  
  ! DM distribution (SHM) ----------------------
  ! Truncated Maxwell-Boltzmann ("MB") distribution.
  
  ! Bulk velocity of the dark matter [km/s] (i.e. the velocity of
  ! the MB rest frame), defined relative to the galactic rest frame.
  REAL*8 :: vbulk(3) = (/ 0d0, 0d0, 0d0 /)
  ! Most probable speed [km/s] in the MB rest frame.
  REAL*8 :: v0 = 235d0
  ! Escape speed [km/s] in the MB rest frame.
  REAL*8 :: vesc = 550d0
  
  ! DM distribution (tabulated) ----------------
  ! Instead of being calculated for SHM above, mean inverse speed
  ! can be given explicitly as a tabulation.
  LOGICAL :: tabulated = .FALSE.
  CHARACTER(LEN=1024) :: eta_file = ''
  INTEGER :: Nvmin = -1
  REAL*8, ALLOCATABLE :: vmin(:)
  REAL*8, ALLOCATABLE :: eta(:)
END TYPE


! Structure to contain detector characteristics and tabulated rates as
! a function of energy.
TYPE, PUBLIC :: DetectorStruct
  
  ! Label --------------------------------------
  ! Label for the experimental result contained in this structure.
  ! Must be at most 12 characters as it will be used in column headers.
  CHARACTER(LEN=12) :: label = ''
  ! More detailed description.
  CHARACTER(LEN=1024) :: description = ''
  
  ! FUTURE IMPLEMENTATION ----------------------
  ! Flag to indicate if array sizes and values are outdated
  ! and need to be reinitialized.
  LOGICAL :: stale = .TRUE.
  
  ! Detector parameters (future use)
  !TYPE(DetectorParametersStruct) :: parameters
  
  ! Detector efficiencies  (future use)
  !TYPE(DetectorEfficiencyStruct) :: efficiency
  
  ! Detector rates  (future use)
  !TYPE(DetectorRateStruct) :: rates
  
  ! Exposure -----------------------------------
  ! Detector fiducial mass [kg]
  REAL*8 :: mass = 118d0
  
  ! Detector exposure time [day]
  REAL*8 :: time = 85.3d0
  
  ! Total detector exposure [kg*day]
  REAL*8 :: exposure = 118d0*85.3d0
  
  ! Events -------------------------------------
  ! Observed number of events
  INTEGER :: Nevents = -1
  
  ! Average expected background events
  REAL*8 :: MuBackground = 0d0
  
  ! Isotopes -----------------------------------
  ! Number of isotopes
  INTEGER :: Niso = -1
  
  ! Detector isotopes, their mass fractions, and nuclear masses [GeV]
  INTEGER, ALLOCATABLE :: Ziso(:)
  INTEGER, ALLOCATABLE :: Aiso(:)
  REAL*8, ALLOCATABLE  :: fiso(:)
  REAL*8, ALLOCATABLE  :: Miso(:)  ! Calculated internally
  
  ! Tabulation ---------------------------------
  ! Number of tabulation points (energies).
  INTEGER :: NE = -1
  
  ! Tabulated energies [keV].  Array of size [1:NE].
  REAL*8, ALLOCATABLE :: E(:)
  
  ! Cached energy tabulation as changes to E array will be made
  ! for internal purposes.
  REAL*8, ALLOCATABLE :: E_cache(:)
  
  ! Efficiencies -------------------------------
  ! Tabulated detection efficiencies.
  ! File containing efficiencies
  CHARACTER(LEN=1024) :: eff_file = ''
  
  ! Number of energies for efficiency tabulation.
  INTEGER :: NEeff = -1
  
  ! Energies for efficiency tabulation [keV].  Array of size [1:NEeff].
  REAL*8, ALLOCATABLE :: Eeff(:)
  
  ! Number of S1 bins/intervals with efficiencies.
  ! Will calculate rates for each bin/interval plus total.
  INTEGER :: Neff = -1
  
  ! Array of size [1:NEeff,0:Neff] with the second index for the S1
  ! bin/interval (zero for full range)
  REAL*8, ALLOCATABLE :: eff(:,:)
  
  ! Indicates if rates for intervals/bins are to also be calculated
  ! in addition to the total rate.  Needed for maximum gap analysis,
  ! but unnecessary for likelihood.  This flag is ignored if the
  ! efficiencies for intervals/bins are not provided.
  LOGICAL :: intervals = .TRUE.
  
  ! Array of size [1:NE,0:Neff] or [1:NE,0:0] containing efficiencies
  ! retabulated to that of the E array.  Sub-interval efficiencies are
  ! dropped if intervals=.FALSE.
  REAL*8, ALLOCATABLE :: eff0(:,:)
  
  ! Form factors -------------------------------
  ! Tabulated spin-independent or spin-dependent form factors combined
  ! with prefactors.  Arrays of size [-1:1,1:NE,1:Niso].  Defined as
  ! [unitless]:
  !   Wsi(+1,:,:) = (1/pi) Z^2 F^2(q)        ! SI proton
  !   Wsi( 0,:,:) = (1/pi) 2*Z*(A-Z) F^2(q)  ! SI crossterm
  !   Wsi(-1,:,:) = (1/pi) (A-Z)^2 F^2(q)    ! SI neutron
  !   Wsd(+1,:,:) = 4/(2J+1) Spp(q)          ! SD proton
  !   Wsd( 0,:,:) = 4/(2J+1) Spn(q)          ! SD crossterm
  !   Wsd(-1,:,:) = 4/(2J+1) Snn(q)          ! SD neutron
  ! The above definitions give for the conventional SI and SD
  ! cross-sections:
  !   \sigma(q) = \mu^2 (hbar c)^2 [W(1)*Gp^2 + W(0)*Gp*Gn + W(-1)*Gn^2]
  ! where Gp and Gn are the effective proton and neutron couplings
  ! in units of [GeV^-2] and \mu is the reduced mass.  In terms of
  ! more commonly used notation:
  !   SI (scalar):        G = 2f
  !   SD (axial-vector):  G = 2\sqrt{2} G_F a
  ! where G, f, and a have 'p' and 'n' subscripts.  While form factors
  ! are often a function of the momentum transfer, we tabulate them
  ! here as a function of recoil energy E = q^2/2M.
  ! NOTE: Need only be calculated once.
  REAL*8, ALLOCATABLE :: Wsi(:,:,:),Wsd(:,:,:)
  
  ! Halo ---------------------------------------
  ! The minimum velocity for producing a recoil of energy E, given
  ! by vmin = sqrt{M E/(2\mu^2)} [km/s].
  ! Array of size [1:NE,1:Niso] that needs to be recalculated when the
  ! WIMP mass changes.
  REAL*8, ALLOCATABLE :: vmin(:,:)
  
  ! Tabulated mean inverse speed (eta) [s/km] at the above vmin.
  REAL*8, ALLOCATABLE :: eta(:,:)
  
  ! Differential rates (reference couplings) ---
  ! Reference differential rates calculated at \sigma = 1 pb for each
  ! coupling and tabulated by energy.  Differential arrays are of size
  ! [-1:1,1:NE] and are given separately for SI and SD couplings.  The
  ! first index is for the proton & neutron components at sigma = 1 pb
  ! in each case and the second index is that of the energy (given by
  ! the E array at the same index).  These represent raw rates prior to
  ! any efficiency cuts. [cpd/kg/keV]
  REAL*8, ALLOCATABLE :: dRdEsi0(:,:),dRdEsd0(:,:)
  
  ! Integrated rates (reference couplings) -----
  ! Reference integrated rates calculated at \sigma = 1 pb for each
  ! coupling, with the efficiency-weighted integral done separately
  ! for each available efficiency curve.  Arrays are of size
  ! [-1:1,0:Neff].  The first index is for the proton & neutron
  ! components at sigma = 1 pb in each case and the second index is for
  ! the S1 bin/interval (0 for full range). [cpd/kg]
  REAL*8, ALLOCATABLE :: Rsi0(:,:),Rsd0(:,:)
  
  ! Differential rates (actual couplings) ------
  ! Differential rates at given couplings, tabulated by energy.
  ! Arrays are of size [1:NE], given separately for SI and SD couplings
  ! as well as the SI+SD total. [cpd/kg/keV]
  REAL*8, ALLOCATABLE :: dRdEsi(:),dRdEsd(:),dRdE(:)
  
  ! Integrated rates (actual couplings) --------
  ! Efficiency-corrected rates at given couplings.  Arrays are of size
  ! [0:Neff] with the index being that of the S1 bin/interval
  ! efficiency curve used in the integral (0 for full range).  Given
  ! separately for SI and SD components as well as the SI+SD total.
  ! [cpd/kg]
  REAL*8, ALLOCATABLE :: Rsi(:),Rsd(:),R(:)
  
  ! Events -------------------------------------
  ! Expected number of signal events at reference couplings.
  ! Arrays of size [-1:1,0:Neff].
  REAL*8, ALLOCATABLE :: MuSignalSI0(:,:),MuSignalSD0(:,:)
  ! Expected number of signal events.  Arrays of size [0:Neff].
  REAL*8, ALLOCATABLE :: MuSignalSI(:),MuSignalSD(:),MuSignal(:)
  
END TYPE


! <<<<THE FOLLOWING STRUCTURES ARE FOR FUTURE USE ONLY>>>>


! Structure to contain tabulated detection efficiencies as a
! function of energy, for the overall analysis range and possibly
! for subintervals/bins.
TYPE, PUBLIC :: DetectorEfficiencyStruct
  
  ! File containing efficiencies
  CHARACTER(LEN=1024) :: file = ''
  
  ! Number of tabulation points (energies).
  INTEGER :: NE = -1
  
  ! Tabulated energies [keV].  Array of size [1:NE].
  REAL*8, ALLOCATABLE :: E(:)
  
  ! Number of S1 bins/sub-intervals with efficiencies (does not
  ! include total interval). Will calculate rates for each bin/interval
  ! plus total.
  INTEGER :: Neff = -1
  
  ! Array of size [1:NE,0:Neff] with the second index for the S1
  ! bin/interval (zero for full range)
  REAL*8, ALLOCATABLE :: eff(:,:)
  
END TYPE


! Structure to contain various detector parameters.
TYPE, PUBLIC :: DetectorParametersStruct
  
  ! Exposure -----------------------------------
  ! Detector fiducial mass [kg]
  REAL*8 :: mass = 118d0
  
  ! Detector exposure time [day]
  REAL*8 :: time = 85.3d0
  
  ! Total detector exposure [kg*day]
  REAL*8 :: exposure = 118d0*85.3d0
  
  ! Isotopes -----------------------------------
  ! Number of isotopes
  INTEGER :: Niso = -1
  
  ! Detector isotopes, their mass fractions, and nuclear masses [GeV]
  INTEGER, ALLOCATABLE :: Ziso(:)
  INTEGER, ALLOCATABLE :: Aiso(:)
  REAL*8, ALLOCATABLE  :: fiso(:)
  REAL*8, ALLOCATABLE  :: Miso(:)  ! Calculated internally
  
END TYPE


! Structure to contain tabulated differential rates dR/dE as a function
! of energy.
TYPE, PUBLIC :: DetectorSpectraStruct
  
  ! Tabulation ---------------------------------
  ! Number of tabulation points (energies).
  ! NOTE: This tabulation is fixed to that used by the efficiency data.
  INTEGER :: NE = -1
  
  ! Tabulated energies [keV].  Array of size [1:NE].
  REAL*8, ALLOCATABLE :: E(:)
  
  ! Efficiencies -------------------------------
  ! Tabulated detection efficiencies.  Here tabulated at desired
  ! E for dR/dE calculations.
  
  ! Number of S1 bins/intervals with efficiencies.
  ! Will calculate rates for each bin/interval plus total.
  INTEGER :: Neff = -1
  
  ! Array of size [1:NE,0:Neff] with the second index for the S1
  ! bin/interval (zero for full range)
  REAL*8, ALLOCATABLE :: eff(:,:)
  
  ! Form factors -------------------------------
  ! Tabulated spin-independent or spin-dependent form factors combined
  ! with prefactors.  Arrays of size [-1:1,1:NE,1:Niso].  Defined as
  ! [unitless]:
  !   Wsi(+1,:,:) = (1/pi) Z^2 F^2(q)        ! SI proton
  !   Wsi( 0,:,:) = (1/pi) 2*Z*(A-Z) F^2(q)  ! SI crossterm
  !   Wsi(-1,:,:) = (1/pi) (A-Z)^2 F^2(q)    ! SI neutron
  !   Wsd(+1,:,:) = 4/(2J+1) Spp(q)          ! SD proton
  !   Wsd( 0,:,:) = 4/(2J+1) Spn(q)          ! SD crossterm
  !   Wsd(-1,:,:) = 4/(2J+1) Snn(q)          ! SD neutron
  ! The above definitions give for the conventional SI and SD
  ! cross-sections:
  !   \sigma(q) = \mu^2 (hbar c)^2 [W(1)*Gp^2 + W(0)*Gp*Gn + W(-1)*Gn^2]
  ! where Gp and Gn are the effective proton and neutron couplings
  ! in units of [GeV^-2] and \mu is the reduced mass.  In terms of
  ! more commonly used notation:
  !   SI (scalar):        G = 2f
  !   SD (axial-vector):  G = 2\sqrt{2} G_F a
  ! where G, f, and a have 'p' and 'n' subscripts.  While form factors
  ! are often a function of the momentum transfer, we tabulate them
  ! here as a function of recoil energy E = q^2/2M.
  ! NOTE: Need only be calculated once.
  REAL*8, ALLOCATABLE :: Wsi(:,:,:),Wsd(:,:,:)
  
  ! Halo ---------------------------------------
  ! The minimum velocity for producing a recoil of energy E, given
  ! by vmin = sqrt{M E/(2\mu^2)} [km/s].
  ! Array of size [1:NE,1:Niso] that needs to be recalculated when the
  ! WIMP mass changes.
  REAL*8, ALLOCATABLE :: vmin(:,:)
  
  ! Tabulated mean inverse speed (eta) [s/km] at the above vmin.
  REAL*8, ALLOCATABLE :: eta(:,:)
  
  ! Differential rates (reference couplings) ---
  ! Reference differential rates calculated at \sigma = 1 pb for each
  ! coupling and tabulated by energy.  Differential arrays are of size
  ! [-1:1,1:NE] and are given separately for SI and SD couplings.  The
  ! first index is for the proton & neutron components at sigma = 1 pb
  ! in each case and the second index is that of the energy (given by
  ! the E array at the same index).  These represent raw rates prior to
  ! any efficiency cuts. [cpd/kg/keV]
  REAL*8, ALLOCATABLE :: dRdEsi0(:,:),dRdEsd0(:,:)
  
  ! Differential rates (actual couplings) ------
  ! Differential rates at given couplings, tabulated by energy.
  ! Arrays are of size [1:NE], given separately for SI and SD couplings
  ! as well as the SI+SD total. [cpd/kg/keV]
  REAL*8, ALLOCATABLE :: dRdEsi(:),dRdEsd(:),dRdE(:)
  
END TYPE

!#### Internal caches of WIMPs, Halos and Detectors for C interface ####
TYPE, PUBLIC :: WIMPStructPtr
  TYPE(WIMPStruct), POINTER :: p
END TYPE
TYPE, PUBLIC :: HaloStructPtr
  TYPE(HaloStruct), POINTER :: p
END TYPE
TYPE, PUBLIC :: DetectorStructPtr
  TYPE(DetectorStruct), POINTER :: p
END TYPE
INTEGER, PUBLIC, PARAMETER :: Max_WIMPs = 1024
INTEGER, PUBLIC, PARAMETER :: Max_Halos = 1024
INTEGER, PUBLIC, PARAMETER :: Max_Detectors = 1024
TYPE(WIMPStructPtr),     DIMENSION(Max_WIMPs), PUBLIC :: WIMPs
TYPE(HaloStructPtr),     DIMENSION(Max_Halos), PUBLIC :: Halos
TYPE(DetectorStructPtr), DIMENSION(Max_Detectors), PUBLIC :: Detectors
INTEGER, PUBLIC :: N_WIMPs = 0, N_Halos = 0, N_Detectors = 0

END MODULE
