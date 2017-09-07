MODULE DDTypes

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! DDTypes
!    Global types for DDCalc
!
! If changed, recompile e.g. via
!   rm -f ddtypes.mod && cd src/ && gfortran -c DDTypes.f90 && rm -f DDTypes.o && mv ddtypes.mod .. && cd .. && make clean && make
! from the DDCalc main directory.
!
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


! WIMP Structure
TYPE, PUBLIC :: WIMPStruct
  REAL*8 :: m   ! WIMP mass [GeV]
  ! WIMP type
  ! Currently the following possibilities are implemented:
  ! 'SIonly':      WIMP with only spin-independent interactions
  ! 'SDonly':      WIMP with only spin-dependent interactions
  ! 'SISD':        WIMP with spin-independent and spin-dependent 
  !                interactions
  ! 'HiggsPortal': Fermionic WIMP with scalar and pseudoscalar
  !                couplings to the Higgs boson
  CHARACTER(LEN=24) :: DMtype = ''
  ! List of WIMP parameters. 
  ! The length and meaning of this list depends on the WIMP type:
  ! 'SIonly':      (fp, fn)
  ! fp: SI DM-proton coupling, fn: SI DM-neutron coupling
  ! 'SDonly':      (ap, an)
  ! ap: SD DM-proton coupling, an: SD DM-neutron coupling
  ! 'SISD':        (fp, fn, ap, an)
  ! fp: SI DM-proton coupling, fn: SI DM-neutron coupling
  ! ap: SD DM-proton coupling, an: SD DM-neutron coupling
  ! 'HiggsPortal': (fsp, fsn, app, apn)
  ! fsp: scalar DM-proton coupling, fsn: scalar DM-neutron coupling
  ! app: pseudoscalar DM-proton coupling, apn: pseudoscalar DM-neutron coupling
  REAL*8, ALLOCATABLE  :: params(:)
  ! Number of parameters
  INTEGER :: Nparams = 1
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
  
  ! Flag that indicates whether the detector has been successfully initialized.
  ! If this remains false after calling SetDetector, it means that some of the 
  ! required input arguments have not been specified.
  LOGICAL :: InitSuccess = .FALSE.

  ! Flag which tells other routines whether this detector is currently set up for
  ! MaxGap, TotalPoisson or BinnedPoisson. This flag should not be edited manually from the outside,
  ! it's value is solely determined by the presence and values of Nevents_tot and Nevents_bin
  ! -1: no correct statistics method defined
  ! 0: TotalPoisson
  ! 1: BinnedPoisson
  ! 2: MaxGap
  INTEGER :: StatisticFlag = -1

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
  
  ! Exposure -----------------------------------
  ! Detector fiducial mass [kg]
  REAL*8 :: mass = -1d0
  
  ! Detector exposure time [day]
  REAL*8 :: time = -1d0
  
  ! Total detector exposure [kg*day]
  REAL*8 :: exposure = -1d0
  
  ! Events -------------------------------------
  ! Observed number of events
  INTEGER, ALLOCATABLE :: Nevents(:) ! Nevents(0) is the total number of observed events, Nevents(1:) the binned observed backgroundevents
  
  ! Average expected background events
  REAL*8, ALLOCATABLE  :: Backgr(:)   ! Backgr(0) is the total number of background events, Nevents(1:) the binned background events
 
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
  
  ! Number of S1 bins/intervals with efficiencies.
  ! Will calculate rates for each bin/interval plus total.
  INTEGER :: Nbins = -1
  
  ! Array of size [1:Niso,1:NE,0:Nbins] with the third index for the S1
  ! bin/interval (zero for full range)
  REAL*8, ALLOCATABLE :: eff(:,:,:)
  
  ! Indicates if rates for intervals/bins are to also be calculated
  ! in addition to the total rate.  Needed for maximum gap analysis,
  ! but unnecessary for likelihood.  This flag is ignored if the
  ! efficiencies for intervals/bins are not provided.
  LOGICAL :: intervals = .FALSE.
  
  
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

  ! structure for the differential rates,
  ! separately for each isotope of the experiment.
  ! These represent rates before efficiency cuts, but they already
  ! do include the mass fraction of the corresponding isotope.
  ! Units of this are cpd/kg/keV
  ! Array is of size [1:NE,1:Niso].
  REAL*8, ALLOCATABLE :: dRdEiso(:,:)
  
  ! Integrated rate --------
  ! Efficiency-corrected rates.  Array is of size
  ! [0:Nbins] with the index being that of the S1 bin/interval
  ! efficiency curve used in the integral (0 for full range).
  ! [cpd/kg]
  REAL*8, ALLOCATABLE :: R(:)
  
  ! Events -------------------------------------
  ! Expected number of signal events.  Arrays of size [0:Nbins].
  REAL*8, ALLOCATABLE :: MuSignal(:) 
  
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
