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

! WIMP Structure
TYPE, PUBLIC :: WIMPStruct
  REAL*8 :: m   ! WIMP mass [GeV]
  ! WIMP type
  ! Currently the following possibilities are implemented:
  ! 'SIonly':      WIMP with only spin-independent interactions
  ! 'SDonly':      WIMP with only spin-dependent interactions
  ! 'SISD':        WIMP with spin-independent and spin-dependent 
  !                interactions
  ! 'SILR':        WIMP with long-range interactions via a light mediator
  ! 'NREffectiveTheory': Effective WIMP defined via the non-relativistic DM-nucleon operators
  ! 'NREFT_CPT':   Effective WIMP defined via non-relativistic operators in chiral perturbation theory
  CHARACTER(LEN=24) :: DMtype = ''
  ! List of WIMP parameters. 
  ! The length and meaning of this list depends on the WIMP type:
  !
  ! 'SIonly':      (fp, fn)
  ! fp: SI DM-proton coupling, fn: SI DM-neutron coupling
  !
  ! 'SDonly':      (ap, an)
  ! ap: SD DM-proton coupling, an: SD DM-neutron coupling
  !
  ! 'SISD':        (fp, fn, ap, an)
  ! fp: SI DM-proton coupling, fn: SI DM-neutron coupling
  ! ap: SD DM-proton coupling, an: SD DM-neutron coupling
  !
  ! 'SILR':        (gp, gn, mmed)
  ! gp: SI DM-proton coupling, gn: SI DM-neutron coupling
  ! mmed: mediator mass
  ! In contrast to WIMP type 'SIonly' the couplings here are assumed to be dimensionless
  ! The correct mass dimension then results from the propagator factor 1/(q^2 + mmed^2)^2
  !
  ! 'NREffectiveTheory': 
  ! (DM spin, O1_0, O1_1, O1q2_0, O1q2_1, O3_0, O3_1, O4_0, O4_1, O4q2_0, O4q2_1, O5_0, O5_1, O6_0, O6_1, ..., O15_0, O15_1, O17_0, O17_1, O18_0, O18_1, alpha_1, ..., alpha_8)
  ! The coefficients of each operator has to be given in units GeV^(-2).
  ! _0 and _1 correspond to isoscalar and isovector coefficients. 
  ! Notice that O2 and O16 do not exist.
  ! Also, for O1 and O4, there also q^2 suppressed operators O1q2 and O4q2.
  ! The final 8 coefficients (alpha_1, ..., alpha_8) indicate which nuclear response functions must be calculated
  ! If alpha_i = 0, the corresponding nuclear response function will be skipped to save computation time
  ! Ideally, the user should not edit the alpha_i directly, but only via NRET_UpdateNRCoefficients.
  !
  ! 'NREFT_CPT': 
  ! (O1_p, O2_p, ..., O23_p, O100_p, O104_p, O1_n, O2_n, ..., O23_n, O100_n, O104_n, DM spin)
  ! The coefficients of each operator has to be given in units GeV^(-2).
  ! The definition of these operators follows DirectDM. For example, O13 = O6/(mpi^2+q^2).
  !
  ! Note that the two bases for NR effective operators are not equivalent.
  ! WIMP type 'NREffectiveTheory' contains a number of additional operators that are relevant only for spin-1 DM.
  ! WIMP type 'NREFT_CPT' includes a number of long-distance (dipole) operators and operators with meson poles.
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
  CHARACTER(LEN=1024) :: g_file = ''
  INTEGER :: Nvmin = -1
  REAL*8, ALLOCATABLE :: vmin(:)
  REAL*8, ALLOCATABLE :: g_vmin(:)
  REAL*8, ALLOCATABLE :: h_vmin(:)
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
  REAL*8, ALLOCATABLE  :: Jiso(:)  ! Calculated internally
  
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
  
  ! Number of bins with efficiencies.
  ! Will calculate rates for each bin plus total.
  INTEGER :: Nbins = -1
  
  ! Array of size [1:Niso,1:NE,0:Nbins] with the third index for the
  ! bin (zero for full range)
  REAL*8, ALLOCATABLE :: eff(:,:,:)  
  
  ! Form factors -------------------------------
  ! Tabulated form factors combined with prefactors.
  ! Arrays of size [0:9,1:4,1:NE,1:Niso].  
  ! The first index denotes the type of nuclear response function:
  !		alpha = 0:	 Helm form factor
  !		alpha = 1 ... 8: Nuclear response functions given by the tables in Wbar/.
  !				 Currently, these are the ones provided by Anand et. al. [1308.6288].
  !		alpha = 9:	 speficic response functions for SD scattering, given via the tables in SDFF/.
  !				 Currently, these are the ones provided by Klos et. al. [1304.7684].
  ! The second index denotes the isospin combination in the order
  ! (00, 01, 10, 11).
  ! The third index corresponds to the energies stored in E(:)
  ! The fourth index corresponds to the isotopes stored in Aiso(:) etc.
  ! NOTE: Need only be calculated once.
  REAL*8, ALLOCATABLE :: WTilde(:,:,:,:)
  
  ! Halo ---------------------------------------
  ! The minimum velocity for producing a recoil of energy E, given
  ! by vmin = sqrt{M E/(2\mu^2)} [km/s].
  ! Array of size [1:NE,1:Niso] that needs to be recalculated when the
  ! WIMP mass changes.
  REAL*8, ALLOCATABLE :: vmin(:,:)
  
  ! Tabulated mean inverse speed (g_vmin) [s/km] at the above vmin.
  REAL*8, ALLOCATABLE :: g_vmin(:,:)

  ! Tabulated mean speed (h_vmin) [km/s] at the above vmin.
  REAL*8, ALLOCATABLE :: h_vmin(:,:)

  ! structure for the differential rates,
  ! separately for each isotope of the experiment.
  ! These represent rates before efficiency cuts, but they already
  ! do include the mass fraction of the corresponding isotope.
  ! Units of this are cpd/kg/keV
  ! Array is of size [1:NE,1:Niso].
  REAL*8, ALLOCATABLE :: dRdEiso(:,:)
  
  ! Integrated rate --------
  ! Efficiency-corrected rates.  Array is of size
  ! [0:Nbins] with the index being that of the bin
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
