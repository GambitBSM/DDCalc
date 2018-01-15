MODULE DDDetectors

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! DDDetectors
!    Routines for initializing and modifying the detector setup,
!    notatably by setting isotopes to use, loading efficiencies from
!    file, and tabulating form factors.
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

USE DDTypes
USE DDCommandLine
USE DDNuclear
USE DDUtils
USE DDInput
USE DDOutput

IMPLICIT NONE

PUBLIC :: DDCalc_GetDetector,DDCalc_SetDetector
INTERFACE DDCalc_GetDetector
  MODULE PROCEDURE GetDetector
END INTERFACE
INTERFACE DDCalc_SetDetector
  MODULE PROCEDURE SetDetector
END INTERFACE

CONTAINS


!=======================================================================
! LOADING EFFICIENCY FILES
!=======================================================================

! ----------------------------------------------------------------------
! Loads and returns tabulated efficiencies from a file.
! 
! Required input argument:
!     file       The file to load efficiencies from.
!     Nbins      Number of intervals/bins expected.
!     NE         Number of recoil energy tabulation points expected.
! Required output arguments:
!     eff        Array to contain efficiency curves for
!                the total range and each interval/bin.  Allocated to
!                size [1:NE,0:Nbins].
!
SUBROUTINE LoadEfficiencyFile(file,Nbins,NE,eff)
  IMPLICIT NONE
  CHARACTER(LEN=*), INTENT(IN) :: file
  INTEGER, INTENT(IN) :: Nbins, NE
  REAL*8, INTENT(OUT) :: eff(:,0:)
  LOGICAL :: status
  INTEGER :: Kcol,Keff,Nrow,Ncol,Nvalid,Neff
  REAL*8, ALLOCATABLE :: data(:,:)
  
  ! Load table from file
  CALL LoadTable(file=file,Nrow=Nrow,Ncol=Ncol,data=data,status=status)
  IF ((.NOT. status) .OR. (Ncol .LT. 2)) THEN
    WRITE(0,*) 'ERROR: Failed to load data from file ' // TRIM(file) // '.'
    STOP
  END IF

  IF ( Nrow .NE. NE ) THEN
    WRITE(0,*) 'ERROR: Number of entries in file ' // TRIM(file) // ' does not agree with expectation.'
    STOP
  END IF
  
  ! Find number of valid efficiency columns
  Nvalid = 0
  DO Kcol = 1,Ncol
    IF (ALL(data(:,Kcol) .GE. 0d0) .AND. ALL(data(:,Kcol) .LE. 1.00001d0)) Nvalid = Nvalid + 1
  END DO
  IF (Nvalid .LE. 0) THEN
    WRITE(0,*) 'ERROR: Failed to find valid data in file ' // TRIM(file) // '.'
    STOP
  END IF
  
  ! Now get efficiencies
  Neff = Nvalid - 1

  IF ( Neff .EQ. Nbins ) THEN
    Keff = 0
    DO Kcol = 1,Ncol
      IF (ALL(data(:,Kcol) .GE. 0d0) .AND. ALL(data(:,Kcol) .LE. 1.00001d0) &
          .AND. (Keff .LE. Neff)) THEN
          eff(:,Keff) = data(:,Kcol)
        Keff = Keff + 1
      END IF
    END DO
  ELSE
    WRITE(0,*) 'ERROR: Number of entries in file ' // TRIM(file) // ' does not agree with expectation.'
    STOP
  END IF

END SUBROUTINE

SUBROUTINE LoadEnergyFile(file,NE,E)
  IMPLICIT NONE
  CHARACTER(LEN=*), INTENT(IN) :: file
  INTEGER, INTENT(OUT) :: NE
  REAL*8, ALLOCATABLE, INTENT(OUT) :: E(:)
  LOGICAL :: status
  INTEGER :: Kcol,Keff,Nrow,Ncol,Nvalid,Neff
  REAL*8, ALLOCATABLE :: data(:,:)
  
  ! Load table from file
  CALL LoadTable(file=file,Nrow=Nrow,Ncol=Ncol,data=data,status=status)
  IF ((.NOT. status) .OR. (Ncol .NE. 1)) THEN
    WRITE(0,*) 'ERROR: Failed to load data from file ' // TRIM(file) // '.'
    STOP
  END IF

  NE = Nrow

  ALLOCATE(E(NE))
  E(:) = data(:,1)
  
END SUBROUTINE


!=======================================================================
! DETECTOR SETUP ROUTINES
!=======================================================================

! ----------------------------------------------------------------------
! Get various detector quantities.
! Note this only allows access to fixed quantities, not quantities that
! are recalculated for each WIMP (see GetRates() for that).
! 
! Optional input argument:
!   D           The DetectorStruct structure to extract detector
!               quantities from.  If not given, a default structure
!               (internally stored) will be used.
! Optional output arguments:
!   mass        Detector fiducial mass [kg]
!   time        Detector exposure time [day]
!   exposure    Detector exposure [kg day]
!   Nevents     Allocatable integer array containing the number of observed 
!               events. Allocated to size [0:Nbins]
!   Backgr      Allocatable integer array containing the average expected 
!               background events. Allocated to size [0:Nbins]
!   Niso        Number of isotopes
!   Ziso        Allocatable integer array to be filled with isotope
!               atomic numbers. Allocated to size [1:Niso].
!   Aiso        Allocatable integer array to be filled with isotope
!               atomic mass numbers. Allocated to size [1:Niso].
!   fiso        Allocatable real array to be filled with isotope
!               mass fractions. Allocated to size [1:Niso].
!   Miso        Allocatable real array to be filled with isotope
!               nuclear masses [GeV]. Allocated to size [1:Niso].
!   NE          Number of tabulated recoil energies E
!   E           Allocatable real array to be filled with recoil energies
!               [keV]. Allocated to size [1:NE].
!   Nbins       Number of subintervals (0 for only total interval)
!   eff         Allocatable dimension=3 array containing efficiencies
!               as a function of recoil energy. Allocated to size
!               [1:Niso,1:NE,0:Nbins], where the first index is over isotopes
!               the second index is over recoil energies and the third
!               index is over the sub-interval number (0 for the total interval).
!   Wsi,Wsd     Allocatable dimension=3 array containing weighted form
!               factors for spin-independent (SI) and spin-dependent
!               (SD) couplings.  Allocated to size [-1:1,1:NE,1:Niso].
!   intervals   LOGICAL indicating if rates for intervals/bins are
!               to be calculated (used for max gap).
! 
SUBROUTINE GetDetector(D,mass,time,exposure,Nevents,background,         &
                       Niso,Ziso,Aiso,fiso,Miso,NE,E,                   &
                       Nbins,eff,Wsi,Wsd,intervals)
  IMPLICIT NONE
  TYPE(DetectorStruct), INTENT(IN) :: D
  LOGICAL, INTENT(OUT), OPTIONAL :: intervals
  INTEGER, INTENT(OUT), OPTIONAL :: Nevents(:),Niso,NE,Nbins
  INTEGER, ALLOCATABLE, INTENT(OUT), OPTIONAL :: Ziso(:),Aiso(:)
  REAL*8, INTENT(OUT), OPTIONAL :: mass,time,exposure,background(:)
  REAL*8, ALLOCATABLE, INTENT(OUT), OPTIONAL :: fiso(:),Miso(:),E(:),   &
          eff(:,:,:),Wsi(:,:,:),Wsd(:,:,:)
  
  IF ( .NOT. D%InitSuccess ) THEN
    WRITE(*,*) 'ERROR: Cannot get information from a detector that has not been correctly initialized.'
    STOP      
  END IF

  ! Exposures
  IF (PRESENT(mass))     mass     = D%mass
  IF (PRESENT(time))     time     = D%time
  IF (PRESENT(exposure)) exposure = D%exposure
  
  ! Observed events and expected background events
  IF (PRESENT(Nevents))    Nevents    = D%Nevents
  IF (PRESENT(background)) background = D%Backgr
  
  ! Isotope data
  IF (PRESENT(Niso)) Niso = D%Niso
  IF (PRESENT(Ziso)) THEN
    ALLOCATE(Ziso(D%Niso))
    Ziso = D%Ziso
  END IF
  IF (PRESENT(Aiso)) THEN
    ALLOCATE(Aiso(D%Niso))
    Aiso = D%Aiso
  END IF
  IF (PRESENT(fiso)) THEN
    ALLOCATE(fiso(D%Niso))
    fiso = D%fiso
  END IF
  IF (PRESENT(Miso)) THEN
    ALLOCATE(Miso(D%Niso))
    Miso = D%Miso
  END IF
  
  ! Recoil energies
  IF (PRESENT(NE)) NE = D%NE
  IF (PRESENT(E)) THEN
    ALLOCATE(E(D%NE))
    E = D%E
  END IF
  
  ! Efficiencies
  IF (PRESENT(Nbins))     Nbins     = D%Nbins
  IF (PRESENT(eff)) THEN
    ALLOCATE(eff(D%Niso,D%NE,0:D%Nbins))
    eff = D%eff
  END IF
  
  ! Weighted form factors
  IF (PRESENT(Wsi)) THEN
    ALLOCATE(Wsi(-1:1,D%NE,0:D%Niso))
    Wsi = D%Wsi
  END IF
  IF (PRESENT(Wsd)) THEN
    ALLOCATE(Wsd(-1:1,D%NE,0:D%Niso))
    Wsd = D%Wsd
  END IF
  
  ! Calculate rates for intervals/bins?
  IF (PRESENT(intervals)) intervals = D%intervals
  
END SUBROUTINE


! ----------------------------------------------------------------------
! Set various detector quantities.  Making changes will reset all
! per-WIMP calculations like rates.
! Note this only allows access to fixed quantities, or quantities that
! need be calculated only once, not quantities that are recalculated
! for each WIMP (see SetRates() for that).
! 
! Optional input/output argument:
!   D           The DetectorStruct structure containing detector
!               data to be modified.  If not given, a default structure
!               (internally stored) will be used.
!

! The following steps are required to define a new detector
!
! 1) Definition of the recoil energy grid used for calculating dR/dE
!   NE          Number of tabulated recoil energies E
!   E           Array of size [1:NE] containing recoil energies [keV].
!
! 2) Definition of bins
!   Nbins       Number of sub-intervals/bins in data (0 if total only).
!
! 3) Specification of the observed event numbers
! EITHER
!   Nevents_tot Total number of observed events.
! OR
!   Nevents_bin Array of size [1:Nbins] containing the number of
!               observed events for each interval.
!
! The statistical method to be used in the analysis depends on the information
! on the observed events supplied by the user. There are three possibilities
! 3a) If Nevents_tot is provided and non-negative, Total Poisson will be used.
! 3b) If Nevents_tot is provided and negative, Maximum Gap will be used and the
!     total number of observed events will be set to Nbins-1. Nbins = 0 results
!     in an error message
! 3c) If Nevents_bin is provided, Binned Poisson will be used.
! Note that it is possible to specify Nevents_tot even if Nbins > 0, in which
! case the information from the individual bins will be ignored.
!
! 4) OPTIONAL: Specification of expected backgrounds.
! EITHER
!   Backgr_tot  Total number of expected background events
! OR
!   Backgr_bin  Array of size [1:Nbins] containing the expected background for
!   		each interval.
!
! Note that information on expected backgrounds will be ignored for Maximum Gap.
! For Total Poisson or Binned Poisson, if no background expectation is specified
! (or if it is specified to be zero), this is taken to mean that there is no
! background model and hence the background expectation is taken as a nuisance
! parameter.
!
! 5) Target elements
! EITHER
!   Niso        Number of isotopes
!   Ziso        Integer array of size [1:Niso] containing atomic
!               numbers.
!   Aiso        Integer array of size [1:Niso] containing atomic
!               mass numbers.
!   fiso        Array of size [1:Niso] containing isotope mass
!               fractions.
! OR
!   Nelem       Number of elements in compound.
!   Zelem       Integer array of size [1:Nelem] containing atomic
!               numbers of compound elements.
!   stoich      Integer array of size [1:Nelem] containing compound
!               stoichiometry.  For example, CF3Cl would be {1,3,1}.
!               Default is 1 for each element.
!
! Note that the target composition can no longer be changed after a detector
! has been successfully initialized.
!
! 6) Detector exposure
! EITHER
!   mass        Detector fiducial mass [kg]
!   time        Detector run-time [day]
! OR
!   exposure    Detector exposure [kg day]
!
! 7) Detector efficiencies
! EITHER
!   eff_all     Array of size [1:NE,0:Nbins] containing efficiencies
!               as a function of recoil energy, where the first index
!               is over recoil energies and the second index is over
!               the sub-interval number (0 for the total interval).
!               The same efficiency will be used for all target elements
!               and isotopes.
! OR
!   eff         Array of size [1:Niso,1:NE,0:Nbins] containing efficiencies
!               as a function of recoil energy as above, but specified
!               separately for every target isotope.
!               If eff is provided together with Nelem and Zelem, it is
!               taken to have the size [1:Nelem,1:NE,0:Nbins] and the same
!               efficiency will be used for all isotopes of each element
! OR
!  eff_file_all File from which efficiencies shoud be read. The efficiencies
!		must have the same format as eff_all and will be interpreted
!		in the same way, i.e. the same efficiency will be used for
!		all target elements and isotopes.
! OR
!  eff_file	Array of files from which efficiencies shoud be read. The efficiencies
!		must have the same format as eff_all and will be interpreted
!		in the same way, i.e. eff_file should have siye [1:Niso] or
!		[1:Nelme] depending on whether Nelem and Zelem have been specified.
!
! 8) OPTIONAL: Energy threshold
!   Emin        If given, sets all efficiencies below the given energy
!               [keV] to zero, removing all contributions from recoils
!               at lower energies.
!

SUBROUTINE SetDetector(D,mass,time,exposure,Nbins,                      &
                       Nevents_tot,Nevents_bin,Backgr_tot,Backgr_bin,   &
                       Niso,Ziso,Aiso,fiso,Nelem,Zelem,stoich,          &
                       NE,E,E_file,eff_file,eff_file_all,eff,eff_all,   &
                       intervals,Emin)
  IMPLICIT NONE
  TYPE(DetectorStruct), INTENT(INOUT) :: D
  CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: eff_file(:)
  CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: eff_file_all,E_file
  LOGICAL, INTENT(IN), OPTIONAL :: intervals
  INTEGER, INTENT(IN), OPTIONAL :: Nevents_tot,Nevents_bin(:),Niso,Nelem,NE,Nbins
  INTEGER, INTENT(IN), OPTIONAL :: Ziso(:),Aiso(:),Zelem(:),stoich(:)
  REAL*8, INTENT(IN), OPTIONAL :: mass,time,exposure,Backgr_tot,Backgr_bin(:),Emin
  REAL*8, INTENT(IN), OPTIONAL :: fiso(:),E(:),eff(:,:,0:),eff_all(:,0:)
  LOGICAL :: E_change,eff_change,intervals_change
  INTEGER :: KE,Kiso,Neff,ind_elem,ind_iso,ind,Niso_temp,eff_files,file_number
  INTEGER, ALLOCATABLE :: stoich0(:)
  REAL*8, ALLOCATABLE :: eff_new(:,:,:),eff_all_new(:,:)

  ! Indicate if quantities changed, indicating need for
  ! array resizing and initialization
  E_change   = .FALSE.
  eff_change = .FALSE.
  intervals_change = .FALSE.

  ! Check whether any unallowed changes to the already initialized detector are attempted...
  IF (D%InitSuccess) THEN
    IF ( (PRESENT(Niso)) .OR. (PRESENT(Ziso)) .OR. (PRESENT(Aiso)) .OR. (PRESENT(fiso)) &
             .OR. (PRESENT(Nelem)) .OR. (PRESENT(Zelem)) .OR. (PRESENT(stoich)) ) THEN
      WRITE(*,*) 'ERROR: Cannot change target composition after detector has been successfully initialized.'
      STOP
    END IF
  END IF

  IF (PRESENT(E_file)) THEN
    IF (E_file .NE. '') THEN
      CALL LoadEnergyFile(E_file,D%NE,D%E)
      E_change = .TRUE.
    END IF
  ELSE IF ( (PRESENT(E)) .AND. (PRESENT(NE)) ) THEN
    IF (NE > 0) THEN
      IF (ALLOCATED(D%E) .AND. (NE .NE. D%NE)) DEALLOCATE(D%E)
      IF (.NOT. ALLOCATED(D%E)) ALLOCATE(D%E(NE))
      D%NE = NE
      D%E = E(1:NE)
      IF (ALLOCATED(D%E_cache)) DEALLOCATE(D%E_cache)
      ALLOCATE(D%E_cache(D%NE))
      D%E_cache = D%E
      D%InitSuccess = .False. ! If NE is changed, a new efficiency array must be provided
      E_change = .TRUE.
    END IF
  END IF

  IF (.NOT. D%NE > 0) THEN
    D%InitSuccess = .False.
    RETURN      
  END IF

  ! Change of Nbins
  ! Notice that if Nbins is changed, D%StatisticFlag is set to -1.
  ! Hence, in order for the detector to be correctly initalized, one also has to provide Nevents_tot or Nevents_bin
  IF ( PRESENT(Nbins) ) THEN
    IF (Nbins .NE. D%Nbins) THEN
      IF (ALLOCATED(D%Nevents)) DEALLOCATE(D%Nevents)
      IF (ALLOCATED(D%Backgr)) DEALLOCATE(D%Backgr)
      D%InitSuccess = .False. ! If Nbins is changed, a new efficiency array must be provided
      ALLOCATE(D%Nevents(0:Nbins),D%Backgr(0:Nbins))
      D%Backgr(0:) = 0
      D%StatisticFlag = -1
    END IF
    D%Nbins = Nbins
  END IF

  ! Read in observed events, which then also defines the StatisticFlag and the intervals flag
  IF ( PRESENT(Nevents_tot) ) THEN
    IF ( PRESENT(Nevents_bin) ) THEN
      WRITE (*,*) 'Error: both Nevents_tot and Nevents_bin are given as arguments.'
      D%StatisticFlag = -1
      D%InitSuccess = .False.
      RETURN      
    END IF
    IF (Nevents_tot .LT. 0) THEN
      IF (D%Nbins .LT. 1) THEN
        WRITE (*,*) 'Error: maximum gap method can only be used for non-zero number of bins'
        D%InitSuccess = .False.
        RETURN      
      END IF
      D%StatisticFlag = 2 ! MaxGap
      D%intervals = .true.
      D%Nevents(0) = D%Nbins - 1
      intervals_change = .true.
    ELSE
      D%StatisticFlag = 0 ! TotalPoisson
      D%intervals = .false.
      D%Nevents(0) = Nevents_tot
      intervals_change = .true.
    END IF
    ! info about binned observed events is irrelevant for MaxGap and TotalPoisson, so just set it to zero
    IF (D%Nbins .GT. 0) THEN
      D%Nevents(1:) = 0
    END IF
  END IF
  IF ( PRESENT(Nevents_bin) ) THEN
    IF (D%Nbins .LT. 1) THEN
      D%InitSuccess = .False.
      RETURN      
    END IF
    IF (D%Nbins .NE. SIZE(Nevents_bin)) THEN
      WRITE (*,*) 'Error: wrong size of Nevents_bin'
      D%InitSuccess = .False.
      RETURN     
    END IF
    D%StatisticFlag = 1 ! BinnedPoisson
    D%intervals = .true.
    D%Nevents(0) = SUM(Nevents_bin(:))
    D%Nevents(1:) = Nevents_bin(:)
    intervals_change = .true.
  END IF

  ! Read in background events
  IF ( PRESENT(Backgr_tot) ) THEN
    IF (PRESENT(Backgr_bin)) THEN
      WRITE (*,*) 'Error: both Backgr_tot and Backgr_bin are given as arguments.'
      D%StatisticFlag = -1
      D%InitSuccess = .False.
      RETURN  
    END IF
    D%Backgr(0) = Backgr_tot
    IF (D%Nbins .GT. 0) THEN
      D%Backgr(1:) = 0
    END IF
  END IF
  IF (PRESENT(Backgr_bin)) THEN
    IF (D%Nbins .LT. 1) THEN
      D%InitSuccess = .False.
      RETURN      
    END IF
    IF (D%Nbins .NE. SIZE(Backgr_bin)) THEN
      WRITE (*,*) 'Error: wrong size of Backgr_bin'
      D%InitSuccess = .False.
      RETURN     
    END IF
    D%Backgr(0) = SUM(Backgr_bin(:))
    D%Backgr(1:) = Backgr_bin(:)
  END IF


  ! Check whether D%StatisticFlag has a valid value
  IF (D%StatisticFlag .LT. 0) THEN
    D%InitSuccess = .False.
    RETURN 
  END IF


  IF ( D%Nbins < 0) THEN
    D%InitSuccess = .False.
    RETURN      
  END IF

  ! Isotope data
  IF (PRESENT(Niso)) THEN
    IF (Niso .EQ. D%Niso) THEN
      IF (PRESENT(Ziso)) THEN
        D%Ziso = Ziso(1:Niso)
      END IF
      IF (PRESENT(Aiso)) THEN
        D%Aiso = Aiso(1:Niso)
      END IF
      IF (PRESENT(fiso)) THEN
        D%fiso = fiso(1:Niso)
      END IF
      D%Miso = IsotopeMass(D%Ziso,D%Aiso)
    ELSE IF (PRESENT(Ziso) .AND. PRESENT(Aiso) .AND. PRESENT(fiso)) THEN
      IF (ALLOCATED(D%Ziso)) DEALLOCATE(D%Ziso)
      IF (ALLOCATED(D%Aiso)) DEALLOCATE(D%Aiso)
      IF (ALLOCATED(D%fiso)) DEALLOCATE(D%fiso)
      IF (ALLOCATED(D%Miso)) DEALLOCATE(D%Miso)
      ALLOCATE(D%Ziso(Niso),D%Aiso(Niso),D%fiso(Niso),D%Miso(Niso))
      D%Niso = Niso
      D%Ziso = Ziso(1:Niso)
      D%Aiso = Aiso(1:Niso)
      D%fiso = fiso(1:Niso)
      D%Miso = IsotopeMass(D%Ziso,D%Aiso)
      D%InitSuccess = .False. ! If Niso is changed, a new efficiency array must be provided
    END IF
  END IF

  IF (PRESENT(Nelem) .AND. PRESENT(Zelem)) THEN
    ALLOCATE(stoich0(Nelem))
    IF (PRESENT(stoich)) THEN
      stoich0 = stoich
    ELSE
      stoich0 = 1
    END IF
    CALL CompoundIsotopeList(Nelem,Zelem,stoich0,                       &
                             D%Niso,D%Ziso,D%Aiso,D%fiso,D%Miso)
    D%InitSuccess = .False. ! If Niso is changed, a new efficiency array must be provided
  END IF
  
  IF (.NOT. D%Niso > 0) THEN
    D%InitSuccess = .False.
    RETURN      
  END IF

  ! Exposures
  IF (PRESENT(mass)) THEN
    D%mass     = mass
    D%exposure = D%mass * D%time
  END IF
  IF (PRESENT(time)) THEN
    D%time     = time
    D%exposure = D%mass * D%time
  END IF
  IF (PRESENT(exposure)) THEN
    IF (PRESENT(mass) .OR. PRESENT(time)) THEN
      D%InitSuccess = .False.
      RETURN           
    END IF
    D%mass     = -1d0
    D%time     = -1d0
    D%exposure = exposure
  END IF

  IF (.NOT. D%exposure > 0d0) THEN
    D%InitSuccess = .False.
    RETURN      
  END IF

  ! Set efficiencies
  ! ...from file (note that right now only one efficiency can be provided for all isotopes)
  IF (PRESENT(eff_file_all)) THEN
    IF (eff_file_all .NE. '') THEN
      ALLOCATE(eff_all_new(D%NE,0:D%Nbins))
      CALL LoadEfficiencyFile(eff_file_all,D%Nbins,D%NE,eff_all_new(:,0:))
      eff_change = .TRUE.
    END IF
  ELSE IF (PRESENT(eff_all)) THEN
    ALLOCATE(eff_all_new(SIZE(eff_all,1),0:(SIZE(eff_all,2)-1)))
    eff_all_new(:,0:) = eff_all(:,0:)
    eff_change = .TRUE.
  END IF
  IF (eff_change) THEN
    ! check whether eff_all has the correct dimensions
    IF ((SIZE(eff_all_new,1) .NE. D%NE) .OR. (SIZE(eff_all_new,2) .NE. (D%Nbins+1))) THEN
      WRITE (*,*) 'Wrong dimension of eff_all'
      D%InitSuccess = .False.
      RETURN  
    END IF
    IF (ALLOCATED(D%eff))  DEALLOCATE(D%eff)
    ALLOCATE(D%eff(D%Niso,D%NE,0:D%Nbins))
    DO ind = 1,D%Niso
      D%eff(ind,1:D%NE,0:D%Nbins) = eff_all_new(1:D%NE,0:D%Nbins)
    END DO
    D%InitSuccess = .TRUE.
  ELSE
    IF (PRESENT(eff_file)) THEN
      eff_files = SIZE(eff_file)
      IF (PRESENT(Nelem) .AND. PRESENT(Zelem)) THEN
        IF(eff_files .NE. Nelem) THEN
          WRITE (*,*) 'Wrong number of efficiency files'
          D%InitSuccess = .False.
          RETURN  
        END IF
      ELSE
        IF(eff_files .NE. D%Niso) THEN
          WRITE (*,*) 'Wrong number of efficiency files'
          D%InitSuccess = .False.
          RETURN  
        END IF
      END IF

      ALLOCATE(eff_new(eff_files,D%NE,0:D%Nbins))

      DO file_number = 1,eff_files
        IF (eff_file(file_number) .NE. '') THEN
          CALL LoadEfficiencyFile(eff_file(file_number),D%Nbins,D%NE,eff_new(file_number,:,0:))
          eff_change = .TRUE.
        END IF
      END DO
    ELSE IF (PRESENT(eff)) THEN
      ALLOCATE(eff_new(SIZE(eff,1),SIZE(eff,2),0:(SIZE(eff,3)-1)))
      eff_new(:,:,0:)=eff(:,:,0:)
      eff_change = .TRUE.
    END IF
    IF (eff_change) THEN
      ! check whether eff has the correct dimensions
      IF ((SIZE(eff_new,2) .NE. D%NE) .OR. (SIZE(eff_new,3) .NE. (D%Nbins+1))) THEN
        WRITE (*,*) 'Wrong dimension of eff'
        D%InitSuccess = .False.
        RETURN  
      END IF
      IF (ALLOCATED(D%eff))  DEALLOCATE(D%eff)
      ALLOCATE(D%eff(D%Niso,D%NE,0:D%Nbins))
      IF (PRESENT(Nelem) .AND. PRESENT(Zelem)) THEN
        IF (SIZE(eff_new,1) .NE. Nelem) THEN
          WRITE (*,*) 'Wrong dimension of eff'
          D%InitSuccess = .False.
          RETURN  
        END IF
        ind = 1
        DO ind_elem = 1,Nelem
          CALL GetNiso(Zelem(ind_elem),Niso_temp) ! this assigns Niso_temp
          DO ind_iso = 1,Niso_temp
            D%eff(ind,1:D%NE,0:D%Nbins) = eff_new(ind_elem,1:D%NE,0:D%Nbins)
            ind = ind + 1
          END DO
        END DO
      ELSE
        IF (SIZE(eff_new,1) .NE. D%Niso) THEN
          WRITE (*,*) 'Wrong dimension of eff'
          D%InitSuccess = .False.
          RETURN
        END IF
        D%eff   = eff_new(1:D%Niso,1:D%NE,0:D%Nbins)
      END IF
      eff_change = .TRUE.
      D%InitSuccess = .TRUE.
    END IF
  END IF

  IF (.NOT. D%InitSuccess) THEN
    RETURN
  END IF

  
  ! Apply threshold cut.
  ! We move all E < Emin tabulation points to Emin.
  IF (PRESENT(Emin)) THEN
    ! First reset to original tabulation
    D%E = MAX(D%E_cache,Emin)
    E_change = .TRUE.
  END IF
  
  ! Calculate weighted form factors (SI)
  IF (ALLOCATED(D%Wsi)) DEALLOCATE(D%Wsi)
  ALLOCATE(D%Wsi(-1:1,D%NE,D%Niso))
  DO Kiso = 1,D%Niso
    CALL CalcWSI(D%Ziso(Kiso),D%Aiso(Kiso),D%NE,                   &
                 EToQ(D%E,D%Miso(Kiso)),D%Wsi(:,:,Kiso))
  END DO

  ! Calculate weighted form factors (SD)
  IF (ALLOCATED(D%Wsd)) DEALLOCATE(D%Wsd)
  ALLOCATE(D%Wsd(-1:1,D%NE,D%Niso))
  DO Kiso = 1,D%Niso
    CALL CalcWSD(D%Ziso(Kiso),D%Aiso(Kiso),D%NE,                   &
                 EToQ(D%E,D%Miso(Kiso)),D%Wsd(:,:,Kiso))
  END DO
  
  ! Resize halo velocity arrays if necessary
  IF (E_change) THEN
    IF (ALLOCATED(D%vmin)) DEALLOCATE(D%vmin)
    ALLOCATE(D%vmin(D%NE,D%Niso))
    IF (ALLOCATED(D%eta))  DEALLOCATE(D%eta)
    ALLOCATE(D%eta(D%NE,D%Niso))
  END IF
  
  ! Number of intervals/bins to do calculations for.
  ! Used for array sizing below.
  IF (D%intervals) THEN
    Neff = D%Nbins
  ELSE
    Neff = 0
  END IF
  
  ! Resize dRdEiso if necessary
  IF (E_change) THEN
    IF (ALLOCATED(D%dRdEiso)) DEALLOCATE(D%dRdEiso)
    ALLOCATE(D%dRdEiso(D%NE,D%Niso))
  END IF

  ! Resize rate arrays if necessary
  IF (E_change .OR. eff_change .OR. intervals_change) THEN
    IF (ALLOCATED(D%R)) DEALLOCATE(D%R)
    ALLOCATE(D%R(0:Neff))
  END IF
  
  ! Resize event arrays if necessary
  IF (eff_change .OR. intervals_change) THEN
    IF (ALLOCATED(D%MuSignal)) DEALLOCATE(D%MuSignal)
    ALLOCATE(D%MuSignal(0:Neff))
  END IF
  
  ! Set all calculable quantities to zero
  IF (E_change .OR. eff_change .OR. intervals_change) THEN
    D%vmin        = 0d0
    D%eta         = 0d0
    D%dRdEiso     = 0d0
    D%R           = 0d0
    D%MuSignal    = 0d0
  END IF

  
END SUBROUTINE
  
END MODULE
