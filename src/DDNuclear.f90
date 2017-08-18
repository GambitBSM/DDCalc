MODULE DDNuclear

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! DDNuclear
!    Routines for determining isotope and nuclear properties, including
!    isotopic compositions, nuclear masses, and form factors (both
!    spin-independent and spin-dependent).
!    [NOTE: LIMITED IMPLEMENTATION FOR SD CASE -- MOSTLY SET TO ZERO.]
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

USE DDConstants

IMPLICIT NONE
PRIVATE

PUBLIC :: IsotopeMass, EtoQ, ElementIsotopeList, GetNiso, CompoundIsotopeList, CalcWSD, CalcWSI

CONTAINS


! ----------------------------------------------------------------------
! For the given element, returns number of isotopes
! 
! Input argument:
!     Z          Atomic number of element
! Output arguments:
!     Niso       Number of isotopes
!
! [added by Sebastian Wild, August 2017] 
! 
SUBROUTINE GetNiso(Z,Niso)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: Z
  INTEGER, INTENT(OUT) :: Niso
  INTEGER, PARAMETER :: NELEMENTS = 92 ! when changing this, change also in ElementIsotopeList
  ! Number of stable isotopes for given element (indexed by Z)
  INTEGER, PARAMETER :: ELEMENT_NISO(NELEMENTS) =                       &
    (/  2,  2,  2,  1,  2,  2,  2,  3,  1,  3,  1,  3,  1,  3,  1,  4,  &
        2,  3,  3,  6,  1,  5,  2,  4,  1,  4,  1,  5,  2,  5,  2,  5,  &
        1,  6,  2,  6,  2,  4,  1,  5,  1,  7,  0,  7,  1,  6,  2,  8,  &
        2, 10,  2,  8,  1,  9,  1,  7,  2,  4,  1,  7,  0,  7,  2,  7,  &
        1,  7,  1,  6,  1,  7,  2,  6,  1,  5,  2,  7,  2,  6,  1,  7,  &
        2,  4,  1,  0,  0,  0,  0,  0,  0,  1,  0,  3 /)

  IF (Z .LE. NELEMENTS) THEN
    Niso = ELEMENT_NISO(Z)
  ELSE
    Niso = 0
  END IF

END SUBROUTINE


! ----------------------------------------------------------------------
! For the given element, fills in allocatable arrays containing isotopic
! atomic numbers (Z), atomic masses (A), mass fractions (f), and
! masses (M).
! 
! Input argument:
!     Z          Atomic number of element
! Output arguments:
!     Niso       Number of isotopes
!     Ziso       Allocatable array (integer) of isotopes' atomic numbers (Z)
!     Aiso       Allocatable array (integer) of isotopes' atomic masses (A)
!     fiso       Allocatable array of isotopes' mass fractions
!     Miso       Allocatable array of isotopes' nuclear masses [GeV]
! 
SUBROUTINE ElementIsotopeList(Z,Niso,Ziso,Aiso,fiso,Miso)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: Z
  INTEGER, INTENT(OUT) :: Niso
  INTEGER, ALLOCATABLE, INTENT(OUT) :: Ziso(:),Aiso(:)
  REAL*8, ALLOCATABLE, INTENT(OUT) :: fiso(:),Miso(:)
  INTEGER :: I1,I2
  ! Isotope data for all elements up to Z=92 (Uranium).
  ! Data is combined into single arrays; the ELEMENT_INDEX indicates
  ! where data for a particular element Z begins in those arrays.
  INTEGER, PARAMETER :: NELEMENTS = 92  ! when changing this, change also in GetNiso
  INTEGER, PARAMETER :: NISOTOPES = 286

  ! Number of stable isotopes for given element (indexed by Z)
  ! ELEMENT_NISO(NELEMENTS) is now defined in GetNiso

  ! First data array index for given element (indexed by Z)
  INTEGER, PARAMETER :: ELEMENT_INDEX(NELEMENTS) =                      &
    (/  1,  3,  5,  7,  8, 10, 12, 14, 17, 18, 21, 22, 25, 26, 29, 30,  &
       34, 36, 39, 42, 48, 49, 54, 56, 60, 61, 65, 66, 71, 73, 78, 80,  &
       85, 86, 92, 94,100,102,106,107,112,113,120,120,127,128,134,136,  &
      144,146,156,158,166,167,176,177,184,186,190,191,198,198,205,207,  &
      214,215,222,223,229,230,237,239,245,246,251,253,260,262,268,269,  &
      276,278,282,283,283,283,283,283,283,283,284,284 /)
  ! Atomic number for an isotope
  INTEGER, PARAMETER :: ISOTOPE_Z(NISOTOPES) =                          &
    (/  1,  1,  2,  2,  3,  3,  4,  5,  5,  6,  6,  7,  7,  8,  8,  8,  &
        9, 10, 10, 10, 11, 12, 12, 12, 13, 14, 14, 14, 15, 16, 16, 16,  &
       16, 17, 17, 18, 18, 18, 19, 19, 19, 20, 20, 20, 20, 20, 20, 21,  &
       22, 22, 22, 22, 22, 23, 23, 24, 24, 24, 24, 25, 26, 26, 26, 26,  &
       27, 28, 28, 28, 28, 28, 29, 29, 30, 30, 30, 30, 30, 31, 31, 32,  &
       32, 32, 32, 32, 33, 34, 34, 34, 34, 34, 34, 35, 35, 36, 36, 36,  &
       36, 36, 36, 37, 37, 38, 38, 38, 38, 39, 40, 40, 40, 40, 40, 41,  &
       42, 42, 42, 42, 42, 42, 42, 44, 44, 44, 44, 44, 44, 44, 45, 46,  &
       46, 46, 46, 46, 46, 47, 47, 48, 48, 48, 48, 48, 48, 48, 48, 49,  &
       49, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 51, 51, 52, 52, 52,  &
       52, 52, 52, 52, 52, 53, 54, 54, 54, 54, 54, 54, 54, 54, 54, 55,  &
       56, 56, 56, 56, 56, 56, 56, 57, 57, 58, 58, 58, 58, 59, 60, 60,  &
       60, 60, 60, 60, 60, 62, 62, 62, 62, 62, 62, 62, 63, 63, 64, 64,  &
       64, 64, 64, 64, 64, 65, 66, 66, 66, 66, 66, 66, 66, 67, 68, 68,  &
       68, 68, 68, 68, 69, 70, 70, 70, 70, 70, 70, 70, 71, 71, 72, 72,  &
       72, 72, 72, 72, 73, 74, 74, 74, 74, 74, 75, 75, 76, 76, 76, 76,  &
       76, 76, 76, 77, 77, 78, 78, 78, 78, 78, 78, 79, 80, 80, 80, 80,  &
       80, 80, 80, 81, 81, 82, 82, 82, 82, 83, 90, 92, 92, 92 /)
  ! Atomic mass number for an isotope
  INTEGER, PARAMETER :: ISOTOPE_A(NISOTOPES) =                          &
    (/  1,  2,  3,  4,  6,  7,  9, 10, 11, 12, 13, 14, 15, 16, 17, 18,  &
       19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34,  &
       36, 35, 37, 36, 38, 40, 39, 40, 41, 40, 42, 43, 44, 46, 48, 45,  &
       46, 47, 48, 49, 50, 50, 51, 50, 52, 53, 54, 55, 54, 56, 57, 58,  &
       59, 58, 60, 61, 62, 64, 63, 65, 64, 66, 67, 68, 70, 69, 71, 70,  &
       72, 73, 74, 76, 75, 74, 76, 77, 78, 80, 82, 79, 81, 78, 80, 82,  &
       83, 84, 86, 85, 87, 84, 86, 87, 88, 89, 90, 91, 92, 94, 96, 93,  &
       92, 94, 95, 96, 97, 98,100, 96, 98, 99,100,101,102,104,103,102,  &
      104,105,106,108,110,107,109,106,108,110,111,112,113,114,116,113,  &
      115,112,114,115,116,117,118,119,120,122,124,121,123,120,122,123,  &
      124,125,126,128,130,127,124,126,128,129,130,131,132,134,136,133,  &
      130,132,134,135,136,137,138,138,139,136,138,140,142,141,142,143,  &
      144,145,146,148,150,144,147,148,149,150,152,154,151,153,152,154,  &
      155,156,157,158,160,159,156,158,160,161,162,163,164,165,162,164,  &
      166,167,168,170,169,168,170,171,172,173,174,176,175,176,174,176,  &
      177,178,179,180,181,180,182,183,184,186,185,187,184,186,187,188,  &
      189,190,192,191,193,190,192,194,195,196,198,197,196,198,199,200,  &
      201,202,204,203,205,204,206,207,208,209,232,234,235,238 /)
  ! Atomic spin for an isotope
  REAL*8, PARAMETER :: ISOTOPE_J(NISOTOPES) =                           &
    (/0.5d0,1.0d0,0.5d0,0.0d0,1.0d0,1.5d0,1.5d0,3.0d0,1.5d0,0.0d0,      &
      0.5d0,1.0d0,0.5d0,0.0d0,2.5d0,0.0d0,0.5d0,0.0d0,1.5d0,0.0d0,      &
      1.5d0,0.0d0,2.5d0,0.0d0,2.5d0,0.0d0,0.5d0,0.0d0,0.5d0,0.0d0,      &
      1.5d0,0.0d0,0.0d0,1.5d0,1.5d0,0.0d0,0.0d0,0.0d0,1.5d0,4.0d0,      &
      1.5d0,0.0d0,0.0d0,3.5d0,0.0d0,0.0d0,0.0d0,3.5d0,0.0d0,2.5d0,      &
      0.0d0,3.5d0,0.0d0,6.0d0,3.5d0,0.0d0,0.0d0,1.5d0,0.0d0,2.5d0,      &
      0.0d0,0.0d0,0.5d0,0.0d0,3.5d0,0.0d0,0.0d0,1.5d0,0.0d0,0.0d0,      &
      1.5d0,1.5d0,0.0d0,0.0d0,2.5d0,0.0d0,0.0d0,1.5d0,1.5d0,0.0d0,      &
      0.0d0,4.5d0,0.0d0,0.0d0,1.5d0,0.0d0,0.0d0,0.5d0,0.0d0,0.0d0,      &
      0.0d0,1.5d0,1.5d0,0.0d0,0.0d0,0.0d0,4.5d0,0.0d0,0.0d0,2.5d0,      &
      1.5d0,0.0d0,0.0d0,4.5d0,0.0d0,0.5d0,0.0d0,2.5d0,0.0d0,0.0d0,      &
      0.0d0,4.5d0,0.0d0,0.0d0,2.5d0,0.0d0,2.5d0,0.0d0,0.0d0,0.0d0,      &
      0.0d0,2.5d0,0.0d0,2.5d0,0.0d0,0.0d0,0.5d0,0.0d0,0.0d0,2.5d0,      &
      0.0d0,0.0d0,0.0d0,0.5d0,0.5d0,0.0d0,0.0d0,0.0d0,0.5d0,0.0d0,      &
      0.5d0,0.0d0,0.0d0,4.5d0,4.5d0,0.0d0,0.0d0,0.5d0,0.0d0,0.5d0,      &
      0.0d0,0.5d0,0.0d0,0.0d0,0.0d0,2.5d0,3.5d0,0.0d0,0.0d0,0.5d0,      &
      0.0d0,0.5d0,0.0d0,0.0d0,0.0d0,2.5d0,0.0d0,0.0d0,0.0d0,0.5d0,      &
      0.0d0,1.5d0,0.0d0,0.0d0,0.0d0,3.5d0,0.0d0,0.0d0,0.0d0,1.5d0,      &
      0.0d0,1.5d0,0.0d0,5.0d0,3.5d0,0.0d0,0.0d0,0.0d0,0.0d0,2.5d0,      &
      0.0d0,3.5d0,0.0d0,3.5d0,0.0d0,0.0d0,0.0d0,0.0d0,3.5d0,0.0d0,      &
      3.5d0,0.0d0,0.0d0,0.0d0,2.5d0,2.5d0,0.0d0,0.0d0,1.5d0,0.0d0,      &
      1.5d0,0.0d0,0.0d0,1.5d0,0.0d0,0.0d0,0.0d0,2.5d0,0.0d0,2.5d0,      &
      0.0d0,3.5d0,0.0d0,0.0d0,0.0d0,3.5d0,0.0d0,0.0d0,0.5d0,0.0d0,      &
      0.0d0,0.5d0,0.0d0,2.5d0,0.0d0,0.0d0,3.5d0,7.0d0,0.0d0,0.0d0,      &
      3.5d0,0.0d0,4.5d0,0.0d0,3.5d0,0.0d0,0.0d0,0.5d0,0.0d0,0.0d0,      &
      2.5d0,2.5d0,0.0d0,0.0d0,0.5d0,0.0d0,1.5d0,0.0d0,0.0d0,1.5d0,      &
      1.5d0,0.0d0,0.0d0,0.0d0,0.5d0,0.0d0,0.0d0,1.5d0,0.0d0,0.0d0,      &
      0.5d0,0.0d0,1.5d0,0.0d0,0.0d0,0.5d0,0.5d0,0.0d0,0.0d0,0.5d0,      &
      0.0d0,4.5d0,0.0d0,0.0d0,3.5d0,0.0d0 /)
  ! Elemental mass fraction for an isotope
  REAL*8, PARAMETER :: ISOTOPE_F(NISOTOPES) =                           &
    (/0.9997d0,   0.0002997d0,1.032d-6,   1.000d0,    0.06578d0,  0.9342d0,   1.000d0,    0.1834d0,    &
      0.8166d0,   0.9880d0,   0.01202d0,  0.9961d0,   0.003920d0, 0.9973d0,   0.0004037d0,0.002250d0,  &
      1.000d0,    0.8964d0,   0.002809d0, 0.1008d0,   1.000d0,    0.7795d0,   0.1028d0,   0.1177d0,    &
      1.000d0,    0.9187d0,   0.04832d0,  0.03295d0,  1.000d0,    0.9475d0,   0.007712d0, 0.04460d0,   &
      0.0002243d0,0.7474d0,   0.2526d0,   0.003030d0, 0.0006006d0,0.9964d0,   0.9294d0,   0.0001196d0, &
      0.07051d0,  0.9666d0,   0.006773d0, 0.001447d0, 0.02292d0,  4.586d-5,   0.002237d0, 1.000d0,     &
      0.07920d0,  0.07298d0,  0.7385d0,   0.05532d0,  0.05405d0,  0.002451d0, 0.9975d0,   0.04174d0,   &
      0.8370d0,   0.09674d0,  0.02453d0,  1.000d0,    0.05646d0,  0.9190d0,   0.02160d0,  0.002925d0,  &
      1.000d0,    0.6720d0,   0.2678d0,   0.01183d0,  0.03834d0,  0.01009d0,  0.6850d0,   0.3150d0,    &
      0.4754d0,   0.2813d0,   0.04196d0,  0.1948d0,   0.006629d0, 0.5942d0,   0.4058d0,   0.1961d0,    &
      0.2704d0,   0.07790d0,  0.3738d0,   0.08184d0,  1.000d0,    0.008332d0, 0.09009d0,  0.07433d0,   &
      0.2346d0,   0.5021d0,   0.09057d0,  0.5007d0,   0.4993d0,   0.003254d0, 0.02174d0,  0.1132d0,    &
      0.1137d0,   0.5708d0,   0.1774d0,   0.7170d0,   0.2830d0,   0.005363d0, 0.09668d0,  0.06943d0,   &
      0.8285d0,   1.000d0,    0.5071d0,   0.1118d0,   0.1728d0,   0.1789d0,   0.02944d0,  1.000d0,     &
      0.1422d0,   0.09055d0,  0.1575d0,   0.1668d0,   0.09647d0,  0.2463d0,   0.1003d0,   0.05257d0,   &
      0.01812d0,  0.1249d0,   0.1246d0,   0.1703d0,   0.3181d0,   0.1914d0,   1.000d0,    0.009768d0,  &
      0.1088d0,   0.2201d0,   0.2720d0,   0.2683d0,   0.1210d0,   0.5138d0,   0.4862d0,   0.01178d0,   &
      0.008543d0, 0.1221d0,   0.1263d0,   0.2402d0,   0.1227d0,   0.2911d0,   0.07723d0,  0.04218d0,   &
      0.9578d0,   0.009144d0, 0.006333d0, 0.003291d0, 0.1420d0,   0.07563d0,  0.2406d0,   0.08604d0,   &
      0.3291d0,   0.04755d0,  0.06043d0,  0.5681d0,   0.4319d0,   0.0008457d0,0.02436d0,  0.008572d0,  &
      0.04603d0,  0.06920d0,  0.1859d0,   0.3181d0,   0.3470d0,   1.000d0,    0.0008966d0,0.0008535d0, &
      0.01861d0,  0.2592d0,   0.04028d0,  0.2117d0,   0.2703d0,   0.1064d0,   0.09168d0,  1.000d0,     &
      0.001003d0, 0.0009701d0,0.02357d0,  0.06476d0,  0.07773d0,  0.1120d0,   0.7200d0,   0.0008935d0, &
      0.9991d0,   0.001794d0, 0.002470d0, 0.8832d0,   0.1126d0,   1.000d0,    0.2676d0,   0.1209d0,    &
      0.2375d0,   0.08339d0,  0.1740d0,   0.05845d0,  0.05821d0,  0.02938d0,  0.1465d0,   0.1106d0,    &
      0.1369d0,   0.07358d0,  0.2703d0,   0.2329d0,   0.4748d0,   0.5252d0,   0.001932d0, 0.02134d0,   &
      0.1458d0,   0.2030d0,   0.1562d0,   0.2495d0,   0.2223d0,   1.000d0,    0.0005757d0,0.0009719d0, &
      0.02303d0,  0.1873d0,   0.2542d0,   0.2497d0,   0.2843d0,   1.000d0,    0.001346d0, 0.01569d0,   &
      0.3324d0,   0.2282d0,   0.2709d0,   0.1515d0,   1.000d0,    0.001262d0, 0.02985d0,  0.1411d0,    &
      0.2169d0,   0.1612d0,   0.3200d0,   0.1297d0,   0.9740d0,   0.02604d0,  0.001559d0, 0.05185d0,   &
      0.1844d0,   0.2720d0,   0.1366d0,   0.3537d0,   1.000d0,    0.001175d0, 0.2623d0,   0.1424d0,    &
      0.3066d0,   0.2876d0,   0.3715d0,   0.6285d0,   0.0001934d0,0.01554d0,  0.01572d0,  0.1313d0,    &
      0.1610d0,   0.2632d0,   0.4130d0,   0.3706d0,   0.6294d0,   0.0001363d0,0.007695d0, 0.3278d0,    &
      0.3381d0,   0.2536d0,   0.07269d0,  1.000d0,    0.001466d0, 0.09840d0,  0.1673d0,   0.2303d0,    &
      0.1321d0,   0.3007d0,   0.06976d0,  0.2932d0,   0.7068d0,   0.01378d0,  0.2396d0,   0.2207d0,    &
      0.5259d0,   1.000d0,    1.000d0,    5.310d-5,   0.007114d0, 0.9928d0 /)

  IF (Z .LE. NELEMENTS) THEN
    CALL GetNiso(Z,Niso)
    ALLOCATE(Ziso(Niso),Aiso(Niso),fiso(Niso),Miso(Niso))
    I1 = ELEMENT_INDEX(Z)
    I2 = I1 + Niso - 1
    Ziso = ISOTOPE_Z(I1:I2)
    Aiso = ISOTOPE_A(I1:I2)
    !Jiso = ISOTOPE_J(I1:I2)
    fiso = ISOTOPE_F(I1:I2)
    Miso = IsotopeMass(Ziso,Aiso)
  ELSE
    Niso = 0
    ! Zero-length arrays (nothing to fill in)
    ALLOCATE(Ziso(Niso),Aiso(Niso),fiso(Niso),Miso(Niso))
  END IF
  
END SUBROUTINE


! ----------------------------------------------------------------------
! For the given compound, specified by a list of element atomic numbers
! and the stoichiometry, fills in allocatable arrays containing isotopic
! atomic numbers (Z), atomic masses (A), mass fractions (f), and
! masses (M).
! 
! Input arguments:
!     N          Number of compound elements.
!     Z          Atomic number of elements, array of size [1:N].
!     stoich     Stoichiometry of the compound elements, array of size
!                [1:N].  For example, CF3Cl would have Z={6,9,17} and
!                stoich={1,3,1}.
! Output arguments:
!     Niso       Number of isotopes
!     Ziso       Allocatable array (integer) of isotopes' atomic numbers (Z)
!     Aiso       Allocatable array (integer) of isotopes' atomic masses (A)
!     fiso       Allocatable array of isotopes' mass fractions
!     Miso       Allocatable array of isotopes' nuclear masses [GeV]
! 
SUBROUTINE CompoundIsotopeList(N,Z,stoich,Niso,Ziso,Aiso,fiso,Miso)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: N
  INTEGER, INTENT(IN) :: Z(N),stoich(N)
  INTEGER, INTENT(OUT) :: Niso
  INTEGER, ALLOCATABLE, INTENT(OUT) :: Ziso(:),Aiso(:)
  REAL*8, ALLOCATABLE, INTENT(OUT) :: fiso(:),Miso(:)
  INTEGER :: tempNiso,K,I1,I2
  REAL*8 :: weight(N)
  INTEGER, ALLOCATABLE :: tempZ(:),tempA(:)
  REAL*8, ALLOCATABLE :: tempf(:),tempM(:)
  
  ! Get number of isotopes.
  Niso = 0
  DO K = 1,N
    CALL ElementIsotopeList(Z(K),tempNiso,tempZ,tempA,tempf,tempM)
    Niso = Niso + tempNiso
    weight(K) = stoich(K) * SUM(tempf*tempM)
  END DO
  
  ! Relative weight (by mass) of each element.
  IF (N .GT. 0) THEN
    weight = weight / SUM(weight)
  ELSE
    weight = 1d0
  END IF
  
  ! Allocate and fill in total arrays.
  ! Reweight isotopes' mass fractions by element's mass fraction.
  I1 = 1
  ALLOCATE(Ziso(Niso),Aiso(Niso),fiso(Niso),Miso(Niso))
  DO K = 1,N
    CALL ElementIsotopeList(Z(K),tempNiso,tempZ,tempA,tempf,tempM)
    I2 = I1 + tempNiso - 1
    Ziso(I1:I2) = tempZ
    Aiso(I1:I2) = tempA
    fiso(I1:I2) = weight(K)*tempf
    Miso(I1:I2) = tempM
    I1 = I2 + 1
  END DO
  
END SUBROUTINE


! ----------------------------------------------------------------------
! Determines the mass of the isotope using the semi-empirical mass
! model.  A few low mass (Z <= 2) isotopes are given explicitly to
! avoid accuracy issues in the mass model in those cases.
! 
! Input arguments:
!     Z          Atomic number
!     A          Mass number
! Returns in [GeV].
! 
ELEMENTAL FUNCTION IsotopeMass(Z,A) RESULT(M)
  IMPLICIT NONE
  REAL*8 :: M
  INTEGER, INTENT(IN) :: Z,A
  REAL*8 :: A0,E
  ! 
  ! The semi-empirical mass formula will be used to calculate nuclear
  ! masses, but a few low-mass cases are given explicitly. [GeV]
  REAL*8, PARAMETER :: H2_MASS  = 1.8756129d0  ! Deuteron
  REAL*8, PARAMETER :: H3_MASS  = 2.8089210d0  ! Tritium ("triton")
  REAL*8, PARAMETER :: HE3_MASS = 2.8083915d0  ! Helium-3 ("helion")
  REAL*8, PARAMETER :: HE4_MASS = 3.7273792d0  ! Alpha
  ! 
  ! Semi-empirical mass formula constants
  !   m = Z*Mp + (A-Z)*Mn                 proton+neutron contributions
  !       - aV*A + aS*A^(2/3)             volume & surface terms
  !       + aC*Z^2/A^(1/3)                Coulomb term
  !       + aA*(A-2Z)^2/A                 Pauli (asymmetry) term
  !       + eps_{Z,A} * aP/A^(1/2)        pairing term
  ! where eps_{Z,A} is 0 if A is odd and (-1)^(Z+1) otherwise.
  ! The values below [GeV] are taken from Rohlf (1994).
  REAL*8, PARAMETER :: SEMF_AV = 0.01575d0
  REAL*8, PARAMETER :: SEMF_AS = 0.0178d0
  REAL*8, PARAMETER :: SEMF_AC = 0.000711d0
  REAL*8, PARAMETER :: SEMF_AA = 0.0237d0
  REAL*8, PARAMETER :: SEMF_AP = 0.01118d0

  ! Bad cases
  IF ((A .LE. 0) .OR. (Z .LT. 0)) THEN
    M = 0d0
    RETURN
  END IF
  
  ! Special cases: all natural and/or commonly used isotopes
  ! with Z <= 2 (i.e. nucleons, hydrogen, helium)
  IF (A .EQ. 1) THEN
    IF (Z .EQ. 0) THEN
      M = NEUTRON_MASS
      RETURN
    ELSE IF (Z .EQ. 1) THEN
      M = PROTON_MASS
      RETURN
    END IF
  ELSE IF ((A .EQ. 2) .AND. (Z .EQ. 1)) THEN
    M = H2_MASS
    RETURN
  ELSE IF (A .EQ. 3) THEN
    IF (Z .EQ. 1) THEN
      M = H3_MASS
      RETURN
    ELSE IF (Z .EQ. 2) THEN
      M = HE3_MASS
      RETURN
    END IF
  ELSE IF ((A .EQ. 4) .AND. (Z .EQ. 2)) THEN
    M = HE4_MASS
    RETURN
  END IF
  
  ! Generic semi-empirical mass formula calculation.
  ! Here, E is binding energy.
  A0 = A      ! type conversion
  E = SEMF_AV*A - SEMF_AS*A0**(2d0/3d0) - SEMF_AC*Z**2/A0**(1d0/3d0)    &
      - SEMF_AA*(A-2*Z)**2/A0
  
  ! Pairing term
  IF (MOD(A,2) .EQ. 0) THEN
    IF (MOD(Z,2) .EQ. 0) THEN
      E = E + SEMF_AP/A0**(1d0/2d0)
    ELSE
      E = E - SEMF_AP/A0**(1d0/2d0)
    END IF
  END IF
  
  M = Z*PROTON_MASS + (A-Z)*NEUTRON_MASS - E
  
END FUNCTION


! ----------------------------------------------------------------------
! Calculates the spin-independent weighted form factor Wsi for the given
! isotope at the given momentum transfer q, defined as:
!   Wsi(+1,:) = (1/pi) Z^2 F^2(:)        ! SI proton
!   Wsi( 0,:) = (1/pi) 2*Z*(A-Z) F^2(:)  ! SI crossterm
!   Wsi(-1,:) = (1/pi) (A-Z)^2 F^2(:)    ! SI neutron
! where F^2 is the standard form factor (':' represents momentum array
! q).  Uses the Helm form factor; see:
!   Lewin & Smith, Astropart. Phys. 6, 87 (1996)  [Eqn 4.7]
! 
! Required input arguments:
!     Z,A        The atomic number and mass number of the isotope.
!     N          Number of momentum q values.
!     q          Array of size [1:N] containing momentum values [GeV].
! Required output argument:
!     W          Array of size [-1,1,1:N] to be filled with weighted
!                form factor values [unitless].
! 
PURE SUBROUTINE CalcWSI(Z,A,N,q,W)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: Z,A,N
  REAL*8, INTENT(IN) :: q(1:N)
  REAL*8, INTENT(OUT) :: W(-1:1,1:N)
  INTEGER :: K
  REAL*8 :: c,rn,qrn,qs,F2,weights(-1:1)
  REAL*8, PARAMETER :: Arn = 0.52d0
  REAL*8, PARAMETER :: S   = 0.9d0
  REAL*8, PARAMETER :: C1  = 1.23d0
  REAL*8, PARAMETER :: C2  = -0.60d0
  
  weights(+1) = Z**2      / PI
  weights( 0) = 2*Z*(A-Z) / PI
  weights(-1) = (A-Z)**2  / PI
  
  c  = C1*A**(1d0/3d0) + C2
  rn = SQRT(c**2 + (7d0/3d0)*PI**2*Arn**2 - 5*S**2)
  
  ! Helm FF:  F^2(q) = 9 [j1(q rn)/(q rn)]^2 exp(-q s)
  DO K = 1, N
    qrn = q(K)*rn / HBARC
    qs  = q(K)*S / HBARC
    ! avoid numerical issues for small q by using a Taylor expansion
    IF (qrn .LE. 0.01d0) THEN
      F2 = (1 - qrn**2*((1d0/5d0) - qrn**2*(3d0/175d0))) * EXP(-qs**2)
    ELSE
      F2 = 9 * (SIN(qrn) - qrn*COS(qrn))**2 / qrn**6 * EXP(-qs**2)
    END IF
    W(:,K) = F2 * weights
  END DO
  
END SUBROUTINE


! ----------------------------------------------------------------------
! Calculates the spin-dependent weighted form factor Wsd for the given
! isotope at the given momentum transfer q, defined as:
!   Wsd(+1,:) = 4/(2J+1) Spp(q)          ! SD proton
!   Wsd( 0,:) = 4/(2J+1) Spn(q)          ! SD crossterm
!   Wsd(-1,:) = 4/(2J+1) Snn(q)          ! SD neutron
! where J is the nuclear spin and Spp/Spn/Snn are the spin structure
! functions.  For a comprehensive review of available spin structure
! functions, see:
!   Bednyakov & Simkovic, Phys. Part. Nucl. 37, S106 (2006)
!     [hep-ph/0608097]
! Note that, in the above review, the quantities "ap" & "an" are
! actually the quantities Gp & Gn as used here, not the quantities
! ap & an as used here or in much of the other direct detection
! literature.
! 
! Xenon form factors implemented by Andre Scaffidi.
! 
! NOTE: ONLY A LIMITED SELECTION OF SD FORM FACTORS ARE CURRENTLY
!       IMPLEMENTED.  THEY ARE SIMPLY SET TO ZERO.
! 
! Required input arguments:
!     Z,A        The atomic number and mass number of the isotope.
!     N          Number of momentum q values.
!     q          Array of size [1:N] containing momentum values [GeV].
! Required output argument:
!     W          Array of size [-1,1,1:N] to be filled with weighted
!                form factor values [unitless].
! 
PURE SUBROUTINE CalcWSD(Z,A,N,q,W)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: Z,A,N
  REAL*8, INTENT(IN) :: q(1:N)
  REAL*8, INTENT(OUT) :: W(-1:1,1:N)
  REAL*8 :: J,b,umax,Spp(1:N),Spn(1:N),Snn(1:N)
  
  ! Initialization
  ! Most isotopes have zero spin (=> zero form factor)
  J   = 0d0
  Spp = 0d0
  Snn = 0d0
  Spn = 0d0
  W   = 0d0
  
  ! Fluorine ----------------------------------
  IF (Z .EQ. 9) THEN
    ! Fluorine 19 --------------
    ! From Klos et al. [1304.7684]; see Table 6.
    ! The S11 and S01 (1b+2b) band means are used.
    ! Functional form is Sij(u) = e^{-u} \sum_{k=0} C_{ij,k} u^k
    ! where u = q^2 b^2 / 2.
    IF (A .EQ. 19) THEN
      ! Harmonic oscillator length [GeV^-1]
      ! NOTE: The length is listed incorrectly in 1304.7684 Table 6.
      ! Correct number found in caption of Figure 12.
      b = 1.7608d0 / HBARC
      ! Nuclear spin
      J = 0.5d0
      ! Fits over 0 < u < 5 (probably).  Limit to this range.
      umax = 5d0
      ! Helper routine for exponential-polynomial Sij forms.
      CALL ExponentialPolynomial(N,Spp,Spn,Snn,b,umax,1d0,1d0,1d0,14,   &
               (/  0.108058d0, -0.143789d0,  0.0680848d0, 4.07415d-4,   &
                  -0.0314817d0, 0.0385933d0,-0.0293716d0, 0.0152264d0,  &
                  -5.52655d-3,  1.41965d-3, -2.56989d-4,  3.20688d-5,   &
                  -2.62562d-6,  1.26950d-7, -2.74719d-9  /),            &
               (/  0.0660281d0,-0.137668d0,  0.161957d0, -0.166004d0,   &
                   0.152569d0, -0.111464d0,  0.0609363d0,-0.0246265d0,  &
                   7.35189d-3, -1.61496d-3,  2.57660d-4, -2.90407d-5,   &
                   2.19205d-6, -9.94286d-8,  2.04873d-9  /),            &
               (/  0.167759d0, -0.286581d0,  0.244497d0, -0.176999d0,   &
                   0.134461d0, -0.0921638d0, 0.0494464d0,-0.0198425d0,  &
                   5.87500d-3, -1.26970d-3,  1.96948d-4, -2.12790d-5,   &
                   1.51602d-6, -6.38405d-8,  1.20003d-9  /) )
      
    ! --------------------------
    ! Other isotopes have zero spin
    ELSE
      W = 0d0
      RETURN
    END IF
  
  ! Sodium ------------------------------------
  ELSE IF (Z .EQ. 11) THEN
    ! Sodium 23 ----------------
    ! From Klos et al. [1304.7684]; see Table 7.
    ! The S11 and S01 (1b+2b) band means are used.
    ! Functional form is Sij(u) = e^{-u} \sum_{k=0} C_{ij,k} u^k
    ! where u = q^2 b^2 / 2.
    IF (A .EQ. 23) THEN
      ! Harmonic oscillator length [GeV^-1]
      b = 1.8032d0 / HBARC
      ! Nuclear spin
      J = 1.5d0
      ! Fits over 0 < u < 5 (probably).  Limit to this range.
      umax = 5d0
      ! Helper routine for exponential-polynomial Sij forms.
      CALL ExponentialPolynomial(N,Spp,Spn,Snn,b,umax,1d0,1d0,1d0,9,    &
               (/  0.0325305d0,-0.0433531d0, 0.0319487d0,-5.68858d-3,   &
                   2.67783d-4,  2.44643d-5, -4.79620d-6,  5.39846d-7,   &
                  -3.24691d-8,  8.09358d-10 /),                         &
               (/  0.0127243d0,-0.0248722d0, 0.0275805d0,-0.0135587d0,  &
                   3.93910d-3, -7.28827d-4,  8.82003d-5, -6.80402d-6,   &
                   3.02566d-7, -5.81500d-9  /),                         &
               (/  0.0404609d0,-0.0677623d0, 0.0660458d0,-0.0269059d0,  &
                   7.22875d-3, -1.45367d-3,  2.08818d-4, -1.89848d-5,   &
                   9.43875d-7, -1.93865d-8  /) )
      
    ! --------------------------
    ! Other isotopes have zero spin
    ELSE
      W = 0d0
      RETURN
    END IF
  
  ! Aluminum ----------------------------------
  ELSE IF (Z .EQ. 13) THEN
    ! Aluminum 27 --------------
    ! From Klos et al. [1304.7684]; see Table 7.
    ! The S11 and S01 (1b+2b) band means are used.
    ! Functional form is Sij(u) = e^{-u} \sum_{k=0} C_{ij,k} u^k
    ! where u = q^2 b^2 / 2.
    IF (A .EQ. 27) THEN
      ! Harmonic oscillator length [GeV^-1]
      b = 1.8405d0 / HBARC
      ! Nuclear spin
      J = 2.5d0
      ! Fits over 0 < u < 5 (probably).  Limit to this range.
      umax = 5d0
      ! Helper routine for exponential-polynomial Sij forms.
      CALL ExponentialPolynomial(N,Spp,Spn,Snn,b,umax,1d0,1d0,1d0,9,    &
               (/  0.0888149d0,-0.117822d0,  0.0631336d0,-9.19554d-3,   &
                   5.84421d-4,  5.54484d-4, -1.15453d-4,  1.40388d-5,   &
                  -9.21830d-7,  2.52336d-8  /),                         &
               (/  0.0334384d0,-0.0710220d0, 0.0805917d0,-0.0514533d0,  &
                   0.0221406d0,-6.14292d-3,  1.08899d-3, -1.17175d-4,   &
                   6.93171d-6, -1.71376d-7  /),                         &
               (/  0.108031d0, -0.182354d0,  0.149131d0, -0.0623283d0,  &
                   0.0196187d0,-4.07283d-3,  6.05456d-4, -5.96107d-5,   &
                   3.28620d-6, -7.69545d-8  /) )
      
    ! --------------------------
    ! Other isotopes have zero spin
    ELSE
      W = 0d0
      RETURN
    END IF
  
  ! Silicon -----------------------------------
  ELSE IF (Z .EQ. 14) THEN
    ! Silicon 29 ---------------
    ! From Klos et al. [1304.7684]; see Table 8.
    ! The S11 and S01 (1b+2b) band means are used.
    ! Functional form is Sij(u) = e^{-u} \sum_{k=0} C_{ij,k} u^k
    ! where u = q^2 b^2 / 2.
    IF (A .EQ. 29) THEN
      ! Harmonic oscillator length [GeV^-1]
      b = 1.8575d0 / HBARC
      ! Nuclear spin
      J = 0.5d0
      ! Fits over 0 < u < 5 (probably).  Limit to this range.
      umax = 5d0
      ! Helper routine for exponential-polynomial Sij forms.
      CALL ExponentialPolynomial(N,Spp,Spn,Snn,b,umax,1d0,1d0,1d0,9,    &
               (/  0.0140647d0,-0.0188522d0, 0.0149891d0,-0.00542122d0, &
                   1.17173d-3, -1.15932d-4,  2.47182d-5, -3.04480d-6,   &
                   2.00549d-7, -5.46011d-9  /),                         &
               (/  5.63416d-3, -0.0121901d0, 0.0156006d0,-0.0110712d0,  &
                   4.85653d-3, -1.28086d-3,  2.08618d-4, -2.04711d-5,   &
                   1.10510d-6, -2.47894d-8  /),                         &
               (/ -0.0176295d0, 0.0301066d0,-0.0308628d0, 0.0175759d0,  &
                  -6.78176d-3,  1.69859d-3, -3.00595d-4,  3.39713d-5,   &
                  -2.08093d-6,  5.28653d-8  /) )
      
    ! --------------------------
    ! Other isotopes have zero spin
    ELSE
      W = 0d0
      RETURN
    END IF
  
  ! Germanium ---------------------------------
  ELSE IF (Z .EQ. 32) THEN
    ! Germanium 73 -------------
    ! From Klos et al. [1304.7684]; see Table 5.
    ! The S11 and S01 (1b+2b) band means are used.
    ! Functional form is Sij(u) = e^{-u} \sum_{k=0} C_{ij,k} u^k
    ! where u = q^2 b^2 / 2.
    IF (A .EQ. 73) THEN
      ! Harmonic oscillator length [GeV^-1]
      b = 2.1058d0 / HBARC
      ! Nuclear spin
      J = 4.5d0
      ! Fits over 0 < u < 5 (probably).  Limit to this range.
      umax = 5d0
      ! Helper routine for exponential-polynomial Sij forms.
      CALL ExponentialPolynomial(N,Spp,Spn,Snn,b,umax,1d0,1d0,1d0,9,    &
               (/  0.215608d0, -0.578786d0,  0.698020d0, -0.372000d0,   &
                   0.107576d0, -0.0182408d0, 2.17108d-3, -2.07981d-4,   &
                   1.65907d-5, -5.95664d-7  /),                         &
               (/  0.0972089d0,-0.308986d0,  0.450727d0, -0.337355d0,   &
                   0.154809d0, -0.0469625d0, 9.71560d-3, -1.33058d-3,   &
                   1.09084d-4, -4.02514d-6  /),                         &
               (/ -0.287562d0,  0.844765d0, -1.133659d0,  0.745494d0,   &
                  -0.296646d0,  0.0788570d0,-0.0147852d0, 1.87401d-3,   &
                  -1.42195d-4,  4.72898d-6  /) )
      
    ! --------------------------
    ! Other isotopes have zero spin
    ELSE
      W = 0d0
      RETURN
    END IF
  
  ! Iodine ------------------------------------
  ELSE IF (Z .EQ. 53) THEN
    ! Iodine 127 ---------------
    ! From Klos et al. [1304.7684]; see Table 5.
    ! The S11 and S01 (1b+2b) band means are used.
    ! Functional form is Sij(u) = e^{-u} \sum_{k=0} C_{ij,k} u^k
    ! where u = q^2 b^2 / 2.
    IF (A .EQ. 127) THEN
      ! Harmonic oscillator length [GeV^-1]
      b = 2.2801d0 / HBARC
      ! Nuclear spin
      J = 2.5d0
      ! Fits over 0 < u < 5 (probably).  Limit to this range.
      umax = 5d0
      ! Helper routine for exponential-polynomial Sij forms.
      CALL ExponentialPolynomial(N,Spp,Spn,Snn,b,umax,1d0,1d0,1d0,9,    &
               (/  0.0928480d0,-0.252496d0,  0.351982d0, -0.260427d0,   &
                   0.118280d0, -0.0319614d0, 4.92618d-3, -4.06546d-4,   &
                   1.55818d-5, -1.64934d-7  /),                         &
               (/  0.0389166d0,-0.119307d0,  0.189835d0, -0.168819d0,   &
                   0.0952229d0,-0.0343338d0, 7.86014d-3, -1.11341d-3,   &
                   8.98377d-5, -3.17792d-6  /),                         &
               (/  0.119382d0, -0.345408d0,  0.515816d0, -0.421111d0,   &
                   0.215622d0, -0.0691557d0, 0.0137850d0,-1.68267d-3,   &
                   1.18375d-4, -3.78243d-6  /) )
      
    ! --------------------------
    ! Other isotopes have zero spin
    ELSE
      W = 0d0
      RETURN
    END IF
  
  ! Xenon -------------------------------------
  ELSE IF (Z .EQ. 54) THEN
    ! Xenon 129 ----------------
    ! From Klos et al. [1304.7684]; see Table 1.
    ! The S11 and S01 (1b+2b) band means are used.
    ! Functional form is Sij(u) = e^{-u} \sum_{k=0} C_{ij,k} u^k
    ! where u = q^2 b^2 / 2.
    IF (A .EQ. 129) THEN
      ! Harmonic oscillator length [GeV^-1]
      b = 2.2853d0 / HBARC
      ! Nuclear spin
      J = 0.5d0
      ! Fits over 0 < u < 10 (probably).  Limit to this range.
      umax = 10d0
      ! Helper routine for exponential-polynomial Sij forms.
      CALL ExponentialPolynomial(N,Spp,Spn,Snn,b,umax,1d0,1d0,1d0,9,    &
               (/  0.0547144d0,-0.146407d0,  0.180603d0, -0.125526d0,   &
                   0.0521484d0,-0.0126363d0, 1.76284d-3, -1.32501d-4,   &
                   4.23423d-6, -1.68052d-9  /),                         &
               (/  0.0289650d0,-0.0867525d0, 0.115723d0, -0.0858610d0,  &
                   0.0384596d0,-0.0105918d0, 1.80025d-3, -1.83841d-4,   &
                   1.03293d-5, -2.44338d-7  /),                         &
               (/ -0.0791167d0, 0.225715d0, -0.293581d0,  0.215439d0,   &
                  -0.0959137d0, 0.0260514d0,-4.33883d-3,  4.32823d-4,   &
                  -2.37266d-5,  5.45635d-7  /) )
      
    ! Xenon 131 ----------------
    ! From Klos et al. [1304.7684]; see Table 1.
    ! The S11 and S01 (1b+2b) band means are used.
    ! Functional form is Sij(u) = e^{-u} \sum_{k=0} C_{ij,k} u^k
    ! where u = q^2 b^2 / 2.
    ELSE IF (A .EQ. 131) THEN
      ! Harmonic oscillator length [GeV^-1]
      b = 2.2905d0 / HBARC
      ! Nuclear spin
      J = 1.5d0
      ! Fits over 0 < u < 10 (probably).  Limit to this range.
      umax = 10d0
      ! Helper routine for exponential-polynomial Sij forms.
      CALL ExponentialPolynomial(N,Spp,Spn,Snn,b,umax,1d0,1d0,1d0,9,    &
               (/  0.0417857d0,-0.111132d0,  0.171306d0, -0.132481d0,   &
                   0.0630161d0,-0.0177684d0, 2.82192d-3, -2.32247d-4,   &
                   7.81471d-6,  1.25984d-9  /),                         &
               (/  0.0219206d0,-0.0642919d0, 0.0957262d0,-0.0727452d0,  &
                   0.0338802d0,-9.86454d-3,  1.75888d-3, -1.83629d-4,   &
                   1.01679d-5, -2.25247d-7  /),                         &
               (/ -0.0602463d0, 0.171349d0, -0.265846d0,  0.211589d0,   &
                  -0.103611d0,  0.0311471d0,-5.60929d-3,  5.81416d-4,   &
                  -3.16217d-5,  6.82201d-7  /) )
      
    !! Xenon 129 ----------------
    !! From Menendez et al. [1208.1094]; see Table 1.
    !! The 1b+2b results are taken where possible.
    !! Functional form is Sij(u) = e^{-u} \sum_{k=0} C_{ij,k} u^k
    !! where u = q^2 b^2 / 2.
    !IF (A .EQ. 129) THEN
    !  ! Harmonic oscillator length [GeV^-1]
    !  b = 2.2853d0 / HBARC
    !  ! Nuclear spin
    !  J = 0.5d0
    !  ! Fits over 0 < u < 3 (probably).  Limit to this range.
    !  umax = 3d0
    !  ! Helper routine for exponential-polynomial Sij forms.
    !  CALL ExponentialPolynomial(N,Spp,Spn,Snn,b,umax,1d0,1d0,1d0,9,    &
    !           (/  0.054731d0, -0.146897d0,  0.182479d0,  -0.128112d0,  &
    !               0.0539978d0,-0.0133335d0, 0.00190579d0,-1.48373d-4,  &
    !               5.11732d-6, -2.06597d-8 /),                          &
    !           (/  0.02933d0,  -0.0905396d0, 0.122783d0,  -0.0912046d0, &
    !               0.0401076d0,-0.010598d0,  0.00168737d0,-1.56768d-4,  &
    !               7.69202d-6, -1.48874d-7 /),                          &
    !           (/ -0.0796645d0, 0.231997d0, -0.304198d0,   0.222024d0,  &
    !              -0.096693d0,  0.0251835d0,-0.00392356d0, 3.53343d-4,  &
    !              -1.65058d-5,  2.88576d-7 /) )
      
    !! Xenon 131 ----------------
    !! From Menendez et al. [1208.1094]; see Table 1.
    !! The 1b+2b results are taken where possible.
    !! Functional form is Sij(u) = e^{-u} \sum_{k=0} C_{ij,k} u^k
    !! where u = q^2 b^2 / 2.
    !ELSE IF (A .EQ. 131) THEN
    !  ! Harmonic oscillator length [GeV^-1]
    !  b = 2.2905d0 / HBARC
    !  ! Nuclear spin
    !  J=1.5d0
    !  ! Fits over 0 < u < 3 (probably).  Limit to this range.
    !  umax = 3d0
    !  ! Helper routine for exponential-polynomial Sij forms.
    !  CALL ExponentialPolynomial(N,Spp,Spn,Snn,b,umax,1d0,1d0,1d0,9,    &
    !           (/  0.0417889d0,-0.111171d0,  0.171966d0,  -0.133219d0,  &
    !               0.0633805d0,-0.0178388d0, 0.00282476d0,-2.31681d-4,  &
    !               7.78223d-6, -4.49287d-10 /),                         &
    !           (/  0.022446d0, -0.0733931d0, 0.110509d0,  -0.0868752d0, &
    !               0.0405399d0,-0.0113544d0, 0.00187572d0,-1.75285d-4,  &
    !               8.40043d-6, -1.53632d-7  /),                         &
    !           (/ -0.0608808d0, 0.181473d0, -0.272533d0,   0.211776d0,  &
    !              -0.0985956d0, 0.027438d0, -0.0044424d0,  3.97619d-4,  &
    !              -1.74758d-5,  2.55979d-7  /) )
      
    ! --------------------------
    ! Other isotopes have zero spin
    ELSE
      W = 0d0
      RETURN
    END IF
  
  ! Zero spin or unimplemented isotope ---------
  ELSE
    W = 0d0
    RETURN
  ENDIF
  
  ! Weighted form factors
  W(+1,:) = (4d0 / (2d0*J + 1d0)) * Spp    ! SD proton
  W( 0,:) = (4d0 / (2d0*J + 1d0)) * Spn    ! SD crossterm
  W(-1,:) = (4d0 / (2d0*J + 1d0)) * Snn    ! SD neutron
  
  
  CONTAINS
  
  ! --------------------------------------------
  ! Calculates form factors of the form:
  !   Sij(u) = e^{-u} A_ij \sum_{k=0} C_{ij,k} u^k
  ! where u = q^2 b^2 / 2 and A = (A00,A11,A01).  Sets to zero
  ! when u > umax.  Fills the arrays S00, S11, & S01 in the
  ! parent routine.
  PURE SUBROUTINE ExponentialPolynomial(N,Spp,Spn,Snn,b,umax,           &
                                        A00,A11,A01,NC,C00,C11,C01)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: N,NC
    REAL*8, INTENT(IN) :: b,umax,A00,A11,A01,C00(0:NC),C11(0:NC),C01(0:NC)
    REAL*8, INTENT(OUT) :: Spp(1:N),Spn(1:N),Snn(1:N)
    INTEGER :: I,K
    REAL*8 :: u,uk,expu,S00,S11,S01
    DO I=1,N
      u  = 0.5d0 * (q(I)*b)**2
      uk = 1
      S00 = 0d0
      S11 = 0d0
      S01 = 0d0
      ! Conservatively set to zero outside range of validity.
      IF (u .LE. umax) THEN
        DO K=0,NC
          S00 = S00 + uk * C00(K)
          S11 = S11 + uk * C11(K)
          S01 = S01 + uk * C01(K)
          uk = u*uk
        END DO
        expu = EXP(-u)
        S00 = expu * A00 * S00
        S11 = expu * A11 * S11
        S01 = expu * A01 * S01
      END IF
      ! Basis transformation: a0=ap+an, a1=ap-an
      Spp(I) = S00 + S11 + S01
      Snn(I) = S00 + S11 - S01
      Spn(I) = 2 * (S00 - S11)
    END DO
    
  END SUBROUTINE ExponentialPolynomial
  
  
END SUBROUTINE


! ----------------------------------------------------------------------
! INTERFACE NAME: EToQ
! Calculates the momentum transfer q [GeV] corresponding to a nuclear
! recoil energy E [keV] for the given isotope mass m [GeV].
! 
! This is the scalar version (single mass and energy).
! 
! Input arguments:
!   E           Recoil energy [keV]
!   Miso        Isotope mass [GeV]
! 
ELEMENTAL FUNCTION EToQ(E,Miso) RESULT(q)
  IMPLICIT NONE
  REAL*8 :: Q
  REAL*8, INTENT(IN) :: E,Miso
  ! Factor of 1d-6 to take keV -> GeV
  q = SQRT(2*Miso*(1d-6*E))
END FUNCTION


END MODULE
