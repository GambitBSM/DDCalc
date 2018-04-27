MODULE DDNuclear

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! DDNuclear
!    Routines for determining isotope and nuclear properties, including
!    isotopic compositions, nuclear masses, and form factors (both
!    spin-independent and spin-dependent).
!    [NOTE: LIMITED IMPLEMENTATION FOR SD CASE -- MOSTLY SET TO ZERO.]
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

USE DDConstants
USE DDInput

IMPLICIT NONE
PRIVATE

PUBLIC :: IsotopeMass, EtoQ, ElementIsotopeList, GetNiso, CompoundIsotopeList, CalcF2, CalcWTilde

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
!     Jiso       Allocatable array of isotopes' toal spin (J)
!     fiso       Allocatable array of isotopes' mass fractions
!     Miso       Allocatable array of isotopes' nuclear masses [GeV]
! 
SUBROUTINE ElementIsotopeList(Z,Niso,Ziso,Aiso,Jiso,fiso,Miso)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: Z
  INTEGER, INTENT(OUT) :: Niso
  INTEGER, ALLOCATABLE, INTENT(OUT) :: Ziso(:),Aiso(:)
  REAL*8, ALLOCATABLE, INTENT(OUT) :: Jiso(:),fiso(:),Miso(:)
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
    ALLOCATE(Ziso(Niso),Aiso(Niso),Jiso(Niso),fiso(Niso),Miso(Niso))
    I1 = ELEMENT_INDEX(Z)
    I2 = I1 + Niso - 1
    Ziso = ISOTOPE_Z(I1:I2)
    Aiso = ISOTOPE_A(I1:I2)
    Jiso = ISOTOPE_J(I1:I2)
    fiso = ISOTOPE_F(I1:I2)
    Miso = IsotopeMass(Ziso,Aiso)
  ELSE
    Niso = 0
    ! Zero-length arrays (nothing to fill in)
    ALLOCATE(Ziso(Niso),Aiso(Niso),Jiso(Niso),fiso(Niso),Miso(Niso))
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
!     Jiso       Allocatable array of isotopes' total spin (J)
!     fiso       Allocatable array of isotopes' mass fractions
!     Miso       Allocatable array of isotopes' nuclear masses [GeV]
! 
SUBROUTINE CompoundIsotopeList(N,Z,stoich,Niso,Ziso,Aiso,Jiso,fiso,Miso)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: N
  INTEGER, INTENT(IN) :: Z(N),stoich(N)
  INTEGER, INTENT(OUT) :: Niso
  INTEGER, ALLOCATABLE, INTENT(OUT) :: Ziso(:),Aiso(:)
  REAL*8, ALLOCATABLE, INTENT(OUT) :: Jiso(:),fiso(:),Miso(:)
  INTEGER :: tempNiso,K,I1,I2
  REAL*8 :: weight(N)
  INTEGER, ALLOCATABLE :: tempZ(:),tempA(:)
  REAL*8, ALLOCATABLE :: tempJ(:),tempf(:),tempM(:)
  
  ! Get number of isotopes.
  Niso = 0
  DO K = 1,N
    CALL ElementIsotopeList(Z(K),tempNiso,tempZ,tempA,tempJ,tempf,tempM)
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
  ALLOCATE(Ziso(Niso),Aiso(Niso),Jiso(Niso),fiso(Niso),Miso(Niso))
  DO K = 1,N
    CALL ElementIsotopeList(Z(K),tempNiso,tempZ,tempA,tempJ,tempf,tempM)
    I2 = I1 + tempNiso - 1
    Ziso(I1:I2) = tempZ
    Aiso(I1:I2) = tempA
    Jiso(I1:I2) = tempJ
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
PURE SUBROUTINE CalcF2(A,N,q,F2)

  IMPLICIT NONE
  INTEGER, INTENT(IN) :: A,N
  REAL*8, INTENT(IN) :: q(1:N)
  REAL*8, INTENT(OUT) :: F2(1:N)
  INTEGER :: K
  REAL*8 :: c,rn,qrn,qs
  REAL*8, PARAMETER :: Arn = 0.52d0
  REAL*8, PARAMETER :: S   = 0.9d0
  REAL*8, PARAMETER :: C1  = 1.23d0
  REAL*8, PARAMETER :: C2  = -0.60d0
    
  c  = C1*A**(1d0/3d0) + C2
  rn = SQRT(c**2 + (7d0/3d0)*PI**2*Arn**2 - 5*S**2)
  
  ! Helm FF:  F^2(q) = 9 [j1(q rn)/(q rn)]^2 exp(-q s)
  DO K = 1, N
    qrn = q(K)*rn / HBARC
    qs  = q(K)*S / HBARC
    ! avoid numerical issues for small q by using a Taylor expansion
    IF (qrn .LE. 0.01d0) THEN
      F2(K) = (1 - qrn**2*((1d0/5d0) - qrn**2*(3d0/175d0))) * EXP(-qs**2)
    ELSE
      F2(K) = 9 * (SIN(qrn) - qrn*COS(qrn))**2 / qrn**6 * EXP(-qs**2)
    END IF
    
  END DO
  
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

SUBROUTINE LoadWbarFile(Z,A,Wbar,success)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: Z,A
  REAL*8, INTENT(INOUT) :: Wbar(8,4,11)
  LOGICAL, INTENT(INOUT) :: success
  INTEGER :: Nrow,Ncol
  REAL*8, ALLOCATABLE :: data(:,:)
  LOGICAL :: status

  character(len=1024) :: filename
  character(len=1024) :: format_string
  
  INTEGER :: Zlength, Alength
  IF (Z < 10) THEN
    Zlength = 1
  ELSE
    Zlength = 2
  END IF

  IF (A < 10) THEN
    Alength = 1
  ELSE IF (A < 100) THEN
    Alength = 2
  ELSE
    Alength = 3
  END IF

  WRITE (format_string,"(A6,I1,A5,I1,A4)") "(A10,I",Zlength,",A1,I",Alength,",A4)"

  WRITE (filename,format_string) "Wbar/",Z,"_",A,".dat"

  ! Load table from file
  CALL LoadTable(file=TRIM(filename),Nrow=Nrow,Ncol=Ncol,data=data,status=status)
  IF (.NOT. status) THEN
    success = .FALSE.
    Wbar(:,:,:) = 0
  ELSE IF ((Ncol .NE. 1) .OR. (Nrow .NE. 352)) THEN
    WRITE(0,*) 'ERROR: Number of entries in file ' // TRIM(filename) // ' does not agree with expectation.'
    STOP
  ELSE
    success = .TRUE.
    Wbar = reshape(data,shape(Wbar),order = (/3,2,1/))
  END IF

END SUBROUTINE




SUBROUTINE CalcWTilde(Z,A,J,NE,qArray,WT)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: Z,A,NE
  REAL*8, INTENT(IN) :: J
  REAL*8, INTENT(IN) :: qArray(NE)
  REAL*8, INTENT(INOUT) :: WT(1:8,1:4,NE)

  REAL*8 :: yArray(NE)
  REAL*8 :: wbar(1:8,1:4,1:11)
  REAL*8 :: F2(1:NE)
  LOGICAL :: success
  INTEGER :: alpha, t_tp, k
 
  CALL LoadWbarFile(Z,A,wbar,success)

  IF (success) THEN
    yArray = 25.6819 * (41.467/(45.0*A**(-1.0/3) - 25.0*A**(-2.0/3))) * qArray**2/4.0 
	! Parameter y following 1308.6288

     DO alpha = 1,8
       DO t_tp = 1,4
          WT(alpha,t_tp,:) = 0.0d0
          DO k = 1,11
             WT(alpha,t_tp,:) = WT(alpha,t_tp,:) + wbar(alpha,t_tp,k)*(yArray**(k-1))
          END DO
          WT(alpha,t_tp,:) = WT(alpha,t_tp,:)*4*PI/(2*J+1) * exp(-2*yArray)
       END DO
     END DO
   ELSE ! Use Helm form factor for SI interaction and set the rest to zero
     WT(:,:,:) = 0
     CALL CalcF2(A,NE,qArray,F2)
     WT(1,1,:) = 1.d0/4.d0 * A**2 * F2(:)
     WT(1,2,:) = 1.d0/4.d0 * A*(2*Z-A) * F2(:)
     WT(1,3,:) = 1.d0/4.d0 * A*(2*Z-A) * F2(:)
     WT(1,4,:) = 1.d0/4.d0 * (2*Z-A)**2 * F2(:)
   END IF

END SUBROUTINE

END MODULE
