MODULE DDConstants

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! DDConstants
!    Constants, global variables and types for DDCalc
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

! MATH CONSTANTS -------------------------------------------------------

REAL*8, PUBLIC, PARAMETER :: PI         = 3.1415926535897932d0   ! Pi
REAL*8, PUBLIC, PARAMETER :: TWOPI      = 6.2831853071795865d0   ! 2*Pi
REAL*8, PUBLIC, PARAMETER :: HALFPI     = 1.5707963267948966d0   ! Pi/2
REAL*8, PUBLIC, PARAMETER :: FOURPI     =12.5663706143591730d0   ! 4*Pi
REAL*8, PUBLIC, PARAMETER :: SQRTPI     = 1.7724538509055160d0   ! Sqrt(Pi)
REAL*8, PUBLIC, PARAMETER :: SQRT2PI    = 2.5066282746310005d0   ! Sqrt(2*Pi)
REAL*8, PUBLIC, PARAMETER :: INVPI      = 0.31830988618379067d0  ! 1/Pi
REAL*8, PUBLIC, PARAMETER :: INV2PI     = 0.15915494309189534d0  ! 1/2Pi
REAL*8, PUBLIC, PARAMETER :: INVSQRTPI  = 0.56418958354775629d0  ! 1/Sqrt(Pi)
REAL*8, PUBLIC, PARAMETER :: INVSQRT2PI = 0.39894228040143268d0  ! 1/Sqrt(2*Pi)

REAL*8, PUBLIC, PARAMETER :: SQRT2  = 1.4142135623730950d0   ! Sqrt(2)
REAL*8, PUBLIC, PARAMETER :: SQRT3  = 1.7320508075688773d0   ! Sqrt(3)


! PHYSICS CONSTANTS ----------------------------------------------------

! Speed of light [m/s]
REAL*8, PUBLIC, PARAMETER :: SPEED_OF_LIGHT = 2.99792458d8

! Planck constant times speed of light [GeV fm]
REAL*8, PUBLIC, PARAMETER :: HBARC = 0.1973269718d0

! Fermi coupling constant, in units of /(hbar c)^3 [GeV^-2]
REAL*8, PUBLIC, PARAMETER :: FERMI_COUPLING_CONSTANT = 1.1663787d-5

! Proton and neutron masses [GeV]
REAL*8, PUBLIC, PARAMETER :: PROTON_MASS    = 0.9382720d0
REAL*8, PUBLIC, PARAMETER :: NEUTRON_MASS   = 0.9395654d0


! CODE CONSTANTS -------------------------------------------------------

! Version of this software
CHARACTER*(*), PUBLIC, PARAMETER :: VERSION_STRING = '1.2.0'

! Verbosity
! Affects the level of output.  Generally, higher magnitude
! means more output; positive will include headers, while
! negative will not.
INTEGER, PUBLIC :: VerbosityLevel = 1


END MODULE
