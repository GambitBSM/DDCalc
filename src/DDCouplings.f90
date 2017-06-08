MODULE DDCouplings

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! DDCouplings
!  * Routines for converting between WIMP-nucleon couplings G, f & a,
!    with different notation corresponding to different normalization.
!  * Routines to convert between WIMP-nucleon couplings G and WIMP-
!    nucleon cross-sections sigma.
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

USE DDConstants

IMPLICIT NONE

CONTAINS


!=======================================================================
! COUPLING CONVERSION ROUTINES
!=======================================================================

!-----------------------------------------------------------------------
! Converts spin-independent coupling f to G, related by:
!   Gp = 2*fp    Gn = 2*fn
! 
! Input arguments:
!   f           SI coupling [GeV^-2]
! Returns in [GeV^-2].
! 
ELEMENTAL FUNCTION FToG(f) RESULT(G)
  IMPLICIT NONE
  REAL*8 :: G
  REAL*8, INTENT(IN) :: f
  REAL*8, PARAMETER :: F_SCALE = 2d0
  G = f * F_SCALE
END FUNCTION


!-----------------------------------------------------------------------
! Converts spin-independent coupling G to f, related by:
!   Gp = 2*fp    Gn = 2*fn
! 
! Input arguments:
!   G           SI coupling [GeV^-2]
! Returns in [GeV^-2].
! 
ELEMENTAL FUNCTION GToF(G) RESULT(f)
  IMPLICIT NONE
  REAL*8 :: f
  REAL*8, INTENT(IN) :: G
  REAL*8, PARAMETER :: F_SCALE = 2d0
  f = G / F_SCALE
END FUNCTION


!-----------------------------------------------------------------------
! Converts spin-dependent coupling a to G, related by:
!   Gp = 2*sqrt{2}*G_F*ap    Gn = 2*sqrt{2}*G_F*an
! 
! Input arguments:
!   a           SD coupling [unitless]
! Returns in [GeV^-2].
! 
ELEMENTAL FUNCTION AToG(a) RESULT(G)
  IMPLICIT NONE
  REAL*8 :: G
  REAL*8, INTENT(IN) :: a
  REAL*8, PARAMETER :: A_SCALE = SQRT(2d0)*2d0*FERMI_COUPLING_CONSTANT
  G = a * A_SCALE
END FUNCTION


!-----------------------------------------------------------------------
! Converts spin-dependent coupling G to a, related by:
!   Gp = 2*sqrt{2}*G_F*ap    Gn = 2*sqrt{2}*G_F*an
! 
! Input arguments:
!   G           SD coupling [GeV^-2]
! Returns in [unitless].
! 
ELEMENTAL FUNCTION GToA(G) RESULT(a)
  IMPLICIT NONE
  REAL*8 :: a
  REAL*8, INTENT(IN) :: G
  REAL*8, PARAMETER :: A_SCALE = SQRT(2d0)*2d0*FERMI_COUPLING_CONSTANT
  a = G / A_SCALE
END FUNCTION



!=======================================================================
! CROSS-SECTION/COUPLING CONVERSION ROUTINES
!=======================================================================

!-----------------------------------------------------------------------
! Converts WIMP-proton spin-independent cross-section to G, related by:
!   sigmapSI = \mu^2/\pi Gp^2
! where \mu is WIMP-proton reduced mass.
! 
! Input arguments:
!   m           WIMP mass [GeV]
!   sigma       SI cross-section [pb].  A negative value means
!               the coupling should be set negative.
! Returns in [GeV^-2].
! 
ELEMENTAL FUNCTION SigmapSIToGp(m,sigma) RESULT(G)
  IMPLICIT NONE
  REAL*8 :: G
  REAL*8, INTENT(IN) :: m,sigma
  REAL*8, PARAMETER :: Mp = PROTON_MASS
  ! Factor of 1d-10 takes pb to fm^2
  IF (sigma .GE. 0d0) THEN
    G = SQRTPI * (m+Mp)/(m*Mp) * SQRT(1d-10*sigma) / HBARC
  ELSE
    G = - SQRTPI * (m+Mp)/(m*Mp) * SQRT(-1d-10*sigma) / HBARC
  END IF
END FUNCTION


!-----------------------------------------------------------------------
! Converts WIMP-neutron spin-independent cross-section to G, related by:
!   sigmanSI = \mu^2/\pi Gn^2
! where \mu is WIMP-neutron reduced mass.
! 
! Input arguments:
!   m           WIMP mass [GeV]
!   sigma       SI cross-section [pb].  A negative value means
!               the coupling should be set negative.
! Returns in [GeV^-2].
! 
ELEMENTAL FUNCTION SigmanSIToGn(m,sigma) RESULT(G)
  IMPLICIT NONE
  REAL*8 :: G
  REAL*8, INTENT(IN) :: m,sigma
  REAL*8, PARAMETER :: Mn = NEUTRON_MASS
  ! Factor of 1d-10 takes pb to fm^2
  IF (sigma .GE. 0d0) THEN
    G = SQRTPI * (m+Mn)/(m*Mn) * SQRT(1d-10*sigma) / HBARC
  ELSE
    G = - SQRTPI * (m+Mn)/(m*Mn) * SQRT(-1d-10*sigma) / HBARC
  END IF
END FUNCTION


!-----------------------------------------------------------------------
! Converts G to WIMP-proton spin-independent cross-section, related by:
!   sigmapSI = \mu^2/\pi Gp^2
! where \mu is WIMP-proton reduced mass.
! 
! Input arguments:
!   m           WIMP mass [GeV]
!   G           SI coupling [GeV^-2]
! Returns in [pb].
! 
ELEMENTAL FUNCTION GpToSigmapSI(m,G) RESULT(sigma)
  IMPLICIT NONE
  REAL*8 :: sigma
  REAL*8, INTENT(IN) :: m,G
  REAL*8, PARAMETER :: Mp = PROTON_MASS
  ! Factor of 1d10 takes fm^2 to pb
  sigma = 1d10 * INVPI * ((m*Mp/(m+Mp))*G*HBARC)**2
END FUNCTION


!-----------------------------------------------------------------------
! Converts G to WIMP-neutron spin-independent cross-section, related by:
!   sigmanSI = \mu^2/\pi Gn^2
! where \mu is WIMP-neutron reduced mass.
! 
! Input arguments:
!   m           WIMP mass [GeV]
!   G           SI coupling [GeV^-2]
! Returns in [pb].
! 
ELEMENTAL FUNCTION GnToSigmanSI(m,G) RESULT(sigma)
  IMPLICIT NONE
  REAL*8 :: sigma
  REAL*8, INTENT(IN) :: m,G
  REAL*8, PARAMETER :: Mn = NEUTRON_MASS
  ! Factor of 1d10 takes fm^2 to pb
  sigma = 1d10 * INVPI * ((m*Mn/(m+Mn))*G*HBARC)**2
END FUNCTION


!-----------------------------------------------------------------------
! Converts WIMP-proton spin-dependent cross-section to G, related by:
!   sigmapSD = 3\mu^2/\pi Gp^2
! where \mu is WIMP-proton reduced mass.
! 
! Input arguments:
!   m           WIMP mass [GeV]
!   sigma       SD cross-section [pb].  A negative value means
!               the coupling should be set negative.
! Returns in [GeV^-2].
! 
ELEMENTAL FUNCTION SigmapSDToGp(m,sigma) RESULT(G)
  IMPLICIT NONE
  REAL*8 :: G
  REAL*8, INTENT(IN) :: m,sigma
  REAL*8, PARAMETER :: Mp = PROTON_MASS
  ! Factor of 1d-10 takes pb to fm^2
  IF (sigma .GE. 0d0) THEN
    G = SQRTPI/SQRT3 * (m+Mp)/(m*Mp) * SQRT(1d-10*sigma) / HBARC
  ELSE
    G = - SQRTPI/SQRT3 * (m+Mp)/(m*Mp) * SQRT(-1d-10*sigma) / HBARC
  END IF
END FUNCTION


!-----------------------------------------------------------------------
! Converts WIMP-neutron spin-dependent cross-section to G, related by:
!   sigmanSD = 3\mu^2/\pi Gn^2
! where \mu is WIMP-neutron reduced mass.
! 
! Input arguments:
!   m           WIMP mass [GeV]
!   sigma       SD cross-section [pb].  A negative value means
!               the coupling should be set negative.
! Returns in [GeV^-2].
! 
ELEMENTAL FUNCTION SigmanSDToGn(m,sigma) RESULT(G)
  IMPLICIT NONE
  REAL*8 :: G
  REAL*8, INTENT(IN) :: m,sigma
  REAL*8, PARAMETER :: Mn = NEUTRON_MASS
  ! Factor of 1d-10 takes pb to fm^2
  IF (sigma .GE. 0d0) THEN
    G = SQRTPI/SQRT3 * (m+Mn)/(m*Mn) * SQRT(1d-10*sigma) / HBARC
  ELSE
    G = - SQRTPI/SQRT3 * (m+Mn)/(m*Mn) * SQRT(-1d-10*sigma) / HBARC
  END IF
END FUNCTION


!-----------------------------------------------------------------------
! Converts G to WIMP-proton spin-dependent cross-section, related by:
!   sigmapSD = 3\mu^2/\pi Gp^2
! where \mu is WIMP-proton reduced mass.
! 
! Input arguments:
!   m           WIMP mass [GeV]
!   G           SD coupling [GeV^-2]
! Returns in [pb].
! 
ELEMENTAL FUNCTION GpToSigmapSD(m,G) RESULT(sigma)
  IMPLICIT NONE
  REAL*8 :: sigma
  REAL*8, INTENT(IN) :: m,G
  REAL*8, PARAMETER :: Mp = PROTON_MASS
  ! Factor of 1d10 takes fm^2 to pb
  sigma = 1d10 * 3*INVPI * ((m*Mp/(m+Mp))*G*HBARC)**2
END FUNCTION


!-----------------------------------------------------------------------
! Converts G to WIMP-neutron spin-dependent cross-section, related by:
!   sigmanSD = 3\mu^2/\pi Gn^2
! where \mu is WIMP-neutron reduced mass.
! 
! Input arguments:
!   m           WIMP mass [GeV]
!   G           SD coupling [GeV^-2]
! Returns in [pb].
! 
ELEMENTAL FUNCTION GnToSigmanSD(m,G) RESULT(sigma)
  IMPLICIT NONE
  REAL*8 :: sigma
  REAL*8, INTENT(IN) :: m,G
  REAL*8, PARAMETER :: Mn = NEUTRON_MASS
  ! Factor of 1d10 takes fm^2 to pb
  sigma = 1d10 * 3*INVPI * ((m*Mn/(m+Mn))*G*HBARC)**2
END FUNCTION

!-----------------------------------------------------------------------
! Composite functions converting couplings to cross sections
ELEMENTAL FUNCTION FpToSigmapSI(m,f) RESULT(sigma)
  IMPLICIT NONE
  REAL*8 :: sigma
  REAL*8, INTENT(IN) :: m,f
  sigma = GpToSigmapSI(m,FtoG(f))
END FUNCTION

ELEMENTAL FUNCTION FnToSigmanSI(m,f) RESULT(sigma)
  IMPLICIT NONE
  REAL*8 :: sigma
  REAL*8, INTENT(IN) :: m,f
  sigma = GnToSigmanSI(m,FtoG(f))
END FUNCTION

ELEMENTAL FUNCTION ApToSigmapSD(m,a) RESULT(sigma)
  IMPLICIT NONE
  REAL*8 :: sigma
  REAL*8, INTENT(IN) :: m,a
  sigma = GpToSigmapSD(m,AtoG(a))
END FUNCTION

ELEMENTAL FUNCTION AnToSigmanSD(m,a) RESULT(sigma)
  IMPLICIT NONE
  REAL*8 :: sigma
  REAL*8, INTENT(IN) :: m,a
  sigma = GnToSigmanSD(m,AtoG(a))
END FUNCTION

!-----------------------------------------------------------------------
! Composite functions converting cross sections to couplings
ELEMENTAL FUNCTION SigmapSItoFp(m,sigma) RESULT(f)
  IMPLICIT NONE
  REAL*8 :: f
  REAL*8, INTENT(IN) :: m,sigma
  f = GtoF(SigmapSIToGp(m,sigma))
END FUNCTION

ELEMENTAL FUNCTION SigmanSItoFn(m,sigma) RESULT(f)
  IMPLICIT NONE
  REAL*8 :: f
  REAL*8, INTENT(IN) :: m,sigma
  f = GtoF(SigmanSIToGn(m,sigma))
END FUNCTION

ELEMENTAL FUNCTION SigmapSDtoAp(m,sigma) RESULT(a)
  IMPLICIT NONE
  REAL*8 :: a
  REAL*8, INTENT(IN) :: m,sigma
  a = GtoA(SigmapSDToGp(m,sigma))
END FUNCTION

ELEMENTAL FUNCTION SigmanSDtoAn(m,sigma) RESULT(a)
  IMPLICIT NONE
  REAL*8 :: a
  REAL*8, INTENT(IN) :: m,sigma
  a = GtoA(SigmanSDToGn(m,sigma))
END FUNCTION

END MODULE
