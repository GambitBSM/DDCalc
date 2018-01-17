MODULE DDHalo

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! DDHalo
!    Routines for initializing and modifying the dark matter halo
!    distribution, as well as performing the mean inverse speed (g)
!    and mean speed (h) calculation needed for recoil rates.
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

USE DDConstants
USE DDTypes
USE DDUtils
USE DDNumerical
USE DDInput

IMPLICIT NONE
PRIVATE

PUBLIC :: DDCalc_GetHalo,DDCalc_SetHalo, &
          DDCalc_InitHalo,               &
          EToVmin, MeanInverseSpeed, C_DDCalc_InitHalo, &
          MeanSpeed
INTERFACE DDCalc_GetHalo
  MODULE PROCEDURE GetHalo
END INTERFACE
INTERFACE DDCalc_SetHalo
  MODULE PROCEDURE SetHalo
END INTERFACE
INTERFACE DDCalc_InitHalo
  MODULE PROCEDURE InitHalo
END INTERFACE

! Parameters describing the dark matter halo.  Only the Standard Halo
! Model (SHM) can be used for the velocity distribution (i.e. Maxwell-
! Boltzmann distribution with a finite cutoff).

! Local dark matter halo density [GeV/cm^3]:
!   0.3 [standard (old)]
! * Catena & Ullio, JCAP 1008, 004 (2010) [arXiv:0907.0018]
!   For spherical halos, not including structure
!     Einasto profile: 0.385 +/- 0.027
!     NFW profile:     0.389 +/- 0.025
! * Weber & de Boer, arXiv:0910.4272
!     0.2 - 0.4 (depending on model)
! * Salucci et al., arXiv:1003.3101
!   Model independent technique?
!     0.430 +/- 0.113 (alpha) +/- 0.096 (r)
! * Pato et al., PRD 82, 023531 (2010) [arXiv:1006.1322]
!   Density at disk may be 1.01 - 1.41 times larger than shell
!   averaged quantity, so above measurements might be underestimates
!   of local density.
! DEFAULT: 0.4 GeV/cm^3

! Sun's peculiar velocity [km/s]:
! motion relative to local standard of rest (LSR)
! * Mignard, Astron. Astrophys. 354, 522 (2000)
! * Schoenrich, Binney & Dehnen, arXiv:0912.3693
! DEFAULT: (11,12,7) km/s

! Disk rotation velocity [km/s]:
! * Kerr & Lynden-Bell, MNRAS 221, 1023 (1986)
!     220 [standard]
!     222 +/- 20 (average over multiple measurements, could be biased
!                 by systematics)
! * Reid et al., Astrophys. J. 700, 137 (2009) [arXiv:0902.3913]
!   Estimate based on masers.
!     254 +/- 16
! * McMillan & Binney, MNRAS 402, 934 (2010) [arXiv:0907.4685]
!   Reanalysis of Reid et al. masers.  Range of estimates based on
!   various models; suggest Sun's velocity with respect to LSR should
!   be modified.
!     200 +/- 20 to 279 +/- 33
! * Bovy, Hogg & Rix, ApJ 704, 1704 (2009) [arXiv:0907.5423]
!     244 +/- 13 (masers only)
!     236 +/- 11 (combined estimate)
! DEFAULT: 235 km/s

! The Local Standard of Rest (LSR) [km/s], which we take to be
! (0,vrot,0).
! DEFAULT: (0,235,0) km/s

! Sun's velocity vector relative to the galactic rest frame [km/s],
! sum of LSR and peculiar velocity, where LSR = (0,vrot,0):
! DEFAULT: (0,235,0) + (11,12,7) km/s

! Sun's speed relative to the galactic rest frame [km/s].
! Equal to magnitude of Sun's velocity.
! DEFAULT: sqrt{11^2 + (235+12)^2 + 7^2} km/s

! Most probable speed (characterizing velocity dispersion) [km/s]:
! Typically similar to rotation velocity.
!     vrms = sqrt(3/2) v0    [rms velocity]
!     vmp  = v0              [most probably velocity]
!     vave = sqrt(4/pi) v0   [mean velocity]
! DEFAULT: 235 km/s

! Local escape velocity [km/s]:
!   650 [standard (old)]
! * Smith et al., MNRAS 379, 755 (2007) [astro-ph/0611671]
!   Note from Fig 7 that the distribution is asymmetric.  The following
!   results assume vrot = 220.
!     544 (mean), 498 - 608 (90% CL)
!     462 - 640 (90% CL when also fitting parameter k)
! DEFAULT: 550 km/s


CONTAINS


! ----------------------------------------------------------------------
! Get various halo quantities.
! 
! Optional output arguments regarding galactic motions:
!   vrot        Local galactic disk rotation speed [km/s].
!   vlsr        Local standard of rest velocity vector (array of size 3)
!               [km/s], defined relative to galactic rest frame.
!   vpec        Sun's peculiar velocity vector (array of size 3) [km/s],
!               defined relative to local standard of rest.
!   vsun        Sun's velocity vector (array of size 3) [km/s], defined
!               relative to galactic rest frame.
! Optional output arguments regarding dark matter density:
!   rho         Local dark matter density [GeV/cm^3].
! Optional output arguments regarding SHM distribution, a truncated
! Maxwell-Boltzmann ("MB"):
!   vbulk       Bulk velocity of dark matter (array of size 3) [km/s],
!               defined relative to galactic rest frame.
!   vobs        Observer/detector's speed (i.e. Sun's speed) [km/s],
!               defined relative to MB rest frame.
!   v0          Most probable speed [km/s] in the MB rest frame.
!   vesc        Galactic/population escape speed [km/s] in the MB rest
!               frame (_galactic_ escape speed only if MB has no bulk
!               motion relative to galactic rest frame).
! Optional output arguments regarding tabulated g(vmin) and h(vmin):
!   tabulated   Indicates if a tabulated g(vmin) is being used
!   g_file      The file tabulated g and h were taken from
!   Nvmin       Number of tabulation points
!   vmin        Allocatable array of vmin [km/s]
!   g_vmin      Allocatable array of mean inverse speeds at vmin [s/km]
!   h_vmin      Allocatable array of mean speeds at vmin [km/s]
! 
SUBROUTINE GetHalo(Halo,vrot,vlsr,vpec,vsun,rho,vbulk,vobs,v0,vesc,    &
                   tabulated,g_file,h_file,Nvmin,vmin,g_vmin,h_vmin)
  IMPLICIT NONE
  TYPE(HaloStruct), INTENT(IN) :: Halo
  CHARACTER(LEN=*), INTENT(OUT), OPTIONAL :: g_file, h_file
  LOGICAL, INTENT(OUT), OPTIONAL :: tabulated
  INTEGER, INTENT(OUT), OPTIONAL :: Nvmin
  REAL*8, INTENT(OUT), OPTIONAL :: vrot,vlsr(3),vpec(3),vsun(3),       &
                                   rho,vbulk(3),vobs,v0,vesc
  REAL*8, ALLOCATABLE, INTENT(OUT), OPTIONAL :: vmin(:),g_vmin(:),h_vmin(:)
  
  IF (PRESENT(vrot))  vrot  = Halo%vrot
  IF (PRESENT(vlsr))  vlsr  = Halo%vlsr
  IF (PRESENT(vpec))  vpec  = Halo%vpec
  IF (PRESENT(vsun))  vsun  = Halo%vsun
  
  IF (PRESENT(rho))   rho   = Halo%rho
  
  IF (PRESENT(vbulk)) vbulk = Halo%vbulk
  IF (PRESENT(vobs))  vobs  = Halo%vobs
  IF (PRESENT(v0))    v0    = Halo%v0
  IF (PRESENT(vesc))  vesc  = Halo%vesc
  
  IF (PRESENT(tabulated)) tabulated = Halo%tabulated
  IF (PRESENT(g_file)) g_file = Halo%g_file
  IF (PRESENT(Nvmin)) Nvmin = Halo%Nvmin
  IF (PRESENT(vmin)) THEN
    ALLOCATE(vmin(Halo%Nvmin))
    vmin = Halo%vmin
  END IF
  IF (PRESENT(g_vmin)) THEN
    ALLOCATE(g_vmin(Halo%Nvmin))
    g_vmin = Halo%g_vmin
  END IF
  IF (PRESENT(h_vmin)) THEN
    ALLOCATE(h_vmin(Halo%Nvmin))
    h_vmin = Halo%h_vmin
  END IF
END SUBROUTINE


! ----------------------------------------------------------------------
! Set various halo quantities.
! 
! Optional input arguments regarding galactic motions:
!   vrot        Local galactic disk rotation speed [km/s].
!   vlsr        Local standard of rest velocity vector (array of size 3)
!               [km/s], defined relative to galactic rest frame.
!   vpec        Sun's peculiar velocity vector (array of size 3) [km/s],
!               defined relative to local standard of rest.
!   vsun        Sun's velocity vector (array of size 3) [km/s], defined
!               relative to galactic rest frame.
! Optional input arguments regarding dark matter density:
!   rho         Local dark matter density [GeV/cm^3].
! Optional input arguments regarding SHM distribution, a truncated
! Maxwell-Boltzmann ("MB"):
!   vbulk       Bulk velocity of dark matter (array of size 3) [km/s],
!               defined relative to galactic rest frame.
!   vobs        Observer/detector's speed (i.e. Sun's speed) [km/s],
!               defined relative to MB rest frame.
!   v0          Most probable speed [km/s] in the MB rest frame.
!   vesc        Galactic/population escape speed [km/s] in the MB rest
!               frame (_galactic_ escape speed only if MB has no bulk
!               motion relative to galactic rest frame).
! Optional tabulated g(vmin) and h(vmin) arguments.  Can be loaded from a 
! given file or explicitly provided.  If provided, Nvmin, vmin and g must
! all be given to take effect. When a tabulation is not provided, the
! mean inverse speed will be calculated explicitly (not tabulated!)
! using the SHM as described by the above parameters. If h (or h_column) is
! not provided, the corresponding entries will be set to zero.
!   tabulated   Indicates if a tabulated g(vmin) is to be used.  Implied
!               by the use of other tabulation arguments, but can be set
!               false to return to the SHM calculation after a tabulation
!               has been loaded.
!   g_file      File from which tabulated g(vmin) should be read
!   g_column    The column in the file to take g from (default is 2)
!   h_column    The column in the file to take h from. If not specified,
!               h will be set to zero.
!   Nvmin       Number of tabulation points
!   vmin        Array of size [1:Nvmin] containing tabulation vmin [km/s]
!   g_vmin      Array of size [1:Nvmin] containing tabulated mean inverse
!               speeds at vmin [s/km]
!   h_vmin      Array of size [1:Nvmin] containing tabulated mean
!               speeds at vmin [km/s]
! 
SUBROUTINE SetHalo(Halo,vrot,vlsr,vpec,vsun,vobs,rho,vbulk,v0,vesc, &
            tabulated,g_file,g_column,h_column,Nvmin,vmin,g_vmin,h_vmin)
  IMPLICIT NONE
  TYPE(HaloStruct), INTENT(INOUT) :: Halo
  CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: g_file
  LOGICAL, INTENT(IN), OPTIONAL :: tabulated
  INTEGER, INTENT(IN), OPTIONAL :: g_column,h_column,Nvmin
  REAL*8, INTENT(IN), OPTIONAL :: vrot,vlsr(3),vpec(3),vsun(3),         &
                                  rho,vbulk(3),vobs,v0,vesc
  REAL*8, INTENT(IN), OPTIONAL :: vmin(:),g_vmin(:),h_vmin(:)
  INTEGER :: K
  
  IF (PRESENT(vrot))  CALL SetDiskRotationSpeed(vrot,Halo)
  IF (PRESENT(vlsr))  CALL SetLocalStandardOfRest(vlsr,Halo)
  IF (PRESENT(vpec))  CALL SetSunPeculiarVelocity(vpec,Halo)
  IF (PRESENT(vsun))  CALL SetSunVelocity(vsun,Halo)
  
  IF (PRESENT(rho))   CALL SetLocalDensity(rho,Halo)
  
  IF (PRESENT(vbulk)) CALL SetBulkVelocity(vbulk,Halo)
  IF (PRESENT(vobs))  CALL SetObserverSpeed(vobs,Halo)
  IF (PRESENT(v0))    CALL SetMostProbableSpeed(v0,Halo)
  IF (PRESENT(vesc))  CALL SetEscapeSpeed(vesc,Halo)
  
  IF (PRESENT(tabulated)) THEN
    IF (Halo%Nvmin .GT. 0) Halo%tabulated = tabulated
  END IF
  IF (PRESENT(Nvmin) .AND. PRESENT(vmin) .AND. PRESENT(g_vmin)) THEN
    IF (Nvmin .GT. 0) THEN
      IF (ALLOCATED(Halo%vmin)) DEALLOCATE(Halo%vmin)
      IF (ALLOCATED(Halo%g_vmin))  DEALLOCATE(Halo%g_vmin)
      IF (ALLOCATED(Halo%h_vmin))  DEALLOCATE(Halo%h_vmin)
      Halo%Nvmin = Nvmin
      Halo%vmin  = vmin
      Halo%g_vmin = g_vmin
      IF (PRESENT(h_vmin)) THEN
        Halo%h_vmin = h_vmin
      ELSE
        ALLOCATE(Halo%h_vmin(Nvmin))
        Halo%h_vmin(:) = 0
      END IF
      Halo%tabulated = .TRUE.
      Halo%g_file  = ''
    END IF
  END IF
  IF (PRESENT(g_file)) THEN
    IF (PRESENT(g_column)) THEN
      K = g_column
    ELSE
      K = 2
    END IF
    CALL LoadArrays(file=g_file,N=Halo%Nvmin,N1=1,C1=Halo%vmin,       &
                    N2=K,C2=Halo%g_vmin)
    IF (PRESENT(h_column)) THEN
      CALL LoadArrays(file=g_file,N=K,N1=h_column,C1=Halo%h_vmin)
    ELSE
        ALLOCATE(Halo%h_vmin(Halo%Nvmin))
        Halo%h_vmin(:) = 0      
    END IF
    Halo%tabulated = .TRUE.
    Halo%g_file  = g_file
  END IF
  
END SUBROUTINE


!-----------------------------------------------------------------------
! Initializes halo.
! Simply sets halo parameters to default values.
! 
FUNCTION InitHalo() RESULT(Halo)

  IMPLICIT NONE
  TYPE(HaloStruct) :: Halo
  
  Halo%vrot  = 235d0
  Halo%vlsr  = (/ 0d0, 235d0, 0d0 /)
  Halo%vpec  = (/ 11d0, 12d0, 7d0 /)
  Halo%vsun  = (/ 0d0, 235d0, 0d0 /) + (/ 11d0, 12d0, 7d0 /)
  
  Halo%rho   = 0.4d0
  
  Halo%vbulk = (/ 0d0, 0d0, 0d0 /)
  Halo%vobs  = SQRT(11d0**2 + (235d0+12d0)**2 + 7d0**2)
  Halo%v0    = 235d0
  Halo%vesc  = 550d0
  
  Halo%tabulated = .FALSE.
  Halo%g_file  = ''
  Halo%Nvmin = 0
  IF (ALLOCATED(Halo%vmin)) DEALLOCATE(Halo%vmin)
  IF (ALLOCATED(Halo%g_vmin))  DEALLOCATE(Halo%g_vmin)
  IF (ALLOCATED(Halo%h_vmin))  DEALLOCATE(Halo%h_vmin)
  
END FUNCTION

!-----------------------------------------------------------------------
! C/C++ wrapper for InitHalo.
!
INTEGER(KIND=C_INT) FUNCTION C_DDCalc_InitHalo() &
 BIND(C,NAME='C_DDHalo_ddcalc_inithalo') 
  USE ISO_C_BINDING, only: C_INT
  IMPLICIT NONE
  N_Halos = N_Halos + 1
  IF (N_Halos .GT. Max_Halos) stop 'DDCalc: Max_Halos exceeded.&
   Please run FreeHalos or modify Max_Halos in DDTypes.f90.'
  ALLOCATE(Halos(N_Halos)%p)
  Halos(N_Halos)%p = InitHalo()
  C_DDCalc_InitHalo = N_Halos
END FUNCTION

! ----------------------------------------------------------------------
! Get/set local galactic disk rotation speed [km/s].
! Modifies the Local Standard of Rest (LSR) and Sun's velocity relative
! to halo rest frame as well as the most probable speed of the velocity
! distribution (v0 = vrot).  The observer speed is updated.
! 
PURE FUNCTION GetDiskRotationSpeed(Halo) RESULT(vrot)
  IMPLICIT NONE
  TYPE(HaloStruct), INTENT(IN) :: Halo
  REAL*8 :: vrot
  vrot = Halo%vrot
END FUNCTION

SUBROUTINE SetDiskRotationSpeed(vrot,Halo)
  IMPLICIT NONE
  TYPE(HaloStruct), INTENT(INOUT) :: Halo
  REAL*8, INTENT(IN) :: vrot
  Halo%vrot = vrot
  Halo%vlsr = (/ 0d0, Halo%vrot, 0d0 /)
  Halo%vsun = Halo%vlsr + Halo%vpec
  Halo%vobs = SQRT(SUM((Halo%vsun - Halo%vbulk)**2))
  Halo%v0   = Halo%vrot
END SUBROUTINE


! ----------------------------------------------------------------------
! Get/set Local Standard Of Rest velocity vector [km/s], defined
! relative to galactic rest frame.  Usually assumed to be (0,vrot,0),
! where vrot is disk rotation speed.  Modifies Sun's velocity relative
! to halo rest frame.  The disk rotation speed and the most probable
! speed of the velocity distribution are set to the y component of this
! velocity vector.  The observer speed is updated.
! 
PURE FUNCTION GetLocalStandardOfRest(Halo) RESULT(vlsr)
  IMPLICIT NONE
  TYPE(HaloStruct), INTENT(IN) :: Halo
  REAL*8 :: vlsr(3)
  vlsr = Halo%vlsr
END FUNCTION

SUBROUTINE SetLocalStandardOfRest(vlsr,Halo)
  IMPLICIT NONE
  TYPE(HaloStruct), INTENT(INOUT) :: Halo
  REAL*8, INTENT(IN) :: vlsr(3)
  Halo%vlsr = vlsr
  Halo%vsun = Halo%vlsr + Halo%vpec
  Halo%vrot = vlsr(2)
  Halo%vobs = SQRT(SUM((Halo%vsun - Halo%vbulk)**2))
  Halo%v0   = ABS(Halo%vrot)
END SUBROUTINE


! ----------------------------------------------------------------------
! Get/set Sun's peculiar velocity vector [km/s], defined relative to
! local standard of rest.  Modifies Sun's velocity relative to halo
! rest frame.  The observer speed is updated.
! 
PURE FUNCTION GetSunPeculiarVelocity(Halo) RESULT(vpec)
  IMPLICIT NONE
  TYPE(HaloStruct), INTENT(IN) :: Halo
  REAL*8 :: vpec(3)
  vpec = Halo%vpec
END FUNCTION

SUBROUTINE SetSunPeculiarVelocity(vpec,Halo)
  IMPLICIT NONE
  TYPE(HaloStruct), INTENT(INOUT) :: Halo
  REAL*8, INTENT(IN) :: vpec(3)
  Halo%vpec = vpec
  Halo%vsun = Halo%vlsr + Halo%vpec
  Halo%vobs = SQRT(SUM((Halo%vsun - Halo%vbulk)**2))
END SUBROUTINE


! ----------------------------------------------------------------------
! Get/set Sun's velocity vector [km/s], defined relative to galactic
! rest frame.  Normally taken to be vlsr + vpec, i.e. the sum of the
! local standard of rest and the Sun's peculiar velocity.  The preferred
! way to set speeds is modifying the disk rotation speed or Sun's
! peculiar velocity, not by setting this velocity directly, as the
! contributing velocities become ill-defined.  If the Sun's velocity is
! set here, the routine will attempt to assign a rotation speed vrot
! and local standard of rest vlsr = (0,vrot,0) that matches the given
! velocity vector, using the current value of the peculiar velocity; if
! not possible, the peculiar motion is set to zero first.  The most
! probable speed of the velocity distribution is updated to the
! resulting vrot and the observer speed is updated.
! 
PURE FUNCTION GetSunVelocity(Halo) RESULT(vsun)
  IMPLICIT NONE
  TYPE(HaloStruct), INTENT(IN) :: Halo
  REAL*8 :: vsun(3)
  vsun = Halo%vsun
END FUNCTION

SUBROUTINE SetSunVelocity(vsun,Halo)
  IMPLICIT NONE
  TYPE(HaloStruct), INTENT(INOUT) :: Halo
  REAL*8, INTENT(IN) :: vsun(3)
  REAL*8 :: vrot2
  Halo%vsun = vsun
  vrot2 = SUM((Halo%vsun - Halo%vpec)**2)
  IF (vrot2 .GE. 0d0) THEN
    Halo%vrot = SQRT(vrot2)
    Halo%vlsr = (/ 0d0, Halo%vrot, 0d0 /)
  ELSE
    Halo%vpec = 0d0
    Halo%vrot = SQRT(SUM(Halo%vsun**2))
    Halo%vlsr = (/ 0d0, Halo%vrot, 0d0 /)
  END IF
  Halo%v0 = ABS(Halo%vrot)
  Halo%vobs = SQRT(SUM((Halo%vsun - Halo%vbulk)**2))
END SUBROUTINE


! ----------------------------------------------------------------------
! Get/set local halo density [GeV/cm^3].
! 
PURE FUNCTION GetLocalDensity(Halo) RESULT(rho)
  IMPLICIT NONE
  TYPE(HaloStruct), INTENT(IN) :: Halo
  REAL*8 :: rho
  rho = Halo%rho
END FUNCTION

SUBROUTINE SetLocalDensity(rho,Halo)
  IMPLICIT NONE
  TYPE(HaloStruct), INTENT(INOUT) :: Halo
  REAL*8, INTENT(IN) :: rho
  Halo%rho = MAX(rho,0d0)
END SUBROUTINE


! ----------------------------------------------------------------------
! Get/set the dark matter population's bulk velocity vector [km/s],
! defined relative to the galactic rest frame.  Modifies the observer
! speed.
! 
PURE FUNCTION GetBulkVelocity(Halo) RESULT(vbulk)
  IMPLICIT NONE
  TYPE(HaloStruct), INTENT(IN) :: Halo
  REAL*8 :: vbulk(3)
  vbulk = Halo%vbulk
END FUNCTION

SUBROUTINE SetBulkVelocity(vbulk,Halo)
  IMPLICIT NONE
  TYPE(HaloStruct), INTENT(INOUT) :: Halo
  REAL*8, INTENT(IN) :: vbulk(3)
  Halo%vbulk = vbulk
  Halo%vobs = SQRT(SUM((Halo%vsun - Halo%vbulk)**2))
END SUBROUTINE


! ----------------------------------------------------------------------
! Get/set observer/detector's speed (i.e. Sun's speed) [km/s], defined
! relative to Maxwell-Boltzmann population rest frame.  Normally taken
! to be |vlsr + vpec - vMB|, i.e. the sum of the local standard of rest
! and the Sun's peculiar velocity less the bulk velocity of the dark
! matter population.  The preferred way to set speeds is modifying the
! disk rotation speed, Sun's peculiar velocity, or the bulk dark matter
! motion, not by setting this speed directly, as the various
! velocities become ill-defined.  If the observer's speed is set here,
! the routine will set the bulk motion of the DM to zero (relative to
! the galactic rest frame) and attempt to assign a rotation speed vrot
! and local standard of rest vlsr = (0,vrot,0) that matches the given
! speed, using the current value of the peculiar velocity; if not
! possible, the peculiar motion is set to zero first.  The most
! probable speed of the velocity distribution is updated to the
! resulting vrot.
! 
PURE FUNCTION GetObserverSpeed(Halo) RESULT(vobs)
  IMPLICIT NONE
  TYPE(HaloStruct), INTENT(IN) :: Halo
  REAL*8 :: vobs
  vobs = Halo%vobs
END FUNCTION

SUBROUTINE SetObserverSpeed(vobs,Halo)
  IMPLICIT NONE
  TYPE(HaloStruct), INTENT(INOUT) :: Halo
  REAL*8, INTENT(IN) :: vobs
  REAL*8 :: vy2
  Halo%vobs  = MAX(vobs,0d0)
  Halo%vbulk = (/ 0d0, 0d0, 0d0 /)
  vy2 = Halo%vobs**2 - Halo%vpec(1)**2 - Halo%vpec(3)**2
  IF (vy2 .GE. 0d0) THEN
    Halo%vrot = SQRT(vy2) - Halo%vpec(2)
    Halo%vlsr = (/ 0d0, Halo%vrot, 0d0 /)
    Halo%vsun = Halo%vlsr + Halo%vpec
  ELSE
    Halo%vpec = 0d0
    Halo%vrot = Halo%vobs
    Halo%vlsr = (/ 0d0, Halo%vrot, 0d0 /)
    Halo%vsun = Halo%vlsr + Halo%vpec
  END IF
  Halo%v0 = ABS(Halo%vrot)
END SUBROUTINE


! ----------------------------------------------------------------------
! Get/set most probable speed v0 [km/s] in the dark matter population's
! rest frame. Related to other speeds characterizing velocity
! distribution by:
!     vrms = sqrt(3/2) v0    [rms velocity]
!     vmp  = v0              [most probably velocity]
!     vave = sqrt(4/pi) v0   [mean velocity]
! 
PURE FUNCTION GetMostProbableSpeed(Halo) RESULT(v0)
  IMPLICIT NONE
  TYPE(HaloStruct), INTENT(IN) :: Halo
  REAL*8 :: v0
  v0 = Halo%v0
END FUNCTION

SUBROUTINE SetMostProbableSpeed(v0,Halo)
  IMPLICIT NONE
  TYPE(HaloStruct), INTENT(INOUT) :: Halo
  REAL*8, INTENT(IN) :: v0
  Halo%v0 = v0
END SUBROUTINE


! ----------------------------------------------------------------------
! Get/set dark matter population escape speed [km/s].  In the case of
! the SHM with no bulk motion relative to the galactic rest frame, this
! is the galactic escape speed.
! 
PURE FUNCTION GetEscapeSpeed(Halo) RESULT(vesc)
  IMPLICIT NONE
  TYPE(HaloStruct), INTENT(IN) :: Halo
  REAL*8 :: vesc
  vesc = Halo%vesc
END FUNCTION

SUBROUTINE SetEscapeSpeed(vesc,Halo)
  IMPLICIT NONE
  TYPE(HaloStruct), INTENT(INOUT) :: Halo
  REAL*8, INTENT(IN) :: vesc
  Halo%vesc = MAX(vesc,0d0)
END SUBROUTINE

! ----------------------------------------------------------------------
! INTERFACE NAME: EToVmin
! Calculate the minimum velocity for producing a recoil of energy E,
! given by vmin = sqrt{M E/(2\mu^2)} [km/s].  Returns as array of
! size [1:N,1:Niso].
! 
! This is the 2D array version (multiple masses and array of energies).
! 
! Input arguments:
!   N           Number of recoil energies
!   E           Array of recoil energies [keV]
!   m           WIMP mass [GeV]
!   Niso        Number of isotopes
!   Miso        Array of isotope masses [GeV]
! 
PURE FUNCTION EToVmin(N,E,m,Niso,Miso) RESULT(vmin)
  IMPLICIT NONE
  REAL*8 :: vmin(N,Niso)
  INTEGER, INTENT(IN) :: N,Niso
  REAL*8, INTENT(IN) :: E(N),m,Miso(Niso)
  INTEGER :: I
  REAL*8 :: mu(Niso)
  REAL*8, PARAMETER :: c = 1d-3*SPEED_OF_LIGHT  ! Speed of light in km/s
  mu = Miso*m / (Miso + m)
  DO I = 1,Niso
    vmin(:,I) = c * SQRT(1d-6*Miso(I)*E/(2*mu(I)**2))
  END DO
END FUNCTION


!-----------------------------------------------------------------------
! INTERFACE NAME: MeanInverseSpeed
! Calculates the mean inverse speed (g) [s/km] for the given 2D
! array of vmin, with g defined as:
!     g(vmin) = \int_{|v|>vmin} d^3v 1/|v| f(v)
! Returns as array of size [1:N1,1:N2].
! 
! This is the 2D array version (2D array of vmin).
! 
! Input arguments:
!   N1,N2       Size of vmin and g arrays, i.e. [1:N1,1:N2]
!   vmin        The minimum speed in the g integral [km/s].
! 
PURE FUNCTION MeanInverseSpeed(N1,N2,vmin,Halo) RESULT(g_vmin)
  IMPLICIT NONE
  TYPE(HaloStruct), INTENT(IN) :: Halo
  REAL*8 :: g_vmin(N1,N2)
  INTEGER, INTENT(IN) :: N1,N2
  REAL*8, INTENT(IN) :: vmin(N1,N2)
  REAL*8 :: v0,vobs,vesc,x(N1,N2),y,z,Nesc
  
  ! If have tabulation, use it
  IF (Halo%tabulated) THEN
    g_vmin = MeanInverseSpeedT(vmin,Halo)
    RETURN
  END IF
  
  ! Easier to use variable names
  v0   = Halo%v0
  vobs = Halo%vobs
  vesc = Halo%vesc
  
  ! Special case: no dispersion
  ! Distribution is delta function
  IF (v0 .EQ. 0) THEN
    IF (vobs .EQ. 0d0) THEN
      g_vmin = 0d0
    ELSE
      WHERE (vmin .LE. vobs)
        g_vmin = 1d0 / vobs
      ELSE WHERE
        g_vmin = 0d0
      END WHERE
    END IF
    RETURN
  END IF
  
  x    = vmin / v0
  y    = vobs / v0
  z    = vesc / v0
  Nesc = ERF(z) - 2*INVSQRTPI*z*EXP(-z**2)
  
  ! Special case: no relative motion by observer
  !   g_vmin = 2/(sqrt(pi) Nesc v0) [e^{-x^2} - e^{-z^2}]
  ! Note: EXP2(a,b) = e^b - e^a
  IF (y .EQ. 0d0) THEN
    WHERE (x .LE. z)
      g_vmin = 2*INVSQRTPI/(Nesc*v0) * EXP2(-z**2,-x**2)
    ELSE WHERE
      g_vmin = 0d0
    END WHERE
    RETURN
  END IF
  
  ! Special case: no finite cutoff (vesc is effectively infinite)
  IF (z .GT. 25d0) THEN
    g_vmin = ERF2(x-y,x+y) / (2*vobs)
    RETURN
  END IF
  
  ! General case.
  ! See e.g. Savage, Freese & Gondolo, PRD 74, 043531 (2006)
  ! [astrop-ph/0607121]; use arxiv version as PRD version has type-
  ! setting issues in the formula.
  ! Note: ERF2(a,b) = ERF(b) - ERF(a)
  ! Separate y < z & y > z cases to make easier use of WHERE statements.
  IF (y .LT. z) THEN
    WHERE (x .LT. z-y)
      g_vmin = 1d0 / (2*Nesc*vobs) * (ERF2(x-y,x+y) - 4*INVSQRTPI*y*EXP(-z**2))
    ELSE WHERE (x .LT. z+y)
      g_vmin = 1d0 / (2*Nesc*vobs) * (ERF2(x-y,z) - 2*INVSQRTPI*(z+y-x)*EXP(-z**2))
    ELSE WHERE
      g_vmin = 0d0
    END WHERE
  ELSE
    WHERE (x .LT. y-z)
      g_vmin = 1d0 / vobs
    ELSE WHERE (x .LT. y+z)
      g_vmin = 1d0 / (2*Nesc*vobs) * (ERF2(x-y,z) - 2*INVSQRTPI*(z+y-x)*EXP(-z**2))
    ELSE WHERE
      g_vmin = 0d0
    END WHERE
  END IF
  
END FUNCTION


!-----------------------------------------------------------------------
! Calculates the mean inverse speed (g) [s/km] for the given vmin,
! with g_vmin define as:
!     g(vmin) = \int_{|v|>vmin} d^3v 1/|v| f(v)
! using the stored tabulation rather than the explicit calculation.
! 
! Input arguments:
!   vmin        The minimum speed in the g integral [km/s]
! 
ELEMENTAL FUNCTION MeanInverseSpeedT(vmin,Halo) RESULT(g_vmin)
  IMPLICIT NONE
  TYPE(HaloStruct), INTENT(IN) :: Halo
  REAL*8 :: g_vmin
  REAL*8, INTENT(IN) :: vmin
  INTEGER :: K
  REAL*8 :: f
  
  IF (.NOT. Halo%tabulated .OR. (Halo%Nvmin .LE. 0)) THEN
    g_vmin = 0d0
    RETURN
  END IF
  
  K = BSearch(Halo%Nvmin,Halo%vmin,vmin)
  
  IF (K .LE. 0) THEN
    g_vmin = Halo%g_vmin(1)
  ELSE IF (K .GE. Halo%Nvmin) THEN
    IF (vmin .EQ. Halo%vmin(Halo%Nvmin)) THEN
      g_vmin = Halo%g_vmin(Halo%Nvmin)
    ELSE
      g_vmin = 0d0
    END IF
  ELSE IF (Halo%vmin(K) .EQ. Halo%vmin(K+1)) THEN
    g_vmin = Halo%g_vmin(K)
  ELSE
    f = (vmin-Halo%vmin(K)) / (Halo%vmin(K+1)-Halo%vmin(K))
    g_vmin = (1-f)*Halo%g_vmin(K) + f*Halo%g_vmin(K+1)
  END IF
  
END FUNCTION

!-----------------------------------------------------------------------
! INTERFACE NAME: MeanSpeed
! Calculates the mean speed (h) [km/s] for the given 2D
! array of vmin, with h defined as:
!     h(vmin) = \int_{|v|>vmin} d^3v |v| f(v)
! Returns as array of size [1:N1,1:N2].
! 
! This is the 2D array version (2D array of vmin).
! 
! Input arguments:
!   N1,N2       Size of vmin and h arrays, i.e. [1:N1,1:N2]
!   vmin        The minimum speed in the h integral [km/s].
! 
FUNCTION MeanSpeed(N1,N2,vmin,Halo) RESULT(h_vmin)
  IMPLICIT NONE
  TYPE(HaloStruct), INTENT(IN) :: Halo
  REAL*8 :: h_vmin(N1,N2)
  INTEGER, INTENT(IN) :: N1,N2
  REAL*8, INTENT(IN) :: vmin(N1,N2)
  REAL*8 :: v0,vobs,vesc,x(N1,N2),y,z,Nesc
  
  ! If have tabulation, use it
  IF (Halo%tabulated) THEN
    h_vmin = MeanSpeedT(vmin,Halo)
    RETURN
  END IF

  ! Easier to use variable names
  v0   = Halo%v0
  vobs = Halo%vobs
  vesc = Halo%vesc

  IF ((v0 .EQ. 0) .OR. (vobs .EQ. 0) .OR. (vesc .EQ. 0)) THEN
    WRITE(*,*) 'ERROR: One of more halo parameters are zero. Cannot calculate mean velocity.'
    STOP   
  END IF

  x    = vmin / v0
  y    = vobs / v0
  z    = vesc / v0
  Nesc = ERF(z) - 2*INVSQRTPI*z*EXP(-z**2)

  IF (z .GT. 25d0) THEN
    WRITE(*,*) 'ERROR: Escape velocity too large. Cannot calculate mean velocity.'
    STOP
  END IF

  IF (z .LT. y) THEN
    WRITE(*,*) 'ERROR: Escape velocity too small. Cannot calculate mean velocity.'
    STOP
  END IF

  WHERE (x .LT. z-y)
    h_vmin = 1.0/Nesc * v0 * (((x-y)/(2*y*SQRTPI) + INVSQRTPI) * exp(-(x-y)**2) & 
        - ((x+y)/(2*y*SQRTPI) - INVSQRTPI) * exp(-(x+y)**2) + (1 + 2 * y**2) / (4*y) * (erf((x+y)) - erf((x-y))) &
        - INVSQRTPI * ( 2 + ( (x+z-(x-y))**3 - (x+z-(x+y))**3 )/(3*y) ) * exp(-z**2) )
  ELSE WHERE (x .LT. z+y)
    h_vmin = 1.0/Nesc * v0 * (((x-y)/(2*y*SQRTPI) + INVSQRTPI) * exp(-(x-y)**2) &
        - (z/(2*y*SQRTPI) - INVSQRTPI) * exp(-z**2) + (1 + 2 * y**2) / (4*y) * (erf(z) - erf((x-y))) &
        - INVSQRTPI * ( 2 + ( (x+z-(x-y))**3 - (x+z-z)**3 )/(3*y) ) * exp(-z**2) )
  ELSE WHERE
    h_vmin = 0d0
  END WHERE
  
END FUNCTION

!-----------------------------------------------------------------------
! Calculates the mean speed (h) [km/s] for the given vmin,
! with h defined as:
!     h(vmin) = \int_{|v|>vmin} d^3v |v| f(v)
! using the stored tabulation rather than the explicit calculation.
! 
! Input arguments:
!   vmin        The minimum speed in the h integral [km/s]
! 
ELEMENTAL FUNCTION MeanSpeedT(vmin,Halo) RESULT(h_vmin)
  IMPLICIT NONE
  TYPE(HaloStruct), INTENT(IN) :: Halo
  REAL*8 :: h_vmin
  REAL*8, INTENT(IN) :: vmin
  INTEGER :: K
  REAL*8 :: f
  
  IF (.NOT. Halo%tabulated .OR. (Halo%Nvmin .LE. 0)) THEN
    h_vmin = 0d0
    RETURN
  END IF
  
  K = BSearch(Halo%Nvmin,Halo%vmin,vmin)
  
  IF (K .LE. 0) THEN
    h_vmin = Halo%h_vmin(1)
  ELSE IF (K .GE. Halo%Nvmin) THEN
    IF (vmin .EQ. Halo%vmin(Halo%Nvmin)) THEN
      h_vmin = Halo%h_vmin(Halo%Nvmin)
    ELSE
      h_vmin = 0d0
    END IF
  ELSE IF (Halo%vmin(K) .EQ. Halo%vmin(K+1)) THEN
    h_vmin = Halo%h_vmin(K)
  ELSE
    f = (vmin-Halo%vmin(K)) / (Halo%vmin(K+1)-Halo%vmin(K))
    h_vmin = (1-f)*Halo%h_vmin(K) + f*Halo%h_vmin(K+1)
  END IF
  
END FUNCTION



END MODULE
