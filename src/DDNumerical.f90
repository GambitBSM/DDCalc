MODULE DDNumerical

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! DDNumerical
!    Numerical routines for use in DDCalc. Includes interpolation, random
!    number generation, special functions and probability distributions.
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

USE DDUtils

IMPLICIT NONE
PRIVATE

PUBLIC :: LOG_SUM, LOGp1, GAMMAI, PoissonLogPDF, ERF2, EXP2

!Interface for interpolation routines.
INTERFACE LinearInterpolate
  MODULE PROCEDURE LinearInterpolate_S,LinearInterpolate_A
END INTERFACE

! The structure here stores state data for the uniform random number
! generator.
TYPE, PRIVATE :: RandState
  LOGICAL :: initialized = .FALSE.
  INTEGER :: I,J      ! I97,J97
  REAL*8 ::  U(1:97)
  REAL*8 ::  C,CD,CM
END TYPE
! Instantiation for default
! OpenMP directive ensures each thread maintains own copy
TYPE(RandState), PRIVATE :: DEFAULT_RandState
!$OMP THREADPRIVATE(DEFAULT_RandState)

CONTAINS


!=======================================================================
! INTERPOLATION
! Interpolation routines are given in two forms: scalar and array
! interpolation location(s).
!=======================================================================

! ----------------------------------------------------------------------
! Determines y(x=x0) using linear interpolation between the points
! specified by the arrays x(N) and y(N).
! 
! Input arguments:
!   N               Length of x & y arrays
!   x,y             Arrays of x and y points defining y(x).
!                   Array in x must be increasing.
!   x0              Value to interpolate at
! Optional input arguments:
!   extrapolate     Indicate if y(x0) should be extrapolated from the
!                   given points if x0 falls outside of the domain of x.
!                   If set to .FALSE., this routine will return 0 if
!                   x0 is outside the domain.  Default is .TRUE.
! 
PURE FUNCTION LinearInterpolate_S(N,x,y,x0,extrapolate) RESULT(y0)
  IMPLICIT NONE
  REAL*8 :: y0
  INTEGER, INTENT(IN) :: N
  REAL*8, INTENT(IN) :: x(N),y(N),x0
  LOGICAL, INTENT(IN), OPTIONAL :: extrapolate
  INTEGER :: I
  REAL*8 :: a,b
  
  y0 = 0d0
  IF (N .EQ. 0) RETURN
  
  IF (PRESENT(extrapolate)) THEN
    IF ((.NOT. extrapolate)                                             &
        .AND. ((x0 .LT. x(1)) .OR. (x0 .GT. x(N)))) THEN
      RETURN
    END IF
  END IF
  
  IF (N .EQ. 1) THEN
    y0 = y(1)
    RETURN
  END IF
  
  ! Find interval to use for linear interpolation.
  ! Index here is for lower point in interval.
  I = BSearch(N,x,x0)
  I = MIN(I,N-1)
  I = MAX(I,1)
  
  IF (x(I) .EQ. x(I+1)) THEN
    y0 = 0.5d0 * (y(I) + y(I+1))
  ELSE
    a  = (y(I+1) - y(I)) / (x(I+1) - x(I))
    b  = y(I) - a*x(I)
    y0 = a*x0 + b
  END IF
  
END FUNCTION


! ----------------------------------------------------------------------
! Determines y(x=x0) using linear interpolation between the points
! specified by the arrays x(N) and y(N).
! 
! Input arguments:
!   N               Length of x & y arrays
!   x,y             Arrays of x and y points defining y(x).
!                   Array in x must be increasing.
!   N0              Length of x0 array
!   x0              Array of values to interpolate at
! Optional input arguments:
!   extrapolate     Indicate if y(x0) should be extrapolated from the
!                   given points if x0 falls outside of the domain of x.
!                   If set to .FALSE., this routine will return 0 if
!                   x0 is outside the domain.  Default is .TRUE.
! 
PURE FUNCTION LinearInterpolate_A(N,x,y,N0,x0,extrapolate) RESULT(y0)
  IMPLICIT NONE
  REAL*8 :: y0(N0)
  INTEGER, INTENT(IN) :: N,N0
  REAL*8, INTENT(IN) :: x(N),y(N),x0(N0)
  LOGICAL, INTENT(IN), OPTIONAL :: extrapolate
  LOGICAL :: extrapolate0
  INTEGER :: I,I0,Is
  REAL*8 :: a,b
  
  extrapolate0 = .TRUE.
  IF (PRESENT(extrapolate)) extrapolate0 = extrapolate
  
  ! Nothing to interpolate
  IF (N .EQ. 0) THEN
    y0 = 0d0
    RETURN
  END IF
  
  ! Single point: no searching necessary
  IF (N .EQ. 1) THEN
    IF (extrapolate0) THEN
      y0 = y(1)
    ELSE
      WHERE(x0 .EQ. x(1))
        y0 = y(1)
      ELSE WHERE
        y0 = 0d0
      END WHERE
    END IF
    RETURN
  END IF
  
  ! Cycle through x0 points
  Is = 0
  IF (extrapolate0) THEN
    DO I0 = 1,N0
      Is = BSearch(N,x,x0(I0),Is)
      Is = MAX(MIN(I,N-1),1)
      IF (x(I0) .EQ. x(I0+1)) THEN
        y0 = 0.5d0 * (y(I0) + y(I0+1))
      ELSE
        a      = (y(I0+1) - y(I0)) / (x(I0+1) - x(I0))
        b      = y(I0) - a*x(I0)
        y0(I0) = a*x0(I0) + b
      END IF
    END DO
  ELSE
    DO I0 = 1,N0
      IF ((x0(I0) .LT. x(1)) .OR. (x0(I0) .GE. x(N))) THEN
        y0(I0) = 0d0
      ELSE
        Is = BSearch(N,x,x0(I0),Is)
        IF (x(I0) .EQ. x(I0+1)) THEN
          y0 = 0.5d0 * (y(I0) + y(I0+1))
        ELSE
          a      = (y(I0+1) - y(I0)) / (x(I0+1) - x(I0))
          b      = y(I0) - a*x(I0)
          y0(I0) = a*x0(I0) + b
        END IF
      END IF
    END DO
  END IF
  
END FUNCTION



!=======================================================================
! RANDOM NUMBER GENERATOR
! The generator here is based on the algorithm of Marsaglia & Zaman
! (and later Tsang).  Some of the routines below are based upon an
! implementation by Christian Walck.  Due to the use of an OpenMP
! THREADPRIVATE directive when declaring the internal state structure,
! the various random number routines below should not have any
! performance issues when using multi-threading through OpenMP (as
! opposed to the RANDOM_NUMBER intrinsic, which may lead to _slower_
! code in multi-threaded cases due to its serial-only implementation
! leading to possible thread pile-up).  Upon entering the first OpenMP
! parallel region, each thread is given its own random seed unless
! explicit seeds are set by calling rand_init() within that parallel
! region; per-thread RNG states are then maintained for later parallel
! regions.  Note the RNG state from the serial region may be replaced
! with that of the master thread.
! 
! See:
!   G. Marsaglia and A. Zaman, "Toward a Universal Random Number
!     Generator," Florida State University Report: FSU-SCRI-87-50
!     (1987).
!   G. Marsaglia, A. Zaman and W.W. Tsang, "Toward a universal
!     random number generator," Statistics & Probability Letters,
!     9, 35, (1990).
!   F. James, "A Review of Pseudorandom Number Generators,"
!     Comput. Phys. Commun. 60, 329 (1990).
! 
!=======================================================================


! ----------------------------------------------------------------------
! Initializes state data for the random number generator of Marsaglia,
! Zaman & Tsang, optionally using the given seed.
! 
! Optional input argument:
!     seed       INTEGER seed used to generate state data.
!                If not given, a random seed will be used.
! Optional output argument:
!     state      A RandState structure containing RNG state
!                data that will be initialized.  If not given,
!                the internal RNG state will be initialized.
! 
SUBROUTINE rand_init(seed,state)
  IMPLICIT NONE
  INTEGER, INTENT(IN), OPTIONAL :: seed
  TYPE(RandState), INTENT(OUT), OPTIONAL :: state
  INTEGER :: seed0,I,J,IF1,IF2,IF3,IC1,M
  TYPE(RandState) :: state0
  REAL*8 :: x,S,T
  INTEGER, PARAMETER :: MAX_SEED = 921350144
  
  IF (PRESENT(seed)) THEN
    seed0 = MODULO(seed,MAX_SEED)
  ELSE
    CALL RANDOM_NUMBER(x)
    seed0 = INT(x*MAX_SEED)
  END IF
  
  state0%initialized = .TRUE.
  state0%I = 97
  state0%J = 33
  
  ! Three Fibonacci generator seeds (2-177) and one congruential
  ! generator seed (0-168)
  IF1 = MOD(seed0/169/176/176,176) + 2
  IF2 = MOD(seed0/169/176,176) + 2
  IF3 = MOD(seed0/169,176) + 2
  IC1 = MOD(seed0,169)
  
  DO I = 1,97
     S = 0.0d0
     T = 0.5d0
     DO J = 1,24
        M = MOD(MOD(IF1*IF2,179)*IF3,179)
        IF1 = IF2
        IF2 = IF3
        IF3 = M
        IC1 = MOD(53*IC1+1,169)
        IF ( MOD(IC1*M,64) .GE. 32 ) S = S + T
        T = 0.5d0 * T
     END DO
     state0%U(I) = S
  END DO
  
  state0%C  =   362436.0d0 / 16777216.0d0
  state0%CD =  7654321.0d0 / 16777216.0d0
  state0%CM = 16777213.0d0 / 16777216.0d0
  
  IF (PRESENT(state)) THEN
    state = state0
  ELSE
    DEFAULT_RandState = state0
  END IF
  
END SUBROUTINE


! ----------------------------------------------------------------------
! Returns a unit random number using the algorithm of Marsaglia, Zaman
! & Tsang, optionally using/updating the given state data.
! 
! Optional input/output argument:
!     state      A RandState structure containing RNG state
!                data that will be used for generating the random
!                number; the state will be updated.  If not given,
!                an internal RNG state will be used.
! 
FUNCTION rand(state) RESULT(x)
  IMPLICIT NONE
  REAL*8 :: x
  TYPE(RandState), INTENT(INOUT), OPTIONAL :: state
  IF (PRESENT(state)) THEN
    IF (.NOT. state%initialized) CALL rand_init(state=state)
    CALL rand_number(state,x)
  ELSE
    IF (.NOT. DEFAULT_RandState%initialized) CALL rand_init(state=DEFAULT_RandState)
    CALL rand_number(DEFAULT_RandState,x)
  END IF
END FUNCTION


! ----------------------------------------------------------------------
! Generates a unit random number using the algorithm of Marsaglia, Zaman
! & Tsang.  Requires a state structure to be given, which _must_ be
! initialized.  The rand() function above is the intended to be the
! main routine for external use.
! 
! Required input/output argument:
!     state      A RandState structure containing RNG state
!                data that will be used for generating the random
!                number; the state will be updated.
! Required output argument:
!     x          The uniform random number.
! 
ELEMENTAL SUBROUTINE rand_number(state,x)
  IMPLICIT NONE
  TYPE(RandState), INTENT(INOUT) :: state
  REAL*8, INTENT(OUT) :: x
  LOGICAL :: valid
  
  valid = .FALSE.
  DO WHILE (.NOT. valid)
    x = state%U(state%I) - state%U(state%J)
    IF (x .LT. 0d0) x = x + 1d0
    state%U(state%I) = x
    
    state%I = state%I - 1
    IF (state%I .EQ. 0) state%I = 97
    state%J = state%J - 1
    IF (state%J .EQ. 0) state%J = 97
    
    state%C = state%C - state%CD
    IF (state%C .LT. 0d0) state%C = state%C + state%CM
    x = x - state%C
    IF (x .LT. 0d0) x = x + 1d0
    
    ! Avoid returning zero
    valid = (x .GT. 0d0) .AND. (x .LT. 1d0)
  END DO
  
END SUBROUTINE


!=======================================================================
! SPECIAL FUNCTIONS
!=======================================================================

!-----------------------------------------------------------------------
! Error function between x1 and x2, i.e. erf(x2)-erf(x1).
! This routine accounts for the case when x1 and x2 are similar
! and loss of precision will result from canceling when an explicit
! subtraction of two calls to the implicit ERF function is performed.
! 
! Implementation here appears to be valid to ~30*epsilon, where
! epsilon is the smallest unit of precision (e.g. 1e-15 for double
! precision).  Precision can get somewhat worse than this when
! x1,x2 > 15 or x1,x2 < -15, but only in cases where erf(x1,x2)
! < 1e-100.  The precision has been determined only by testing various
! cases using Mathematica and has not been formally proven.
! 
! Possibly useful identities:
!   erf(x)  = GammaP(1/2,x^2)
!   erfc(x) = GammaQ(1/2,x^2)
! 
ELEMENTAL FUNCTION ERF2(x1,x2) RESULT(z)
  IMPLICIT NONE
  REAL*8 :: z
  REAL*8, INTENT(IN) :: x1,x2
  REAL*8 :: xc,delx
  REAL*8, PARAMETER :: SQRTPI = 1.7724538509055160d0
  ! 1-eps is approximate level of cancelation at which to handle
  ! the difference more carefully
  REAL*8, PARAMETER :: EPS = 0.03d0
  
  ! Opposite sign: no canceling to worry about here
  IF (x1*x2 .LE. 0d0) THEN
    z = erf(x2) - erf(x1)
    RETURN
  END IF
  
  xc   = 0.5d0 * (x1+x2)
  delx = 0.5d0 * (x2-x1)
  
  ! Smaller arguments:
  !   |xc| < 1    --> |x1|,|x2| < 2
  ! Canceling is significant if |delx| < eps*|xc|
  IF ((ABS(xc) .LE. 1d0) .AND. (ABS(delx) .GT. EPS*ABS(xc))) THEN
    ! erf(x2) - erf(x1)
    z = erf(x2) - erf(x1)
    RETURN
    
  ! At least one argument is "large": 
  !   |xc| > 1    --> |x1| > 1 or |x2| > 1
  ! Canceling is significant if |4*xc*delx| < eps
  ELSE IF ((ABS(xc) .GT. 1d0) .AND. (ABS(4*xc*delx) .GT. EPS)) THEN
    IF (xc .GT. 0d0) THEN
      ! Difference of complementary error function gives better
      ! precision here:
      ! erf(x2) - erf(x1) = erfc(x1) - erfc(x2)
      z = -(erfc(x2) - erfc(x1))
    ELSE
      ! Difference of complementary error function gives better
      ! precision here (use symmetry of erf(x)):
      ! erf(x2) - erf(x1) = erf(-x1) - erf(-x2) = erfc(-x2) - erfc(-x1)
      z = erfc(-x2) - erfc(-x1)
    END IF
    RETURN
  END IF
  
  ! If we reached this point, there is significant canceling.
  ! For these cases, the integrand in the error function does not
  ! change much over x1 < x < x2.  Taylor expand the integrand
  ! about xc = (x1+x2)/2 and integrate.  The following keeps terms
  ! up through tenth order.
  z = 4 * delx * EXP(-xc**2) / SQRTPI                                   &
        * (1 + (2*xc**2 - 1)*delx**2 / 3                                &
             + (4*xc**4 - 12*xc**2 + 3)*delx**4 / 30                    &
             + (8*xc**6 - 60*xc**4 + 90*xc**2 - 15)*delx**6 / 630       &
             + (16*xc**8 - 224*xc**6 + 840*xc**4 - 840*xc**2            &
                         + 105)*delx**8 / 22680                         &
          )
  
END FUNCTION


!-----------------------------------------------------------------------
! The logarithm of the error function.  Requires x > 0.
! Accounts for large x case to avoid loss of precision.
! 
ELEMENTAL FUNCTION LOG_ERF(x) RESULT(z)
  IMPLICIT NONE
  REAL*8 :: z
  REAL*8, INTENT(IN) :: x
  
  ! Invalid input
  IF (x .LE. 0d0) THEN
    z = -HUGE(z)
    RETURN
  END IF
  
  ! Negligible loss of precision for smaller x
  IF (x .LE. 1d0) THEN
    z = LOG(ERF(x))
  ! Work with smaller complementary error function
  ELSE
    z = LOGp1(-ERFC(x))
  END IF
  
END FUNCTION


!-----------------------------------------------------------------------
! The logarithm of the complementary error function.
! Accounts for large x case to avoid loss of precision.
! 
ELEMENTAL FUNCTION LOG_ERFC(x) RESULT(z)
  IMPLICIT NONE
  REAL*8 :: z
  REAL*8, INTENT(IN) :: x
  REAL*8 :: y,w
  REAL*8, PARAMETER :: SQRTPI = 1.7724538509055160d0  ! Sqrt(Pi)
  
  ! Work with smaller error function for smaller |x|
  IF (ABS(x) .LE. 0.1d0) THEN
    z = LOGp1(-ERF(x))
  ! Negligible loss of precision for x not too large
  ELSE IF (x .LE. 25d0) THEN
    z = LOG(ERFC(x))
  ! Use asymptotic expansion:
  !   w = sqrt(pi) x e^{x^2} erfc(x) - 1
  !       ->  \sum_{k=1} (-1)^k (2k-1)!!/(2x^2)^k
  ! For x > 25, double precision in eight terms (k <= 8)
  ELSE
    y = 1 / (2*x**2)
    w = -y*(1-3*y*(1-5*y*(1-7*y*(1-9*y*(1-11*y*(1-13*y*(1-15*y)))))))
    z = -x**2 - LOGp1(w) - LOG(SQRTPI*x)
  END IF
  
END FUNCTION


!-----------------------------------------------------------------------
! The inverse error function, finding x given y such that y = ERF(x).
! 
! This implementation uses Halley's method.  It seems to be reasonable
! in terms of speed and precision and has had a cursory check for
! robustness, but there may be more optimal algorithms than this.
! Precision seems to be within ~10*EPSILON (usually even closer).
! 
ELEMENTAL FUNCTION ERFINV(y) RESULT(x)
  IMPLICIT NONE
  REAL*8 :: x
  REAL*8, INTENT(IN) :: y
  INTEGER :: K
  REAL*8 :: y0,y2,w,fk,delta
  REAL*8, PARAMETER :: QUARTERPI  = 0.78539816339744831d0  ! Pi/4
  REAL*8, PARAMETER :: INVSQRTPI  = 0.56418958354775629d0  ! 1/Sqrt(Pi)
  
  y0 = ABS(y)
  y2 = y*y
  
  ! Initial guess.  Taken from B.A. Popov, ACM SIGSAM 34, 25 (2000).
  IF (y0 .EQ. 0d0) THEN
    x = 0d0
    RETURN
  ELSE IF (y0 .LE. 0.97314979d0) THEN
    x = y*(-0.95493118d0+y2*(0.53160534d0+y2*0.23343441d0)) / (-1.0977154d0+y2)
  ELSE IF (y0 .LE. 0.99767065d0) THEN
    x = y*(1.6200516d0+y2*(-4.9295187d0+y2*3.2890636d0)) / (-1.0083317d0+y2)
  ELSE IF (y0 .LE. 0.99978842d0) THEN
    x = y*(29.849915d0+y2*(-61.896833d0+y2*32.044810d0)) / (-1.0007339d0+y2)
  ELSE IF (y0 .LT. 1d0) THEN
    w = -LOG(1-y2)
    x = SIGN(SQRT(w - 0.5d0*LOG(QUARTERPI*w)),y)
  ELSE
    x = SIGN(HUGE(x),y)
    RETURN
  END IF
  
  ! Halley's method.  Checked by hand that double precision is achieved
  ! by three iterations for all tested cases (epsilon, 1-epsilon, +/-,
  ! etc).  No safety checks as the algorithm appears to be convergent
  ! in all cases (though not proved).
  DO K=1,3
    fk    = ERF(x) - y
    ! Newton's
    !delta = - fk / (2*INVSQRTPI*EXP(-x**2))
    ! Halley's
    delta = fk / (x*fk - 2*INVSQRTPI*EXP(-x**2))
    x     = x + delta
  END DO
  
END FUNCTION


!-----------------------------------------------------------------------
! The inverse complementary error function, finding x given y such that
! y = ERFC(x).
! 
! This is a quick implementation that has not been optimized for speed
! and has not been tested for precision (though it should be close to
! double precision in most cases).
! 
ELEMENTAL FUNCTION ERFCINV(y) RESULT(x)
  IMPLICIT NONE
  REAL*8 :: x
  REAL*8, INTENT(IN) :: y
  INTEGER :: K
  REAL*8 :: fk,delta
  REAL*8, PARAMETER :: INVPI      = 0.31830988618379067d0  ! 1/Pi
  REAL*8, PARAMETER :: INVSQRTPI  = 0.56418958354775629d0  ! 1/Sqrt(Pi)

  ! Use ERFINV via the relation erfc^-1(1-z) = erf^-1(z).
  ! Below algorith (Halley's method, 5 iterations) works fine for
  ! y < 0.5 and possibly to higher y values, but ERFINV is faster
  ! (by x2-5), so we use that where loss of precision is minimal.
  IF (y .GT. 0.01d0) THEN
    x = ERFINV(1d0-y)
    RETURN
  ELSE IF (y .LE. 0d0) THEN
    x = HUGE(x)
    RETURN
  END IF

  ! Only small, positive y at this point.

  ! Initial guess: invert first term of asymptotic expansion
  !   y = erfc(x) = e^{-x^2} / (x\sqrt{\pi}) [1 + O(1/x^2)]
  x = SQRT(0.5d0*LAMBERTW(2*INVPI/y**2))

  ! Halley's method.  Checked by hand that double precision is achieved
  ! by five iterations for tested cases (epsilon, ~ 0.5).  No safety
  ! checks as the algorithm appears to be convergent in all cases
  ! (though not proved).
  DO K=1,5
    fk    = ERFC(x) - y
    ! Newton's
    !delta = fk / (2*INVSQRTPI*EXP(-x**2))
    ! Halley's
    delta = fk / (x*fk + 2*INVSQRTPI*EXP(-x**2))
    x     = x + delta
  END DO

END FUNCTION


!-----------------------------------------------------------------------
! The value of exp(x2)-exp(x1).
! Accounts for the case when x2-x1 is small and loss of precision
! due to canceling of the two terms can occur.
ELEMENTAL FUNCTION EXP2(x1,x2) RESULT(z)
  IMPLICIT NONE
  REAL*8 :: z
  REAL*8, INTENT(IN) :: x1,x2
  REAL*8 :: x,yk
  INTEGER :: K
  x = x2-x1
  ! Good to ~ 100*EPSILON precision (matched to sum terms below)
  IF (ABS(x) .GT. 0.007d0) THEN
    z = EXP(x2) - EXP(x1)
    RETURN
  END IF
  IF (x1 .EQ. x2) THEN
    z = 0d0
    RETURN
  END IF
  ! Can write: e^x2 - e^x1 = e^x1 [e^(x2-x1) - 1]
  ! Here, we find e^(x2-x1) - 1
  K = 1
  yk = x
  z = yk
  ! Nearly full precision, but does check limit optimization?
  !DO WHILE (ABS(yk) .GT. EPSILON(1d0)*ABS(z))
  !  K = K+1
  !  yk = yk*x/K
  !  z = z + yk
  !END DO
  ! Optimization: might be faster to calculate a fixed number
  !               of terms rather than check after each term if
  !               further terms are required.
  ! Possibilities:
  !     Good to ~   20*EPSILON in 7 terms for |x| < 0.04.
  !     Good to ~   50*EPSILON in 6 terms for |x| < 0.02.
  !     Good to ~  100*EPSILON in 5 terms for |x| < 0.007.
  !     Good to ~  500*EPSILON in 4 terms for |x| < 0.002.
  !     Good to ~ 4000*EPSILON in 3 terms for |x| < 0.00025.
  DO K = 2,5
    yk = yk*x/K
    z = z + yk
  END DO
  z = EXP(x1)*z
END FUNCTION


!-----------------------------------------------------------------------
! The value of exp(x)-1.
! Accounts for the case when x is small and loss of precision due
! to canceling of the two terms can occur.
ELEMENTAL FUNCTION EXPm1(x) RESULT(z)
  IMPLICIT NONE
  REAL*8 :: z
  REAL*8, INTENT(IN) :: x
  z = EXP2(0d0,x)
END FUNCTION


! ----------------------------------------------------------------------
! Function to calculate the quantity:
!   ln(1+x)
! Accounts for the case where x is small.
! 
! Precision is ~ 100*EPSILON.
! 
ELEMENTAL FUNCTION LOGp1(x) RESULT(z)
  IMPLICIT NONE
  REAL*8 :: z
  REAL*8, INTENT(IN) :: x
  REAL*8 :: xabs,xk
  INTEGER :: k,klast
  
  xabs = ABS(x)
  
  ! If x is not too small, we can just evaluate explicitly
  IF (xabs .GT. 0.01d0) THEN
    z = LOG(1+x)
    RETURN
  END IF
  
  ! We will use an expansion about zero, keeping terms in the
  ! expansion up to x^klast/klast.
  ! klast is chosen below to give almost double precision.
  ! Precision can be as poor as ~ 100*EPSILON, but is better
  ! for smaller x.  The expansion is:
  !     \sum_{k=1}^{\infty} (-1)^{k+1} x^k / k
  ! The precision is approximately:
  !     x^klast / (klast+1)
  ! where klast is the last term in the sum.
  
  IF (xabs .LT. 1d-4) THEN
    klast = 4
  ELSE
    klast = 8
  END IF
  ! Go to klast=12 for |x| < 0.05
  ! Go to klast=16 for |x| < 0.1
  
  ! Use expansion about zero
  xk = x
  z  = xk
  DO k = 2,klast
    xk      = -x*xk
    z = z + xk/k
  END DO
  
END FUNCTION LOGp1


! ----------------------------------------------------------------------
! Function to calculate the quantity ln(a+b) given ln(a) and ln(b).
! Precision is ~ EPSILON.
! 
ELEMENTAL FUNCTION LOG_SUM(lna,lnb) RESULT(z)
  IMPLICIT NONE
  REAL*8 :: z
  REAL*8, INTENT(IN) :: lna,lnb
  REAL*8 :: lnx,r
  
  ! Write a+b as x(1+r) with x = max(a,b) and r = min(b/a,a/b).
  IF (lna .GE. lnb) THEN
    lnx = lna
    r   = EXP(lnb-lna)
  ELSE
    lnx = lnb
    r   = EXP(lna-lnb)
  END IF
  
  ! Below, we use an expansion of ln(1+r) if r is small.
  ! The sum is terminated at approximately double precision.
  IF (r .EQ. 0d0) THEN
    z = lnx
  ELSE IF (r .GT. 0.01d0) THEN
    z = lnx + LOG(1d0 + r)
  ELSE IF (r .GT. 0.001d0) THEN
    z = lnx + r*(1d0                                                    &
                 -r*((1d0/2d0)                                          &
                     -r*((1d0/3d0)                                      &
                         -r*((1d0/4d0)                                  &
                             -r*((1d0/5d0)                              &
                                 -r*((1d0/6d0)                          &
                                     -r*((1d0/7d0)                      &
                                         -r*((1d0/8d0)                  &
                                             -r*(1d0/9d0)               &
                ))))))))
  ELSE IF (r .GT. 0.00001d0) THEN
    z = lnx + r*(1d0                                                    &
                 -r*((1d0/2d0)                                          &
                     -r*((1d0/3d0)                                      &
                         -r*((1d0/4d0)                                  &
                             -r*((1d0/5d0)                              &
                                 -r*(1d0/6d0)                           &
                )))))
  ELSE
    z = lnx + r*(1d0                                                    &
                 -r*((1d0/2d0)                                          &
                     -r*((1d0/3d0)                                      &
                         -r*(1d0/4d0)                                   &
                )))
  END IF
  
END FUNCTION


! ----------------------------------------------------------------------
! Function to calculate the quantity ln(a-b) given ln(a) and ln(b).
! Requires a > b.  Precision is ~ EPSILON.
! 
ELEMENTAL FUNCTION LOG_DIFF(lna,lnb) RESULT(z)
  IMPLICIT NONE
  REAL*8 :: z
  REAL*8, INTENT(IN) :: lna,lnb
  REAL*8 :: lnx,r
  
  ! Bad case
  IF (lnb .GE. lna) THEN
    z = -HUGE(z)
    RETURN
  END IF
  
  ! Write a-b as x(1-r) with x = a and r = b/a.
  lnx = lna
  r   = EXP(lnb-lna)
  
  ! Below, we use an expansion of ln(1-r) if r is small.
  ! The sum is terminated at approximately double precision.
  IF (r .EQ. 0d0) THEN
    z = lnx
  ELSE IF (r .GT. 0.01d0) THEN
    z = lnx + LOG(1d0 - r)
  ELSE IF (r .GT. 0.001d0) THEN
    z = lnx - r*(1d0                                                    &
                 +r*((1d0/2d0)                                          &
                     +r*((1d0/3d0)                                      &
                         +r*((1d0/4d0)                                  &
                             +r*((1d0/5d0)                              &
                                 +r*((1d0/6d0)                          &
                                     +r*((1d0/7d0)                      &
                                         +r*((1d0/8d0)                  &
                                             +r*(1d0/9d0)               &
                ))))))))
  ELSE IF (r .GT. 0.00001d0) THEN
    z = lnx - r*(1d0                                                    &
                 +r*((1d0/2d0)                                          &
                     +r*((1d0/3d0)                                      &
                         +r*((1d0/4d0)                                  &
                             +r*((1d0/5d0)                              &
                                 +r*(1d0/6d0)                           &
                )))))
  ELSE
    z = lnx - r*(1d0                                                    &
                 +r*((1d0/2d0)                                          &
                     +r*((1d0/3d0)                                      &
                         +r*(1d0/4d0)                                   &
                )))
  END IF
  
END FUNCTION


!-----------------------------------------------------------------------
! Gamma function of integer argument [double precision]
! NOTE: Intrinsic routine of real argument available in Fortran 2008;
!       the speed of the intrinsic routine may be comparable to this
!       tabulated case.
! 
ELEMENTAL FUNCTION GAMMAI(n) RESULT(z)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: n
  REAL*8 :: z
  REAL*8 :: x,xinv
  INTEGER, PARAMETER :: NMAX  = 171    ! n > NMAX exceeds HUGE(1d0)
  ! Tabulated values
  ! NOTE: number of continuation lines required to place all tabulated
  !       values into one table exceeds F95 spec (>39).
  !       Both Intel and GNU Fortran compilers allow more continuation
  !       lines, but we use multiple tables for portability.
  INTEGER, PARAMETER :: NMAX1 = 100
  INTEGER, PARAMETER :: NMAX2 = 171
  ! Table for n=1-100
  REAL*8, PARAMETER, DIMENSION(1:NMAX1) :: GVALS1 = &
    (/  1.0d0,                  1.0d0,                  2.0d0,                  6.0d0,                  & ! 1-4
        24.0d0,                 120.0d0,                720.0d0,                5040.0d0,               & ! 5-8
        40320.0d0,              362880.0d0,             3.6288000000000000d6,   3.9916800000000000d7,   & ! 9-12
        4.7900160000000000d8,   6.2270208000000000d9,   8.7178291200000000d10,  1.3076743680000000d12,  & ! 13-16
        2.0922789888000000d13,  3.5568742809600000d14,  6.4023737057280000d15,  1.2164510040883200d17,  & ! 17-20
        2.4329020081766400d18,  5.1090942171709440d19,  1.1240007277776077d21,  2.5852016738884977d22,  & ! 21-24
        6.2044840173323944d23,  1.5511210043330986d25,  4.0329146112660564d26,  1.0888869450418352d28,  & ! 25-28
        3.0488834461171386d29,  8.8417619937397020d30,  2.6525285981219106d32,  8.2228386541779228d33,  & ! 29-32
        2.6313083693369353d35,  8.6833176188118865d36,  2.9523279903960414d38,  1.0333147966386145d40,  & ! 33-36
        3.7199332678990122d41,  1.3763753091226345d43,  5.2302261746660111d44,  2.0397882081197443d46,  & ! 37-40
        8.1591528324789773d47,  3.3452526613163807d49,  1.4050061177528799d51,  6.0415263063373836d52,  & ! 41-44
        2.6582715747884488d54,  1.1962222086548019d56,  5.5026221598120889d57,  2.5862324151116818d59,  & ! 45-48
        1.2413915592536073d61,  6.0828186403426756d62,  3.0414093201713378d64,  1.5511187532873823d66,  & ! 49-52
        8.0658175170943879d67,  4.2748832840600256d69,  2.3084369733924138d71,  1.2696403353658276d73,  & ! 53-56
        7.1099858780486345d74,  4.0526919504877217d76,  2.3505613312828786d78,  1.3868311854568984d80,  & ! 57-60
        8.3209871127413901d81,  5.0758021387722480d83,  3.1469973260387938d85,  1.9826083154044401d87,  & ! 61-64
        1.2688693218588416d89,  8.2476505920824707d90,  5.4434493907744306d92,  3.6471110918188685d94,  & ! 65-68
        2.4800355424368306d96,  1.7112245242814131d98,  1.1978571669969892d100, 8.5047858856786232d101, & ! 69-72
        6.1234458376886087d103, 4.4701154615126843d105, 3.3078854415193864d107, 2.4809140811395398d109, & ! 73-76
        1.8854947016660503d111, 1.4518309202828587d113, 1.1324281178206298d115, 8.9461821307829753d116, & ! 77-80
        7.1569457046263802d118, 5.7971260207473680d120, 4.7536433370128417d122, 3.9455239697206587d124, & ! 81-84
        3.3142401345653533d126, 2.8171041143805503d128, 2.4227095383672732d130, 2.1077572983795277d132, & ! 85-88
        1.8548264225739844d134, 1.6507955160908461d136, 1.4857159644817615d138, 1.3520015276784030d140, & ! 89-92
        1.2438414054641307d142, 1.1567725070816416d144, 1.0873661566567431d146, 1.0329978488239059d148, & ! 93-96
        9.9167793487094969d149, 9.6192759682482120d151, 9.4268904488832477d153, 9.3326215443944153d155 /) ! 97-100
  ! Table for n=101-171
  REAL*8, PARAMETER, DIMENSION(NMAX1+1:NMAX2) :: GVALS2 = &
    (/  9.3326215443944153d157, 9.4259477598383594d159, 9.6144667150351266d161, 9.9029007164861804d163, & ! 101-104
        1.0299016745145628d166, 1.0813967582402909d168, 1.1462805637347084d170, 1.2265202031961379d172, & ! 105-108
        1.3246418194518290d174, 1.4438595832024936d176, 1.5882455415227429d178, 1.7629525510902447d180, & ! 109-112
        1.9745068572210740d182, 2.2311927486598136d184, 2.5435597334721876d186, 2.9250936934930157d188, & ! 113-116
        3.3931086844518982d190, 3.9699371608087209d192, 4.6845258497542907d194, 5.5745857612076059d196, & ! 117-120
        6.6895029134491271d198, 8.0942985252734437d200, 9.8750442008336014d202, 1.2146304367025330d205, & ! 121-124
        1.5061417415111409d207, 1.8826771768889261d209, 2.3721732428800469d211, 3.0126600184576595d213, & ! 125-128
        3.8562048236258042d215, 4.9745042224772874d217, 6.4668554892204737d219, 8.4715806908788205d221, & ! 129-132
        1.1182486511960043d224, 1.4872707060906857d226, 1.9929427461615189d228, 2.6904727073180505d230, & ! 133-136
        3.6590428819525487d232, 5.0128887482749917d234, 6.9177864726194885d236, 9.6157231969410890d238, & ! 137-140
        1.3462012475717525d241, 1.8981437590761710d243, 2.6953641378881628d245, 3.8543707171800728d247, & ! 141-144
        5.5502938327393048d249, 8.0479260574719919d251, 1.1749972043909108d254, 1.7272458904546389d256, & ! 145-148
        2.5563239178728656d258, 3.8089226376305697d260, 5.7133839564458546d262, 8.6272097742332404d264, & ! 149-152
        1.3113358856834525d267, 2.0063439050956824d269, 3.0897696138473509d271, 4.7891429014633939d273, & ! 153-156
        7.4710629262828944d275, 1.1729568794264144d278, 1.8532718694937348d280, 2.9467022724950383d282, & ! 157-160
        4.7147236359920613d284, 7.5907050539472187d286, 1.2296942187394494d289, 2.0044015765453026d291, & ! 161-164
        3.2872185855342962d293, 5.4239106661315888d295, 9.0036917057784374d297, 1.5036165148649990d300, & ! 165-168
        2.5260757449731984d302, 4.2690680090047053d304, 7.2574156153079990d306  /)                        ! 169-171
  ! Coefficients for asymptotic expansion
  ! Not all terms are necessary or used here
  INTEGER, PARAMETER :: NC = 20
  REAL*8, PARAMETER, DIMENSION(0:NC) :: C = &
    (/  1.0d0,                                                  &
        0.083333333333333333d0,     0.0034722222222222222d0,    &
       -0.0026813271604938272d0,   -0.00022947209362139918d0,   &
        0.00078403922172006663d0,   0.000069728137583658578d0,  &
       -0.00059216643735369388d0,  -0.000051717909082605922d0,  &
        0.00083949872067208728d0,   0.000072048954160200106d0,  &
       -0.0019144384985654775d0,   -0.00016251626278391582d0,   &
        0.0064033628338080698d0,    0.00054016476789260452d0,   &
       -0.029527880945699121d0,    -0.0024817436002649977d0,    &
        0.17954011706123486d0,      0.015056113040026424d0,     &
       -1.3918010932653375d0,      -0.1165462765994632d0       /)
  
  IF (n .LE. 0) THEN
    z = -HUGE(z)
    !STOP 'ERROR: gamma cannot be called with non-positive argument'
  ELSE IF (n .GT. NMAX) THEN
    z = HUGE(z)
  ELSE IF (n .LE. NMAX1) THEN
    z = GVALS1(n)
  ELSE IF (n .LE. NMAX2) THEN
    z = GVALS2(n)
  ELSE
    ! This case should not occur
    z = HUGE(z)
  END IF
  ! Algorithm here unnecessary, but shown for reference
  !! Asymptotic expansion (Stirling's)
  !! Error will be ~ C[5]/x^5
  !! Good to ~ 1e-15 for n > 171
  !! Result would be HUGE(1d0) for n > 171
  !x = n
  !xinv = 1d0/x
  !z = x**(x-0.5d0) * EXP(-x) * SQRT2PI                               &
  !    * (C(0) +                                                      &
  !       xinv*(C(1) +                                                &
  !             xinv*(C(2) +                                          &
  !                   xinv*(C(3) +                                    &
  !                         xinv*(C(4) +                              &
  !                               xinv*C(5) ) ) ) ) )
  !! Alternatively, could use LogGamma function for large n
  !!x = n
  !!z = EXP(LOG_GAMMA(x))
END FUNCTION


!-----------------------------------------------------------------------
! Lower gamma function of real arguments [double precision]
!   \gamma(s,x) = \int_0^x dt t^{s-1} e^{-t}
! 
ELEMENTAL FUNCTION LOWER_GAMMA(s,x) RESULT(z)
  IMPLICIT NONE
  REAL*8, INTENT(IN) :: s,x
  REAL*8 :: z
  REAL*8 :: P,Q
  CALL GAMMA_PQ(s,x,P,Q)
  z = GAMMA(s) * P
END FUNCTION


!-----------------------------------------------------------------------
! Upper gamma function of real arguments [double precision]
!   \Gamma(s,x) = \int_x^\infty dt t^{s-1} e^{-t}
! 
ELEMENTAL FUNCTION UPPER_GAMMA(s,x) RESULT(z)
  IMPLICIT NONE
  REAL*8, INTENT(IN) :: s,x
  REAL*8 :: z
  REAL*8 :: P,Q
  CALL GAMMA_PQ(s,x,P,Q)
  z = GAMMA(s) * Q
END FUNCTION


!-----------------------------------------------------------------------
! Regularized gamma function P of real arguments [double precision]
!   P(s,x) = \gamma(s,x) / \Gamma(s)
! 
ELEMENTAL FUNCTION GAMMA_P(s,x) RESULT(z)
  IMPLICIT NONE
  REAL*8, INTENT(IN) :: s,x
  REAL*8 :: z
  REAL*8 :: P,Q
  CALL GAMMA_PQ(s,x,P,Q)
  z = P
END FUNCTION


!-----------------------------------------------------------------------
! Regularized gamma function Q of real arguments [double precision]
!   Q(s,x) = \Gamma(s,x) / \Gamma(s)
! 
ELEMENTAL FUNCTION GAMMA_Q(s,x) RESULT(z)
  IMPLICIT NONE
  REAL*8, INTENT(IN) :: s,x
  REAL*8 :: z
  REAL*8 :: P,Q
  CALL GAMMA_PQ(s,x,P,Q)
  z = Q
END FUNCTION


!-----------------------------------------------------------------------
! Regularized incomplete gamma functions P, Q of real arguments
! [double precision].  These are defined as:
!   P(s,x) = \gamma(s,x) / \Gamma(s)
!   Q(s,x) = \Gamma(s,x) / \Gamma(s)
! where \gamma(s,x) & \Gamma(s,x) are the lower & upper incomplete
! gamma functions:
!   \gamma(s,x) = \int_0^x dt t^{s-1} e^{-t}
!   \Gamma(s,x) = \int_x^\infty dt t^{s-1} e^{-t}
! The following relations hold:
!   P(s,x) + Q(s,x) = 1
!   \gamma(s,x) + \Gamma(s,x) = \Gamma(s)
! The below routine calculates the approximately smaller of P or Q
! and uses P+Q=1 to determine the other.
! 
! NOTE: Valid only for s > 0 and x >= 0.
! 
ELEMENTAL SUBROUTINE GAMMA_PQ(s,x,P,Q)
  IMPLICIT NONE
  REAL*8, INTENT(IN) :: s,x
  REAL*8, INTENT(OUT) :: P,Q
  
  ! Special cases
  IF (x .EQ. 0d0) THEN
    P = 0d0
    Q = 1d0
    RETURN
  END IF
  
  ! Bad cases
  IF ((x .LT. 0d0) .OR. (s .LT. 0d0)) THEN
    P = -HUGE(1d0)
    Q = HUGE(1d0)
    RETURN
  END IF
  
  ! Calculate P,Q using uniform asymptotic expansion
  IF ((s .GE. 10d0) .AND. (x-10d0 .GT. 0.302d0*(s-10d0)) &
      .AND. (x-10d0 .LT. 2.357d0*(s-10d0))) THEN
    CALL GAMMA_PQ_UA(s,x,P,Q)
  ! Calculate P using Taylor series
  ELSE IF (x .LE. s) THEN
    P = GAMMA_P_TS(s,x)
    Q = 1-P
  ! Calculate Q using continued fraction
  ELSE
    Q = GAMMA_Q_CF(s,x)
    P = 1-Q
  END IF
  
  CONTAINS
  
  !---------------------------------------------
  ! Calculates gamma function quantity P(s,x) using a Taylor
  ! series expansion.
  !   P(s,x) = \gamma(s,x) / \Gamma(s)
  !   \gamma(s,x) = \int_0^x dt t^{s-1} e^{-t}
  ! NOTE: Requires 0 < s, 0 <= x.  Should only be used for x
  !       not (much) larger than s or convergence may be slow.
  !       Also best for s not large (timing scales as maybe
  !       sqrt(s)).
  PURE FUNCTION GAMMA_P_TS(s,x) RESULT(P)
    IMPLICIT NONE
    REAL*8, INTENT(IN) :: s,x
    REAL*8 :: P
    INTEGER :: K
    REAL*8 :: zk,sum
    ! P(s,x) = x^s e^{-x} \Sum_{k=0}^\infty x^k / \Gamma(s+k+1)
    K   = 0
    zk  = 1d0
    sum = zk
    DO WHILE (zk .GT. EPSILON(zk))
      K   = K+1
      zk  = zk * x / (s+K)
      sum = sum + zk
      IF (K .GE. 10000) EXIT
    END DO
    P = sum * EXP(s*LOG(x) - x - LOG_GAMMA(s+1))
  END FUNCTION
  
  !---------------------------------------------
  ! Calculates gamma function quantity Q(s,x) using a
  ! continued fraction.
  !   Q(s,x) = \Gamma(s,x) / \Gamma(s)
  !   \Gamma(s,x) = \int_x^{\infty} dt t^{s-1} e^{-t}
  ! NOTE: Requires 0 < x,s.  Should only be used for x not
  !       (much) smaller than s or convergence may be slow
  !       (or simply fail).  Also best for s not large
  !       (timing scales as maybe sqrt(s)).
  PURE FUNCTION GAMMA_Q_CF(s,x) RESULT(Q)
    IMPLICIT NONE
    REAL*8, INTENT(IN) :: s,x
    REAL*8 :: Q
    INTEGER :: K
    REAL*8 :: xs1,fk,Ck,Dk,Deltak
    ! Continued fraction (see Numerical Recipes):
    !   Q(s,x) = x^s e^{-x} / \Gamma(s)
    !            / (x-s+1 + K_k(-k(k-s),-s+2k+x+1)_{1}^{\infty})
    ! where
    !   K_k(a_k,b_k)_1 = a_1/b1+ a2/b2+ a3/b3+ ...
    ! is the Gauss notation for a continued fraction.
    xs1 = x-s+1
    K  = 0
    Ck = xs1
    Dk = 0
    fk = xs1
    Deltak = HUGE(x)
    DO WHILE (ABS(Deltak-1) .GE. EPSILON(x))
      K  = K+1
      Ck = (xs1+2*K) - K*(K-s)/Ck
      IF (Ck .EQ. 0d0) Ck = 1d-30
      Dk = (xs1+2*K) - K*(K-s)*Dk
      IF (Dk .EQ. 0d0) Dk = 1d-30
      Dk = 1/Dk
      Deltak = Ck*Dk
      fk = Deltak*fk
      IF (K .GE. 10000) EXIT
    END DO
    Q = EXP(s*LOG(x) - x - LOG_GAMMA(s)) / fk
  END FUNCTION
  
  !---------------------------------------------
  ! Calculates gamma function quantities P(s,x) and Q(s,x)
  ! using a uniform asymptotic expansion (see 1306.1754).
  ! Nearly constant evaluation time for any s & x, unlike
  ! the Taylor series or continued fraction algorithms.
  ! NOTE: Intended for s not small (larger than ~10) and
  !       0.30 < x/s < 2.35; accuracy starts to drop
  !       outside this region.
  PURE SUBROUTINE GAMMA_PQ_UA(s,x,P,Q)
    IMPLICIT NONE
    REAL*8, INTENT(IN) :: s,x
    REAL*8, INTENT(OUT) :: P,Q
    INTEGER :: K
    REAL*8 :: eta,lambda,u,Sa,Ra,sum
    REAL*8 :: beta(0:32)
    INTEGER, PARAMETER :: KMAX = 25
    REAL*8, PARAMETER :: D(0:31) = &
      (/  1.0000000000000000d+00,-3.3333333333333333d-01, 8.3333333333333333d-02, &
         -1.4814814814814815d-02, 1.1574074074074074d-03, 3.5273368606701940d-04, &
         -1.7875514403292181d-04, 3.9192631785224378d-05,-2.1854485106799922d-06, &
         -1.8540622107151600d-06, 8.2967113409530860d-07,-1.7665952736826079d-07, &
          6.7078535434014986d-09, 1.0261809784240308d-08,-4.3820360184533532d-09, &
          9.1476995822367902d-10,-2.5514193994946250d-11,-5.8307721325504251d-11, &
          2.4361948020667416d-11,-5.0276692801141756d-12, 1.1004392031956135d-13, &
          3.3717632624009854d-13,-1.3923887224181621d-13, 2.8534893807047443d-14, &
         -5.1391118342425726d-16,-1.9752288294349443d-15, 8.0995211567045613d-16, &
         -1.6522531216398162d-16, 2.5305430097478884d-18, 1.1686939738559577d-17, &
         -4.7700370498204848d-18, 9.6991260590562371d-19 /)
    REAL*8, PARAMETER :: TWOPI = 6.2831853071795865d0   ! 2*Pi
    
    ! Define ratio of x & s
    lambda = x/s
    
    ! eta is defined as 1/2 eta^2 = lambda - 1 - ln(lambda) with
    ! the same sign as lambda-1.  Use an expansion for lambda ~ 1
    ! to avoid cancellation issues.  Rapid convergent of below
    ! series in eta requires |eta| < 1, which corresponds to
    ! 0.3017 < lambda < 2.3577.
    u = lambda - 1
    IF (ABS(u) .GT. 0.01d0) THEN
      eta = SQRT(2*(lambda-1-LOG(lambda)))
    ELSE
      eta = u*(1+u*(-1d0/3d0+u*(7d0/36d0+u*(-73d0/540d0+u*(1331d0/12960d0)))))
    END IF
    IF (x .LT. s) eta = -eta
    
    ! Use erfc in eta*sqrt{s/2} as first approximation and calculate
    ! correction term:
    !   R_a(\eta) = e^{-1/2 s \eta^2} S_a(\eta) / sqrt{2\pi s}
    ! where Sa(eta) can be given by an expansion in eta:
    !   S_a(\eta) \approx a/(a+\beta_1) \Sum_{k=0}^N \beta_k \eta^k
    ! The beta terms are dependent upon s and can be calculated using
    ! a recursion relation.
    beta(KMAX+2) = 0d0
    beta(KMAX+1) = 0d0
    sum = 0d0
    DO K = KMAX,0,-1
      beta(K) = (K+2)/s * beta(K+2) + D(K+1)
      sum = beta(K) + eta*sum
    END DO
    Sa = s/(s+beta(1)) * sum
    Ra = EXP(-0.5d0*s*eta**2) * Sa / SQRT(TWOPI*s)
    
    ! Formulas for P & Q always valid, but just calculate the
    ! smaller one and use P+Q=1 for the other.
    IF (x .LE. s) THEN
      P = 0.5d0*ERFC(-eta*SQRT(0.5d0*s)) - Ra
      Q = 1-P
    ELSE
      Q = 0.5d0*ERFC(+eta*SQRT(0.5d0*s)) + Ra
      P = 1-Q
    END IF
    
  END SUBROUTINE
  
END SUBROUTINE


!-----------------------------------------------------------------------
! Logarithm of regularized gamma function P of real arguments
! [double precision]
!   P(s,x) = \gamma(s,x) / \Gamma(s)
! 
ELEMENTAL FUNCTION LOG_GAMMA_P(s,x) RESULT(z)
  IMPLICIT NONE
  REAL*8, INTENT(IN) :: s,x
  REAL*8 :: z
  REAL*8 :: lnP,lnQ
  CALL LOG_GAMMA_PQ(s,x,lnP,lnQ)
  z = lnP
END FUNCTION


!-----------------------------------------------------------------------
! Logarithm of regularized gamma function Q of real arguments
! [double precision]
!   Q(s,x) = \Gamma(s,x) / \Gamma(s)
! 
ELEMENTAL FUNCTION LOG_GAMMA_Q(s,x) RESULT(z)
  IMPLICIT NONE
  REAL*8, INTENT(IN) :: s,x
  REAL*8 :: z
  REAL*8 :: lnP,lnQ
  CALL LOG_GAMMA_PQ(s,x,lnP,lnQ)
  z = lnQ
END FUNCTION


!-----------------------------------------------------------------------
! Logarithm of regularized incomplete gamma functions P, Q of real
! arguments [double precision].  These are defined as:
!   P(s,x) = \gamma(s,x) / \Gamma(s)
!   Q(s,x) = \Gamma(s,x) / \Gamma(s)
! where \gamma(s,x) & \Gamma(s,x) are the lower & upper incomplete
! gamma functions:
!   \gamma(s,x) = \int_0^x dt t^{s-1} e^{-t}
!   \Gamma(s,x) = \int_x^\infty dt t^{s-1} e^{-t}
! The following relations hold:
!   P(s,x) + Q(s,x) = 1
!   \gamma(s,x) + \Gamma(s,x) = \Gamma(s)
! The below routine calculates the approximately smaller of P or Q
! and uses P+Q=1 to determine the other.
! 
! NOTE: Valid only for s > 0 and x >= 0.
! 
ELEMENTAL SUBROUTINE LOG_GAMMA_PQ(s,x,lnP,lnQ)
  IMPLICIT NONE
  REAL*8, INTENT(IN) :: s,x
  REAL*8, INTENT(OUT) :: lnP,lnQ
  
  ! Special cases
  IF (x .EQ. 0d0) THEN
    lnP = -HUGE(1d0)
    lnQ = 0d0
    RETURN
  END IF
  
  ! Bad cases
  IF ((x .LT. 0d0) .OR. (s .LE. 0d0)) THEN
    lnP = -HUGE(1d0)
    lnQ = HUGE(1d0)
    RETURN
  END IF
  
  ! Calculate P,Q using uniform asymptotic expansion
  IF ((s .GE. 10d0) .AND. (x-10d0 .GT. 0.302d0*(s-10d0)) &
      .AND. (x-10d0 .LT. 2.357d0*(s-10d0))) THEN
    CALL LOG_GAMMA_PQ_UA(s,x,lnP,lnQ)
  ! Calculate P using Taylor series
  ELSE IF (x .LE. s) THEN
    lnP = LOG_GAMMA_P_TS(s,x)
    lnQ = LOGp1(-EXP(lnP))
  ! Calculate Q using continued fraction
  ELSE
    lnQ = LOG_GAMMA_Q_CF(s,x)
    lnP = LOGp1(-EXP(lnQ))
  END IF
  
  CONTAINS
  
  !---------------------------------------------
  ! Calculates gamma function quantity P(s,x) using a Taylor
  ! series expansion.
  !   P(s,x) = \gamma(s,x) / \Gamma(s)
  !   \gamma(s,x) = \int_0^x dt t^{s-1} e^{-t}
  ! NOTE: Requires 0 < s, 0 <= x.  Should only be used for x
  !       not (much) larger than s or convergence may be slow.
  !       Also best for s not large (timing scales as maybe
  !       sqrt(s)).
  PURE FUNCTION LOG_GAMMA_P_TS(s,x) RESULT(lnP)
    IMPLICIT NONE
    REAL*8, INTENT(IN) :: s,x
    REAL*8 :: lnP
    INTEGER :: K
    REAL*8 :: zk,sum
    ! P(s,x) = x^s e^{-x} \Sum_{k=0}^\infty x^k / \Gamma(s+k+1)
    ! sum here excludes first term (which is 1)
    K   = 1
    zk  = x / (s+1)
    sum = zk
    DO WHILE (zk .GT. EPSILON(zk))
      K   = K+1
      zk  = zk * x / (s+K)
      sum = sum + zk
      IF (K .GE. 10000) EXIT
    END DO
    lnP = s*LOG(x) - x - LOG_GAMMA(s+1) + LOGp1(sum)
  END FUNCTION
  
  !---------------------------------------------
  ! Calculates gamma function quantity Q(s,x) using a
  ! continued fraction.
  !   Q(s,x) = \Gamma(s,x) / \Gamma(s)
  !   \Gamma(s,x) = \int_x^{\infty} dt t^{s-1} e^{-t}
  ! NOTE: Requires 0 < x,s.  Should only be used for x not
  !       (much) smaller than s or convergence may be slow
  !       (or simply fail).  Also best for s not large
  !       (timing scales as maybe sqrt(s)).
  PURE FUNCTION LOG_GAMMA_Q_CF(s,x) RESULT(lnQ)
    IMPLICIT NONE
    REAL*8, INTENT(IN) :: s,x
    REAL*8 :: lnQ
    INTEGER :: K
    REAL*8 :: xs1,fk,Ck,Dk,Deltak
    ! Continued fraction (see Numerical Recipes):
    !   Q(s,x) = x^s e^{-x} / \Gamma(s)
    !            / (x-s+1 + K_k(-k(k-s),-s+2k+x+1)_{1}^{\infty})
    ! where
    !   K_k(a_k,b_k)_1 = a_1/b1+ a2/b2+ a3/b3+ ...
    ! is the Gauss notation for a continued fraction.
    xs1 = x-s+1
    K  = 0
    Ck = xs1
    Dk = 0
    fk = xs1
    Deltak = HUGE(x)
    DO WHILE (ABS(Deltak-1) .GE. EPSILON(x))
      K  = K+1
      Ck = (xs1+2*K) - K*(K-s)/Ck
      IF (Ck .EQ. 0d0) Ck = 1d-30
      Dk = (xs1+2*K) - K*(K-s)*Dk
      IF (Dk .EQ. 0d0) Dk = 1d-30
      Dk = 1/Dk
      Deltak = Ck*Dk
      fk = Deltak*fk
      IF (K .GE. 10000) EXIT
    END DO
    lnQ = s*LOG(x) - x - LOG_GAMMA(s) - LOG(fk)
  END FUNCTION
  
  !---------------------------------------------
  ! Calculates gamma function quantities P(s,x) and Q(s,x)
  ! using a uniform asymptotic expansion (see 1306.1754).
  ! Nearly constant evaluation time for any s & x, unlike
  ! the Taylor series or continued fraction algorithms.
  ! NOTE: Intended for s not small (larger than ~10) and
  !       0.30 < x/s < 2.35; accuracy starts to drop
  !       outside this region.
  ! This routine is about ~ 2 slower than the non-
  ! logarithmic version.
  PURE SUBROUTINE LOG_GAMMA_PQ_UA(s,x,lnP,lnQ)
    IMPLICIT NONE
    REAL*8, INTENT(IN) :: s,x
    REAL*8, INTENT(OUT) :: lnP,lnQ
    INTEGER :: K,sgnRa
    REAL*8 :: eta,lambda,u,lnSa,lnRa,sum
    REAL*8 :: beta(0:32)
    INTEGER, PARAMETER :: KMAX = 25
    REAL*8, PARAMETER :: D(0:31) = &
      (/  1.0000000000000000d+00,-3.3333333333333333d-01, 8.3333333333333333d-02, &
         -1.4814814814814815d-02, 1.1574074074074074d-03, 3.5273368606701940d-04, &
         -1.7875514403292181d-04, 3.9192631785224378d-05,-2.1854485106799922d-06, &
         -1.8540622107151600d-06, 8.2967113409530860d-07,-1.7665952736826079d-07, &
          6.7078535434014986d-09, 1.0261809784240308d-08,-4.3820360184533532d-09, &
          9.1476995822367902d-10,-2.5514193994946250d-11,-5.8307721325504251d-11, &
          2.4361948020667416d-11,-5.0276692801141756d-12, 1.1004392031956135d-13, &
          3.3717632624009854d-13,-1.3923887224181621d-13, 2.8534893807047443d-14, &
         -5.1391118342425726d-16,-1.9752288294349443d-15, 8.0995211567045613d-16, &
         -1.6522531216398162d-16, 2.5305430097478884d-18, 1.1686939738559577d-17, &
         -4.7700370498204848d-18, 9.6991260590562371d-19 /)
    REAL*8, PARAMETER :: TWOPI = 6.2831853071795865d0   ! 2*Pi
    REAL*8, PARAMETER :: LN2   = 0.69314718055994531d0  ! ln(2)
    
    ! Define ratio of x & s
    lambda = x/s
    
    ! eta is defined as 1/2 eta^2 = lambda - 1 - ln(lambda) with
    ! the same sign as lambda-1.  Use an expansion for lambda ~ 1
    ! to avoid cancellation issues.  Rapid convergent of below
    ! series in eta requires |eta| < 1, which corresponds to
    ! 0.3017 < lambda < 2.3577.
    u = lambda - 1
    IF (ABS(u) .GT. 0.01d0) THEN
      eta = SQRT(2*(lambda-1-LOG(lambda)))
    ELSE
      eta = u*(1+u*(-1d0/3d0+u*(7d0/36d0+u*(-73d0/540d0+u*(1331d0/12960d0)))))
    END IF
    IF (x .LT. s) eta = -eta
    
    ! Use erfc in eta*sqrt{s/2} as first approximation and calculate
    ! correction term:
    !   R_a(\eta) = e^{-1/2 s \eta^2} S_a(\eta) / sqrt{2\pi s}
    ! where Sa(eta) can be given by an expansion in eta:
    !   S_a(\eta) \approx a/(a+\beta_1) \Sum_{k=0}^N \beta_k \eta^k
    ! The beta terms are dependent upon s and can be calculated using
    ! a recursion relation.
    beta(KMAX+2) = 0d0
    beta(KMAX+1) = 0d0
    sum = 0d0
    DO K = KMAX,0,-1
      beta(K) = (K+2)/s * beta(K+2) + D(K+1)
      sum = beta(K) + eta*sum
    END DO
    IF (sum .GE. 0d0) THEN
      sgnRa = +1
    ELSE
      sgnRa = -1
    END IF
    ! Assuming b_1 > -s, which seems to be the case.
    ! Have not checked if LOG(sum) loses accuracy here....
    lnSa = LOG(ABS(sum)) - LOGp1(beta(1)/s)
    lnRa = -0.5d0*s*eta**2 + lnSa - 0.5d0*LOG(TWOPI*s)
    
    ! Formulas for P & Q always valid, but just calculate the
    ! smaller one and use P+Q=1 for the other.
    IF (x .LE. s) THEN
      IF (sgnRa .GT. 0) THEN
        lnP = LOG_DIFF(LOG_ERFC(-eta*SQRT(0.5d0*s))-LN2,lnRa)
      ELSE
        lnP = LOG_SUM(LOG_ERFC(-eta*SQRT(0.5d0*s))-LN2,lnRa)
      END IF
      lnQ = LOGp1(-EXP(lnP))
    ELSE
      IF (sgnRa .GT. 0) THEN
        lnQ = LOG_SUM(LOG_ERFC(+eta*SQRT(0.5d0*s))-LN2,lnRa)
      ELSE
        lnQ = LOG_DIFF(LOG_ERFC(+eta*SQRT(0.5d0*s))-LN2,lnRa)
      END IF
      lnP = LOGp1(-EXP(lnQ))
    END IF
    
  END SUBROUTINE
  
END SUBROUTINE


!-----------------------------------------------------------------------
! Binomial coefficient (n k) = n!/k!(n-k)! [integer]
! Various techniques are used to provide a fast calculation in
! different cases, but no check for integer overflow is performed
! for the final result.
! 
ELEMENTAL FUNCTION BINOMIAL_COEFF(N,K) RESULT(Cnk)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: N,K
  INTEGER :: Cnk
  INTEGER :: K0,J
  ! Array of binomial coefficients C(n,k).
  ! For k <= n/2, coefficient is at index:
  !    J = floor[(n+1)/2]*floor[(n+2)/2] + k + 1
  ! For k > n/2, we will use C(n,k) = C(n,n-k).
  INTEGER, PARAMETER :: NMAX_TAB = 30    ! Coefficients for N in [0,NMAX_TAB]
  INTEGER, PARAMETER :: NC       = 256
  INTEGER, PARAMETER :: BC_ARRAY(1:NC) =                                &
    (/  1,1,1,2,1,3,1,4,6,1,5,10,1,6,15,20,1,7,21,35,1,8,28,56,70,1,9,  &  ! n=1-9
        36,84,126,1,10,45,120,210,252,1,11,55,165,330,462,1,12,66,220,  &  ! n=9-12
        495,792,924,1,13,78,286,715,1287,1716,1,14,91,364,1001,2002,    &  ! n=12-14
        3003,3432,1,15,105,455,1365,3003,5005,6435,1,16,120,560,1820,   &  ! n=14-16
        4368,8008,11440,12870,1,17,136,680,2380,6188,12376,19448,24310, &  ! n=16-17
        1,18,153,816,3060,8568,18564,31824,43758,48620,1,19,171,969,    &  ! n=18-19
        3876,11628,27132,50388,75582,92378,1,20,190,1140,4845,15504,    &  ! n=19-20
        38760,77520,125970,167960,184756,1,21,210,1330,5985,20349,      &  ! n=20-21
        54264,116280,203490,293930,352716,1,22,231,1540,7315,26334,     &  ! n=21-22
        74613,170544,319770,497420,646646,705432,1,23,253,1771,8855,    &  ! n=22-23
        33649,100947,245157,490314,817190,1144066,1352078,1,24,276,     &  ! n=23-24
        2024,10626,42504,134596,346104,735471,1307504,1961256,2496144,  &  ! n=24
        2704156,1,25,300,2300,12650,53130,177100,480700,1081575,        &  ! n=24-25
        2042975,3268760,4457400,5200300,1,26,325,2600,14950,65780,      &  ! n=25-26
        230230,657800,1562275,3124550,5311735,7726160,9657700,10400600, &  ! n=26
        1,27,351,2925,17550,80730,296010,888030,2220075,4686825,        &  ! n=27
        8436285,13037895,17383860,20058300,1,28,378,3276,20475,98280,   &  ! n=27-28
        376740,1184040,3108105,6906900,13123110,21474180,30421755,      &  ! n=28
        37442160,40116600,1,29,406,3654,23751,118755,475020,1560780,    &  ! n=28-29
        4292145,10015005,20030010,34597290,51895935,67863915,77558760,  &  ! n=29
        1,30,435,4060,27405,142506,593775,2035800,5852925,14307150,     &  ! n=30
        30045015,54627300,86493225,119759850,145422675,155117520 /)        ! n=30
  ! Largest N at each K0 = min[K,N-K] for which recursive method
  ! below will not lead to overflow for 4-byte integers.
  INTEGER, PARAMETER :: K0MAX = 14
  INTEGER, PARAMETER :: NMAX_AT_K0(0:K0MAX) =                           &
    (/  HUGE(1),HUGE(1),46341,1626,338,140,82,58,46,39,35,33,31,30,30 /)
  ! Bad cases
  IF ((N .LT. K) .OR. (K .LT. 0)) THEN
    Cnk = -HUGE(N)
    !STOP 'ERROR: invalid binomial coefficient arguments (0 <= k <= n)'
    RETURN
  END IF
  ! Use (n k) = (n n-k)
  K0 = MIN(K,N-K)
  ! Use tabulated values for small n
  IF (N .LE. NMAX_TAB) THEN
    ! Parentheses are important: act as floor function
    J = ((N+1)/2) * ((N+2)/2) + K0 + 1
    Cnk = BC_ARRAY(J)
    RETURN
  ! Use recursion: (n k) = n/k (n-1 k-1)
  ! Note overflow for large enough k and/or n.
  ! Below check should avoid all potential overflow cases for
  ! 4-byte integers.
  ELSE IF (K0 .LE. K0MAX) THEN
    IF (N .LE. NMAX_AT_K0(K0)) THEN
      Cnk = 1
      DO J=1,K0
        Cnk = (Cnk*(N-K0+J))/J
      END DO
      RETURN
    END IF
  END IF
  ! Use LogGamma for large arguments.
  ! Seems to give full INTEGER*4 precision when not in an
  ! overflow condition.
  Cnk = NINT(EXP(LOG_GAMMA(N+1d0) - LOG_GAMMA(K+1d0) - LOG_GAMMA(N-K+1d0)))
END FUNCTION


!-----------------------------------------------------------------------
! Logarithm of binomial coefficient (n k) = n!/k!(n-k)! [double precision]
! 
ELEMENTAL FUNCTION LOG_BINOMIAL_COEFF(N,K) RESULT(lnCnk)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: N,K
  REAL*8 :: lnCnk
  INTEGER :: Cnk,K0,J
  INTEGER, PARAMETER :: NMAX_TAB = 30
  ! Largest N at each K0 = min[K,N-K] for which recursive method
  ! below will not lead to overflow for 4-byte integers.
  INTEGER, PARAMETER :: K0MAX = 14
  INTEGER, PARAMETER :: NMAX_AT_K0(0:K0MAX) =                           &
    (/  HUGE(1),HUGE(1),46341,1626,338,140,82,58,46,39,35,33,31,30,30 /)
  ! Bad cases
  IF ((N .LT. K) .OR. (K .LT. 0)) THEN
    lnCnk = -HUGE(N)
    !STOP 'ERROR: invalid binomial coefficient arguments (0 <= k <= n)'
    RETURN
  END IF
  ! Use (n k) = (n n-k)
  K0 = MIN(K,N-K)
  ! Get tabulated value from other routine
  IF (N .LE. NMAX_TAB) THEN
    lnCnk = LOG(1d0*BINOMIAL_COEFF(N,K))
    RETURN
  ! Use recursion: (n k) = n/k (n-1 k-1)
  ! Note overflow for large enough k and/or n.
  ! Below check should avoid all potential overflow cases for
  ! 4-byte integers.
  ELSE IF (K0 .LE. K0MAX) THEN
    IF (N .LE. NMAX_AT_K0(K0)) THEN
      Cnk = 1
      DO J=1,K0
        Cnk = (Cnk*(N-K0+J))/J
      END DO
      lnCnk = LOG(1d0*Cnk)
      RETURN
    END IF
  END IF
  ! Use LogGamma for large arguments.
  lnCnk = LOG_GAMMA(N+1d0) - LOG_GAMMA(K+1d0) - LOG_GAMMA(N-K+1d0)
END FUNCTION


!-----------------------------------------------------------------------
! Calculates W(x) to double precision, where W is the Lambert W
! function defined as the solution w = W(x) to the equation:
!     x = w e^w
! Returns the principal branch value (the larger of two solutions
! over -1/e < x < 0; there is only one solution for x > 0).  The
! W function is undefined for x < -1/e.
! 
! Valid input:  -1/e < x < \infty
! 
ELEMENTAL FUNCTION LAMBERTW(x) RESULT(w)
  IMPLICIT NONE
  REAL*8 :: w
  REAL*8, INTENT(IN) :: x
  REAL*8 :: epsk,zk,qk,wk,wk1,p,num,den,a,b,ia
  REAL*8, PARAMETER :: e    = 2.7182818284590452d0
  REAL*8, PARAMETER :: einv = 0.36787944117144232d0
  
  ! No solution for x < 1/e
  IF (x .LT. -einv) THEN
    !STOP 'Error in lambertw: argument must be larger than -EXP(-1)'
    w = -HUGE(w)
    RETURN
  END IF
  
  ! Could use Newton's or Halley's iteration method to find the
  ! solution, but the Lambert function has a faster method,
  ! Fritsch's iteration:
  !    W_{k+1} = W_k * (1 + eps_k)
  ! with
  !    eps_k   = z_k/(1 + W_k) * (q_k - z_k)/(q_k - 2*z_k)
  !    z_k     = ln(x/W_k) - W_k
  !    q_k     = 2 (1 + W_k) (1 + W_k + (2/3)*z_k)
  ! If eps_k is the error in W_k, then the error in W_{k+1} is of
  ! order (eps_k)^4, a faster convergent that Halley's method.
  ! For a first guess accurate to order O(10^-4), double precision
  ! can be achieved in only a single iteration.  We use estimates
  ! for the initial guess as determined by Veberic.
  ! 
  ! For further information, see:
  !   D. Veberic, arXiv:1003.1628.
  !   F.N. Fritsch, R.E. Shafer and W.P. Crowley,
  !       Commun. ACM 16, 123 (1973).
  
  ! Initial estimate by Veberic
  IF (x .LT. -0.32358170806015724d0) THEN
    ! branch point expansion
    p = SQRT(2*(1+e*x))
    wk = -1 + p*(1 + p*(-1d0/3 + p*(11d0/72 + p*(-43d0/540              &
              + p*(769d0/17280 + p*(-221d0/8505                         &
              + p*(680863d0/43545600 + p*(-1963d0/204120                &
              + p*(226287557d0/37623398400d0) ))))))))
  ELSE IF (x .LT. 0.14546954290661823d0) THEN
    ! rational fit
    num = x * (1 + x *                                                  &
              (5.931375839364438d0 + x *                                &
              (11.392205505329132d0+ x *                                &
              (7.338883399111118d0 + x * 0.6534490169919599d0) )))
    den = 1 + x *                                                       &
             (6.931373689597704d0 + x *                                 &
             (16.82349461388016d0 + x *                                 &
             (16.43072324143226d0 + x * 5.115235195211697d0) ))
    wk = num / den
  ELSE IF (x .LT. 8.706658967856612d0) THEN
    ! rational fit
    num = x * (1 + x *                                                  &
              (2.4450530707265568d0 + x *                               &
              (1.3436642259582265d0 + x *                               &
              (0.14844005539759195d0+ x * 8.047501729129999d-4) )))
    den = 1 + x *                                                       &
             (3.4447089864860025d0 + x *                                &
             (3.2924898573719523d0 + x *                                &
             (0.9164600188031222d0 + x * 0.05306864044833221d0) ))
    wk = num / den
  ELSE
    ! asymptotic expansion
    a = LOG(x)
    b = LOG(a)
    ia = 1/a
    wk = a - b + b * ia *                                               &
                (1 + ia *                                               &
                (0.5d0*(-2 + b) + ia *                                  &
                (1/6d0*(6 + b*(-9 + b*2)) + ia *                        &
                (1/12d0*(-12 + b*(36 + b*(-22 + b*3))) + ia *           &
                 1/60d0*(60 + b*(-300 + b*(350 + b*(-125 + b*12))))     &
                ))))
  END IF
  
  ! Special cases:
  ! For x equal to 0 or -1/e, the Fritsch iteration does
  ! not work as some of the terms go to infinity.  However,
  ! for x sufficiently near 0 or -1/e, the above first
  ! approximation is already nearly double precision.
  IF ((ABS(x) .LT. 1d-7) .OR. (x .LT. -einv+1d-6) ) THEN
    w = wk
    RETURN
  END IF
  
  ! Now apply Fritsch iteration
  wk1  = wk + 1
  zk   = LOG(x/wk) - wk
  qk   = 2 * wk1 * (wk1 + (2d0/3)*zk)
  epsk = zk * (qk - zk) / (wk1 * (qk - 2*zk))
  wk   = wk * (1 + epsk)
  ! In most cases, no further iterations will be necessary
  DO WHILE (ABS(epsk) .GT. 1d-5)
    wk1  = wk + 1
    zk   = LOG(x/wk) - wk
    qk   = 2 * wk1 * (wk1 + (2d0/3)*zk)
    epsk = zk * (qk - zk) / (wk1 * (qk - 2*zk))
    wk   = wk * (1 + epsk)
  END DO
  
  w = wk
  
END FUNCTION


!-----------------------------------------------------------------------
! Calculates W(x) to double precision, where W is the Lambert W
! function defined as the solution w = W(x) to the equation:
!     x = w e^w
! Returns the secondary branch value (the smaller of two solutions
! over -1/e < x < 0).  This branch is defined only for
! -1/e < x < 0.
! 
! Valid input:  -1/e < x < 0
! 
ELEMENTAL FUNCTION LAMBERTW2(x) RESULT(w2)
  IMPLICIT NONE
  REAL*8 :: w2
  REAL*8, INTENT(IN) :: x
  INTEGER :: k
  REAL*8 :: epsk,zk,qk,wk,wk1,p,num,den,a
  REAL*8, PARAMETER :: e    = 2.7182818284590452d0
  REAL*8, PARAMETER :: einv = 0.36787944117144232d0
  
  ! No solution for x < 1/e
  IF (x .LT. -einv) THEN
    !STOP 'Error in lambertw2: argument must be larger than -EXP(-1)'
    w2 = HUGE(w2)
    RETURN
  END IF
  
  ! No second branch for x > 0
  IF (x .GT. 0) THEN
    !STOP 'Error in lambertw2: argument must be smaller than 0'
    w2 = HUGE(w2)
    RETURN
  END IF
  
  ! Limit as x ->0 is -\infty
  IF (x .EQ. 0) THEN
    w2 = -HUGE(w2)
    RETURN
  END IF
  
  ! Could use Newton's or Halley's iteration method to find the
  ! solution, but the Lambert function has a faster method,
  ! Fritsch's iteration:
  !    W_{k+1} = W_k * (1 + eps_k)
  ! with
  !    eps_k   = z_k/(1 + W_k) * (q_k - z_k)/(q_k - 2*z_k)
  !    z_k     = ln(x/W_k) - W_k
  !    q_k     = 2 (1 + W_k) (1 + W_k + (2/3)*z_k)
  ! If eps_k is the error in W_k, then the error in W_{k+1} is of
  ! order (eps_k)^4, a faster convergent that Halley's method.
  ! For a first guess accurate to order O(10^-4), double precision
  ! can be achieved in only a single iteration.  We use estimates
  ! for the initial guess as determined by Veberic.
  ! 
  ! For further information, see:
  !   D. Veberic, arXiv:1003.1628.
  !   F.N. Fritsch, R.E. Shafer and W.P. Crowley,
  !       Commun. ACM 16, 123 (1973).
  
  ! Initial estimate by Veberic
  IF (x .LT. -0.30298541769d0) THEN
    ! branch point expansion
    p = -SQRT(2*(1+e*x))
    wk = -1 + p*(1 + p*(-1d0/3 + p*(11d0/72 + p*(-43d0/540              &
              + p*(769d0/17280 + p*(-221d0/8505                         &
              + p*(680863d0/43545600 + p*(-1963d0/204120                &
              + p*(226287557d0/37623398400d0) ))))))))
  ELSE IF (x .LT. -0.051012917658221676d0) THEN
    ! rational fit
    num = -7.814176723907436d0 + x *                                    &
            (253.88810188892484d0 + x * 657.9493176902304d0)
    den = 1 + x *                                                       &
             (-60.43958713690808d0+ x *                                 &
             (99.98567083107612d0 + x *                                 &
             (682.6073999909428d0 + x *                                 &
             (962.1784396969866d0 + x * 1477.9341280760887d0) )))
    wk = num / den
  ELSE
    ! continued logarithm
    a  = LOG(-x)
    wk = a
    DO k=1,9
      wk = a - LOG(-wk)
    END DO
  END IF
  
  ! Special cases:
  ! For x equal to -1/e, the Fritsch iteration does not
  ! not work as some of the terms go to infinity.  However,
  ! for x sufficiently near -1/e, the above first
  ! approximation is already nearly double precision.
  IF (x .LT. -einv+1d-6) THEN
    w2 = wk
    RETURN
  END IF
  
  ! Now apply Fritsch iteration
  wk1  = wk + 1
  zk   = LOG(x/wk) - wk
  qk   = 2 * wk1 * (wk1 + (2d0/3)*zk)
  epsk = zk * (qk - zk) / (wk1 * (qk - 2*zk))
  wk   = wk * (1 + epsk)
  ! In most cases, no further iterations will be necessary
  DO WHILE (ABS(epsk) .GT. 1d-5)
    wk1  = wk + 1
    zk   = LOG(x/wk) - wk
    qk   = 2 * wk1 * (wk1 + (2d0/3)*zk)
    epsk = zk * (qk - zk) / (wk1 * (qk - 2*zk))
    wk   = wk * (1 + epsk)
  END DO
  
  w2 = wk
  
END FUNCTION


!-----------------------------------------------------------------------
! Calculates W(x=e^ln(x)) to double precision, where W is the Lambert W
! function defined as the solution w = W(x) to the equation:
!     x = w e^w
! Returns the principal branch value (the larger of two solutions
! over -1/e < x < 0; there is only one solution for x > 0).  The
! W function is undefined for x < -1/e.
! 
! This version takes ln(x) as the input to allow for cases where x
! is large.
! 
! Valid input:  -\infty < lnx < \infty
! 
ELEMENTAL FUNCTION LAMBERTWLN(lnx) RESULT(w)
  IMPLICIT NONE
  REAL*8 :: w
  REAL*8, INTENT(IN) :: lnx
  REAL*8 :: epsk,zk,qk,wk,wk1,a,b,ia
  REAL*8, PARAMETER :: e    = 2.7182818284590452d0
  REAL*8, PARAMETER :: einv = 0.36787944117144232d0
  
  ! Here, we only calculate W(x) for very large x.  If x is a
  ! not very large, we use the lambertw routine.
  IF (lnx .LT. 300d0) THEN
    w = LAMBERTW(EXP(lnx))
    RETURN
  END IF
  
  ! Could use Newton's or Halley's iteration method to find the
  ! solution, but the Lambert function has a faster method,
  ! Fritsch's iteration:
  !    W_{k+1} = W_k * (1 + eps_k)
  ! with
  !    eps_k   = z_k/(1 + W_k) * (q_k - z_k)/(q_k - 2*z_k)
  !    z_k     = ln(x/W_k) - W_k
  !    q_k     = 2 (1 + W_k) (1 + W_k + (2/3)*z_k)
  ! If eps_k is the error in W_k, then the error in W_{k+1} is of
  ! order (eps_k)^4, a faster convergent that Halley's method.
  ! For a first guess accurate to order O(10^-4), double precision
  ! can be achieved in only a single iteration.  We use estimates
  ! for the initial guess as determined by Veberic.
  ! 
  ! For further information, see:
  !   D. Veberic, arXiv:1003.1628.
  !   F.N. Fritsch, R.E. Shafer and W.P. Crowley,
  !       Commun. ACM 16, 123 (1973).
  
  ! Initial estimate by Veberic
  ! asymptotic expansion
  a = lnx
  b = LOG(a)
  ia = 1/a
  wk = a - b + b * ia *                                                 &
              (1 + ia *                                                 &
              (0.5d0*(-2 + b) + ia *                                    &
              (1/6d0*(6 + b*(-9 + b*2)) + ia *                          &
              (1/12d0*(-12 + b*(36 + b*(-22 + b*3))) + ia *             &
               1/60d0*(60 + b*(-300 + b*(350 + b*(-125 + b*12))))       &
              ))))
  
  ! Now apply Fritsch iteration
  wk1  = wk + 1
  zk   = lnx - LOG(wk) - wk
  qk   = 2 * wk1 * (wk1 + (2d0/3)*zk)
  epsk = zk * (qk - zk) / (wk1 * (qk - 2*zk))
  wk   = wk * (1 + epsk)
  ! In most cases, no further iterations will be necessary
  DO WHILE (ABS(epsk) .GT. 1d-5)
    wk1  = wk + 1
    zk   = lnx - LOG(wk) - wk
    qk   = 2 * wk1 * (wk1 + (2d0/3)*zk)
    epsk = zk * (qk - zk) / (wk1 * (qk - 2*zk))
    wk   = wk * (1 + epsk)
  END DO
  
  w = wk
  
END FUNCTION


!=======================================================================
! PROBABILITY DISTRIBUTIONS
!=======================================================================

!----------UNIFORM DISTRIBUTION-----------------------------------------
! Uniform distribution has probability:
!   P(x|xmin,xmax) = 1/(xmax-xmin)    xmin < x < xmax
!                  = 0                otherwise
! The range of the distribution [xmin,xmax] has default [0,1].
!-----------------------------------------------------------------------

!----------------------------------------
! Probability for uniform distribution
!   x           Value distributed normally
!   xmin        Minimum of uniform values (optional, default is 0)
!   xmax        Maximum of uniform values (optional, default is 1)
ELEMENTAL FUNCTION UniformPDF(x,xmin,xmax) RESULT(pdf)
  IMPLICIT NONE
  REAL*8 :: pdf
  REAL*8, INTENT(IN) :: x
  REAL*8, INTENT(IN), OPTIONAL :: xmin,xmax
  REAL*8 :: xmin0,xmax0
  IF (PRESENT(xmin)) THEN
    xmin0 = xmin
  ELSE
    xmin0 = 0d0
  END IF
  IF (PRESENT(xmax)) THEN
    xmax0 = xmax
  ELSE
    xmax0 = 1d0
  END IF
  IF ((x .GE. xmin0) .AND. (x .LE. xmax0)) THEN
    pdf = 1d0 / (xmax0-xmin0)
  ELSE
    pdf = 0d0
  END IF
END FUNCTION UniformPDF

!----------------------------------------
! Log of probability for uniform distribution
!   x           Value distributed normally
!   xmin        Minimum of uniform values (optional, default is 0)
!   xmax        Maximum of uniform values (optional, default is 1)
ELEMENTAL FUNCTION UniformLogPDF(x,xmin,xmax) RESULT(lnpdf)
  IMPLICIT NONE
  REAL*8 :: lnpdf
  REAL*8, INTENT(IN) :: x
  REAL*8, INTENT(IN), OPTIONAL :: xmin,xmax
  REAL*8 :: xmin0,xmax0
  IF (PRESENT(xmin)) THEN
    xmin0 = xmin
  ELSE
    xmin0 = 0d0
  END IF
  IF (PRESENT(xmax)) THEN
    xmax0 = xmax
  ELSE
    xmax0 = 1d0
  END IF
  IF ((x .GE. xmin0) .AND. (x .LE. xmax0)) THEN
    lnpdf = LOG(1d0 / (xmax0-xmin0))
  ELSE
    lnpdf = -HUGE(1d0)
  END IF
END FUNCTION UniformLogPDF

!----------------------------------------
! CDF for uniform distribution
!   x           Value distributed normally
!   xmin        Minimum of uniform values (optional, default is 0)
!   xmax        Maximum of uniform values (optional, default is 1)
ELEMENTAL FUNCTION UniformCDF(x,xmin,xmax) RESULT(cdf)
  IMPLICIT NONE
  REAL*8 :: cdf
  REAL*8, INTENT(IN) :: x
  REAL*8, INTENT(IN), OPTIONAL :: xmin,xmax
  REAL*8 :: xmin0,xmax0
  IF (PRESENT(xmin)) THEN
    xmin0 = xmin
  ELSE
    xmin0 = 0d0
  END IF
  IF (PRESENT(xmax)) THEN
    xmax0 = xmax
  ELSE
    xmax0 = 1d0
  END IF
  IF (x .LE. xmin0) THEN
    cdf = 0d0
  ELSE IF (x .GE. xmax0) THEN
    cdf = 1d0
  ELSE
    cdf = (x-xmin0) / (xmax0-xmin0)
  END IF
END FUNCTION UniformCDF

!----------------------------------------
! 1-CDF for uniform distribution
!   x           Value distributed normally
!   xmin        Minimum of uniform values (optional, default is 0)
!   xmax        Maximum of uniform values (optional, default is 1)
ELEMENTAL FUNCTION Uniform1mCDF(x,xmin,xmax) RESULT(p)
  IMPLICIT NONE
  REAL*8 :: p
  REAL*8, INTENT(IN) :: x
  REAL*8, INTENT(IN), OPTIONAL :: xmin,xmax
  REAL*8 :: xmin0,xmax0
  IF (PRESENT(xmin)) THEN
    xmin0 = xmin
  ELSE
    xmin0 = 0d0
  END IF
  IF (PRESENT(xmax)) THEN
    xmax0 = xmax
  ELSE
    xmax0 = 1d0
  END IF
  IF (x .LE. xmin0) THEN
    p = 1d0
  ELSE IF (x .GE. xmax0) THEN
    p = 0d0
  ELSE
    p = (xmax-x) / (xmax0-xmin0)
  END IF
END FUNCTION Uniform1mCDF

!----------------------------------------
! Random number for uniform distribution
!   xmin        Minimum of uniform values (optional, default is 0)
!   xmax        Maximum of uniform values (optional, default is 1)
FUNCTION RandUniform(xmin,xmax) RESULT(x)
  IMPLICIT NONE
  REAL*8 :: x
  REAL*8, INTENT(IN), OPTIONAL :: xmin,xmax
  REAL*8 :: urand
  urand = rand()
  IF (PRESENT(xmin)) THEN
    IF (PRESENT(xmax)) THEN
      x = xmin + urand*(xmax-xmin)
    ELSE
      x = xmin + urand*(1d0-xmin)
    END IF
  ELSE
    IF (PRESENT(xmax)) THEN
      x = urand*xmax
    ELSE
      x = urand
    END IF
  END IF
END FUNCTION RandUniform


!----------NORMAL DISTRIBUTION------------------------------------------
! Normal distribution has probability:
!   P(x|mu,sigma) = EXP(-(x-mu)^2/(2 sigma^2)) / SQRT(2 PI sigma^2)
! The average mu has default 0 and the standard deviation sigma
! has default 1.
!-----------------------------------------------------------------------

!----------------------------------------
! Probability for normal distribution
!   x           Value distributed normally
!   mu          Average (optional, default is 0)
!   sigma       Standard deviation (optional, default is 1)
ELEMENTAL FUNCTION NormalPDF(x,mu,sigma) RESULT(pdf)
  IMPLICIT NONE
  REAL*8 :: pdf
  REAL*8, INTENT(IN) :: x
  REAL*8, INTENT(IN), OPTIONAL :: mu,sigma
  REAL*8 :: z
  REAL*8, PARAMETER :: SQRT2PI = 2.5066282746310005d0
  IF (PRESENT(mu)) THEN
    z = x - mu
  ELSE
    z = x
  END IF
  IF (PRESENT(sigma)) THEN
    pdf = EXP(-z**2/(2*sigma**2)) / (SQRT2PI*sigma)
  ELSE
    pdf = EXP(-z**2/2) / SQRT2PI
  END IF
END FUNCTION NormalPDF

!----------------------------------------
! Log of probability for normal distribution
!   x           Value distributed normally
!   mu          Average (optional, default is 0)
!   sigma       Standard deviation (optional, default is 1)
ELEMENTAL FUNCTION NormalLogPDF(x,mu,sigma) RESULT(lnpdf)
  IMPLICIT NONE
  REAL*8 :: lnpdf
  REAL*8, INTENT(IN) :: x
  REAL*8, INTENT(IN), OPTIONAL :: mu,sigma
  REAL*8 :: z
  REAL*8, PARAMETER :: LOGSQRT2PI = 0.91893853320467274d0
  IF (PRESENT(mu)) THEN
    z = x - mu
  ELSE
    z = x
  END IF
  IF (PRESENT(sigma)) THEN
    lnpdf = -z**2/(2*sigma**2) - LOGSQRT2PI - LOG(sigma)
  ELSE
    lnpdf = -z**2/2 - LOGSQRT2PI
  END IF
END FUNCTION NormalLogPDF

!----------------------------------------
! CDF for normal distribution
!   x           Value distributed normally
!   mu          Average (optional, default is 0)
!   sigma       Standard deviation (optional, default is 1)
ELEMENTAL FUNCTION NormalCDF(x,mu,sigma) RESULT(cdf)
  IMPLICIT NONE
  REAL*8 :: cdf
  REAL*8, INTENT(IN) :: x
  REAL*8, INTENT(IN), OPTIONAL :: mu,sigma
  REAL*8 :: z
  REAL*8, PARAMETER :: SQRT2 = 1.4142135623730950d0
  IF (PRESENT(mu)) THEN
    z = x - mu
  ELSE
    z = x
  END IF
  IF (PRESENT(sigma)) THEN
    z = z / (SQRT2*sigma)
  ELSE
    z = z / SQRT2
  END IF
  cdf = 0.5d0 * ERFC(-z)
END FUNCTION NormalCDF

!----------------------------------------
! 1-CDF for normal distribution
!   x           Value distributed normally
!   mu          Average (optional, default is 0)
!   sigma       Standard deviation (optional, default is 1)
ELEMENTAL FUNCTION Normal1mCDF(x,mu,sigma) RESULT(p)
  IMPLICIT NONE
  REAL*8 :: p
  REAL*8, INTENT(IN) :: x
  REAL*8, INTENT(IN), OPTIONAL :: mu,sigma
  REAL*8 :: z
  REAL*8, PARAMETER :: SQRT2 = 1.4142135623730950d0
  IF (PRESENT(mu)) THEN
    z = x - mu
  ELSE
    z = x
  END IF
  IF (PRESENT(sigma)) THEN
    z = z / (SQRT2*sigma)
  ELSE
    z = z / SQRT2
  END IF
  p = 0.5d0 * ERFC(z)
END FUNCTION Normal1mCDF

!----------------------------------------
! Random number for normal distribution
!   mu          Average (optional, default is 0)
!   sigma       Standard deviation (optional, default is 1)
FUNCTION RandNormal(mu,sigma) RESULT(x)
  IMPLICIT NONE
  REAL*8 :: x
  REAL*8, INTENT(IN), OPTIONAL :: mu,sigma
  REAL*8 :: u,v,r2,w,urand(2),rand1
  REAL*8, SAVE :: rand2
  LOGICAL, SAVE :: needs_calc = .TRUE.
  ! Algorithm calculates two random normal numbers at a time.
  ! Store second for next call to this function.
  IF (needs_calc) THEN
    r2 = 1d0
    DO WHILE (r2 .GE. 1d0)
      urand(1) = rand()
      urand(2) = rand()
      u = 2d0*urand(1) - 1d0
      v = 2d0*urand(2) - 1d0
      r2 = u**2 + v**2
    END DO
    w = SQRT(-2d0 * LOG(r2) / r2)
    rand1 = u * w
    rand2 = v * w
    x = rand1
    needs_calc = .FALSE.
  ELSE
    x = rand2
    needs_calc = .TRUE.
  END IF
  IF (PRESENT(sigma)) x = sigma*x
  IF (PRESENT(mu)) x = mu + x
END FUNCTION RandNormal


!----------LOG NORMAL DISTRIBUTION--------------------------------------
! Normal distribution has probability:
!   P(x|mu,sigma) = EXP(-(ln(x)-mu)^2/(2 sigma^2)) / SQRT(2 PI sigma^2) / x
! The location mu has default 0 and the scale sigma has default 1.
! Note the parameters mu and sigma are not the mean and standard
! deviation.
!   mean:      muL      = exp(mu + sigma^2/2)
!   variance:  sigmaL^2 = (exp(sigma^2)-1) exp(2*mu+sigma^2)
!       =>     mu      = -1/2 log[(1+sigmaL^2/muL^2)/muL^2]
!              sigma^2 = log[1+sigmaL^2/muL^2]
! A normal distribution with parameters sigmaN << muN is asymptotically
! similar to a log normal distribution with parameters:
!   mu    = ln(muN)
!   sigma = sigmaN/muN
! Specifically, setting muL = muN and sigmaL = sigmaN:
!   mu      = log(muN) - 1/2 log(1+sigmaN^2/muN^2)
!           = log(muN) + O[sigmaN^2/muN^2]
!   sigma^2 = log(1+sigmaN^2/muN^2)
!           = sigmaN^2/muN^2 (1 + O[sigmaN^2/muN^2])
!-----------------------------------------------------------------------

!----------------------------------------
! Probability for log normal distribution
!   x           Value distributed log-normally
!   mu          Location (optional, default is 0)
!   sigma       Scale (optional, default is 1)
ELEMENTAL FUNCTION LogNormalPDF(x,mu,sigma) RESULT(pdf)
  IMPLICIT NONE
  REAL*8 :: pdf
  REAL*8, INTENT(IN) :: x
  REAL*8, INTENT(IN), OPTIONAL :: mu,sigma
  REAL*8 :: z
  REAL*8, PARAMETER :: SQRT2PI = 2.5066282746310005d0
  ! Bad case
  IF (x .LE. 0d0) THEN
    pdf = 0d0
    RETURN
  END IF
  IF (PRESENT(mu)) THEN
    z = LOG(x) - mu
  ELSE
    z = LOG(x)
  END IF
  IF (PRESENT(sigma)) THEN
    pdf = EXP(-z**2/(2*sigma**2)) / (SQRT2PI*sigma*x)
  ELSE
    pdf = EXP(-z**2/2) / (SQRT2PI*x)
  END IF
END FUNCTION LogNormalPDF

!----------------------------------------
! Log of probability for log normal distribution
!   x           Value distributed normally
!   mu          Average (optional, default is 0)
!   sigma       Standard deviation (optional, default is 1)
ELEMENTAL FUNCTION LogNormalLogPDF(x,mu,sigma) RESULT(lnpdf)
  IMPLICIT NONE
  REAL*8 :: lnpdf
  REAL*8, INTENT(IN) :: x
  REAL*8, INTENT(IN), OPTIONAL :: mu,sigma
  REAL*8 :: z
  REAL*8, PARAMETER :: LOGSQRT2PI = 0.91893853320467274d0
  ! Bad case
  IF (x .LE. 0d0) THEN
    lnpdf = -HUGE(1d0)
    RETURN
  END IF
  IF (PRESENT(mu)) THEN
    z = LOG(x) - mu
  ELSE
    z = LOG(x)
  END IF
  IF (PRESENT(sigma)) THEN
    lnpdf = -z**2/(2*sigma**2) - LOGSQRT2PI - LOG(x*sigma)
  ELSE
    lnpdf = -z**2/2 - LOGSQRT2PI - LOG(x)
  END IF
END FUNCTION LogNormalLogPDF

!----------------------------------------
! CDF for log normal distribution
!   x           Value distributed log-normally
!   mu          Location (optional, default is 0)
!   sigma       Scale (optional, default is 1)
ELEMENTAL FUNCTION LogNormalCDF(x,mu,sigma) RESULT(cdf)
  IMPLICIT NONE
  REAL*8 :: cdf
  REAL*8, INTENT(IN) :: x
  REAL*8, INTENT(IN), OPTIONAL :: mu,sigma
  REAL*8 :: z
  REAL*8, PARAMETER :: SQRT2 = 1.4142135623730950d0
  ! Bad case
  IF (x .LE. 0d0) THEN
    cdf = 0d0
    RETURN
  END IF
  IF (PRESENT(mu)) THEN
    z = LOG(x) - mu
  ELSE
    z = LOG(x)
  END IF
  IF (PRESENT(sigma)) THEN
    z = LOG(z) / (SQRT2*sigma)
  ELSE
    z = LOG(z) / SQRT2
  END IF
  cdf = 0.5d0 * ERFC(-z)
END FUNCTION LogNormalCDF

!----------------------------------------
! 1-CDF for log normal distribution
!   x           Value distributed log-normally
!   mu          Location (optional, default is 0)
!   sigma       Scale (optional, default is 1)
ELEMENTAL FUNCTION LogNormal1mCDF(x,mu,sigma) RESULT(p)
  IMPLICIT NONE
  REAL*8 :: p
  REAL*8, INTENT(IN) :: x
  REAL*8, INTENT(IN), OPTIONAL :: mu,sigma
  REAL*8 :: z
  REAL*8, PARAMETER :: SQRT2 = 1.4142135623730950d0
  ! Bad case
  IF (x .LE. 0d0) THEN
    p = 1d0
    RETURN
  END IF
  IF (PRESENT(mu)) THEN
    z = LOG(x) - mu
  ELSE
    z = LOG(x)
  END IF
  IF (PRESENT(sigma)) THEN
    z = LOG(z) / (SQRT2*sigma)
  ELSE
    z = LOG(z) / SQRT2
  END IF
  p = 0.5d0 * ERFC(z)
END FUNCTION LogNormal1mCDF

!----------------------------------------
! Random number for log normal distribution
!   mu          Location (optional, default is 0)
!   sigma       Scale (optional, default is 1)
FUNCTION RandLogNormal(mu,sigma) RESULT(x)
  IMPLICIT NONE
  REAL*8 :: x
  REAL*8, INTENT(IN), OPTIONAL :: mu,sigma
  REAL*8 :: xn
  xn = RandNormal()
  IF (PRESENT(sigma)) xn = sigma*xn
  IF (PRESENT(mu)) xn = mu + xn
  x = EXP(xn)
END FUNCTION RandLogNormal


!----------EXPONENTIAL DISTRIBUTION-------------------------------------
! Exponential distribution has probability:
!   P(x|xmin,xmax,xs) = A EXP(-x/xs)    xmin < x < xmax
!                     = 0               otherwise
! where the constant 'A' is a normalization factor.
! The range of the distribution [xmin,xmax] has default [0,infinity].
! The scale of the exponential xs has default 1.
!-----------------------------------------------------------------------

!----------------------------------------
! Probability for exponential distribution
!   x           Value distributed exponentially
!   xmin        Minimum value (optional, default is 0)
!   xmax        Maximum value (optional, default is infinity)
!   xs          The exponential scale (optional, default is 1)
ELEMENTAL FUNCTION ExponentialPDF(x,xmin,xmax,xs) RESULT(pdf)
  IMPLICIT NONE
  REAL*8 :: pdf
  REAL*8, INTENT(IN) :: x
  REAL*8, INTENT(IN), OPTIONAL :: xmin,xmax,xs
  REAL*8 :: xmin0,xmax0
  IF (PRESENT(xmin)) THEN
    xmin0 = xmin
  ELSE
    xmin0 = 0d0
  END IF
  IF (PRESENT(xmax)) THEN
    xmax0 = xmax
  ELSE
    xmax0 = HUGE(1d0)
  END IF
  IF ((x .LT. xmin0) .OR. (x .GT. xmax0)) THEN
    pdf = 0d0
    RETURN
  END IF
  IF (PRESENT(xs)) THEN
    IF (xs .EQ. 0d0) THEN
      pdf = 1d0 / (xmax0-xmin0)
    ELSE
      pdf = EXP((xmin0-x)/xs) / (1d0 - EXP((xmin0-xmax0)/xs)) / xs
    END IF
  ELSE
    pdf = EXP(xmin0-x) / (1d0 - EXP(xmin0-xmax0))
  END IF
END FUNCTION ExponentialPDF

!----------------------------------------
! Log of probability for exponential distribution
!   x           Value distributed exponentially
!   xmin        Minimum value (optional, default is 0)
!   xmax        Maximum value (optional, default is infinity)
!   xs          The exponential scale (optional, default is 1)
ELEMENTAL FUNCTION ExponentialLogPDF(x,xmin,xmax,xs) RESULT(lnpdf)
  IMPLICIT NONE
  REAL*8 :: lnpdf
  REAL*8, INTENT(IN) :: x
  REAL*8, INTENT(IN), OPTIONAL :: xmin,xmax,xs
  REAL*8 :: xmin0,xmax0
  IF (PRESENT(xmin)) THEN
    xmin0 = xmin
  ELSE
    xmin0 = 0d0
  END IF
  IF (PRESENT(xmax)) THEN
    xmax0 = xmax
  ELSE
    xmax0 = HUGE(1d0)
  END IF
  IF ((x .LT. xmin0) .OR. (x .GT. xmax0)) THEN
    lnpdf = -HUGE(1d0)
    RETURN
  END IF
  IF (PRESENT(xs)) THEN
    IF (xs .EQ. 0d0) THEN
      lnpdf = LOG(1d0 / (xmax0-xmin0))
    ELSE
      lnpdf = (xmin0-x)/xs - LOGp1(-EXP((xmin0-xmax0)/xs)) - LOG(xs)
    END IF
  ELSE
    lnpdf = (xmin0-x) - LOGp1(-EXP(xmin0-xmax0))
  END IF
END FUNCTION ExponentialLogPDF

!----------------------------------------
! CDF for exponential distribution
!   x           Value distributed exponentially
!   xmin        Minimum value (optional, default is 0)
!   xmax        Maximum value (optional, default is infinity)
!   xs          The exponential scale (optional, default is 1)
ELEMENTAL FUNCTION ExponentialCDF(x,xmin,xmax,xs) RESULT(cdf)
  IMPLICIT NONE
  REAL*8 :: cdf
  REAL*8, INTENT(IN) :: x
  REAL*8, INTENT(IN), OPTIONAL :: xmin,xmax,xs
  REAL*8 :: xmin0,xmax0
  IF (PRESENT(xmin)) THEN
    xmin0 = xmin
  ELSE
    xmin0 = 0d0
  END IF
  IF (PRESENT(xmax)) THEN
    xmax0 = xmax
  ELSE
    xmax0 = HUGE(1d0)
  END IF
  IF (x .LE. xmin0) THEN
    cdf = 0d0
  ELSE IF (x .GE. xmax0) THEN
    cdf = 1d0
  ELSE
    IF (PRESENT(xs)) THEN
      IF (xs .EQ. 0d0) THEN
        cdf = (x-xmin0) / (xmax0-xmin0)
      ELSE
        !cdf = (1d0 - EXP((xmin0-x)/xs)) / (1d0 - EXP((xmin0-xmax0)/xs))
        cdf = EXPm1((xmin0-x)/xs) / EXPm1((xmin0-xmax0)/xs)
      END IF
    ELSE
      !cdf = (1d0 - EXP(xmin0-x)) / (1d0 - EXP(xmin0-xmax0))
      cdf = EXPm1(xmin0-x) / EXPm1(xmin0-xmax0)
    END IF
  END IF
END FUNCTION ExponentialCDF

!----------------------------------------
! 1-CDF for exponential distribution
!   x           Value distributed exponentially
!   xmin        Minimum value (optional, default is 0)
!   xmax        Maximum value (optional, default is infinity)
!   xs          The exponential scale (optional, default is 1)
ELEMENTAL FUNCTION Exponential1mCDF(x,xmin,xmax,xs) RESULT(p)
  IMPLICIT NONE
  REAL*8 :: p
  REAL*8, INTENT(IN) :: x
  REAL*8, INTENT(IN), OPTIONAL :: xmin,xmax,xs
  REAL*8 :: xmin0,xmax0
  IF (PRESENT(xmin)) THEN
    xmin0 = xmin
  ELSE
    xmin0 = 0d0
  END IF
  IF (PRESENT(xmax)) THEN
    xmax0 = xmax
  ELSE
    xmax0 = HUGE(1d0)
  END IF
  IF (x .LE. xmin0) THEN
    p = 1d0
  ELSE IF (x .GE. xmax0) THEN
    p = 0d0
  ELSE
    IF (PRESENT(xs)) THEN
      IF (xs .EQ. 0d0) THEN
        p = (xmax0-x) / (xmax0-xmin0)
      ELSE
        !p = (EXP((xmax0-x)/xs) - 1d0) / (EXP((xmax0-xmin0)/xs) - 1d0)
        p = EXPm1((xmax0-x)/xs) / EXPm1((xmax0-xmin0)/xs)
      END IF
    ELSE
      !p = (EXP(xmax0-x) - 1d0) / (EXP(xmax0-xmin0) - 1d0)
      p = EXPm1(xmax0-x) / EXPm1(xmax0-xmin0)
    END IF
  END IF
END FUNCTION Exponential1mCDF

!----------------------------------------
! Random number for exponential distribution
!   xmin        Minimum value (optional, default is 0)
!   xmax        Maximum value (optional, default is infinity)
!   xs          The exponential scale (optional, default is 1)
FUNCTION RandExponential(xmin,xmax,xs) RESULT(x)
  IMPLICIT NONE
  REAL*8 :: x
  REAL*8, INTENT(IN), OPTIONAL :: xmin,xmax,xs
  REAL*8 :: xmin0,xmax0
  REAL*8 :: urand
  urand = rand()
  IF (PRESENT(xmin)) THEN
    xmin0 = xmin
  ELSE
    xmin0 = 0d0
  END IF
  IF (PRESENT(xmax)) THEN
    xmax0 = xmax
  ELSE
    xmax0 = HUGE(1d0)
  END IF
  IF (PRESENT(xs)) THEN
    IF (xs .EQ. 0d0) THEN
      x = xmin0 + urand*(xmax0-xmin0)
    ELSE
      x = xmin0 - xs*LOG(1d0 - urand*(1d0-EXP((xmin0-xmax0)/xs)))
    END IF
  ELSE
    x = xmin0 - LOG(1d0 - urand*(1d0-EXP(xmin0-xmax0)))
  END IF
END FUNCTION RandExponential


!----------CHI-SQUARE DISTRIBUTION--------------------------------------
! Chi-square distribution has probability:
!   P(x|k) = x^(k/2 - 1) EXP(-x/2) / 2^(k/2) / GAMMA(k/2)
!-----------------------------------------------------------------------

!----------------------------------------
! Probability for chi-square distribution
!   x           Value chi-square distributed
!   k           Degrees of freedom
ELEMENTAL FUNCTION ChiSquarePDF(x,k) RESULT(pdf)
  IMPLICIT NONE
  REAL*8 :: pdf
  REAL*8, INTENT(IN) :: x
  INTEGER, INTENT(IN) :: k
  IF (x .GE. 0d0) THEN
    pdf = x**(0.5d0*k-1) * EXP(-0.5d0*x) / (2**(0.5d0*k)*GAMMA(0.5d0*k))
  ELSE
    pdf = 0d0
  END IF
END FUNCTION ChiSquarePDF

!----------------------------------------
! Log of probability for chi-square distribution
!   x           Value chi-square distributed
!   k           Degrees of freedom
ELEMENTAL FUNCTION ChiSquareLogPDF(x,k) RESULT(lnpdf)
  IMPLICIT NONE
  REAL*8 :: lnpdf
  REAL*8, INTENT(IN) :: x
  INTEGER, INTENT(IN) :: k
  REAL*8, PARAMETER :: LOG2 = 0.69314718055994531d0
  IF (x .GE. 0d0) THEN
    lnpdf = (0.5d0*k-1)*LOG(x) - 0.5d0*x                                &
            - (0.5d0*k)*LOG2 - LOG_GAMMA(0.5d0*k)
  ELSE
    lnpdf = 0d0
  END IF
END FUNCTION ChiSquareLogPDF

!----------------------------------------
! CDF for chi-square distribution
!   x           Value chi-square distributed
!   k           Degrees of freedom
ELEMENTAL FUNCTION ChiSquareCDF(x,k) RESULT(cdf)
  IMPLICIT NONE
  REAL*8 :: cdf
  REAL*8, INTENT(IN) :: x
  INTEGER, INTENT(IN) :: k
  IF (x .GE. 0d0) THEN
    ! Regularized gamma function P(k/2,x/2)
    cdf = GAMMA_P(0.5d0*k,0.5d0*x)
  ELSE
    cdf = 0d0
  END IF
END FUNCTION ChiSquareCDF

!----------------------------------------
! 1-CDF for chi-square distribution
!   x           Value chi-square distributed
!   k           Degrees of freedom
ELEMENTAL FUNCTION ChiSquare1mCDF(x,k) RESULT(p)
  IMPLICIT NONE
  REAL*8 :: p
  REAL*8, INTENT(IN) :: x
  INTEGER, INTENT(IN) :: k
  IF (x .GE. 0d0) THEN
    ! Regularized gamma function Q(k/2,x/2)
    p = GAMMA_Q(0.5d0*k,0.5d0*x)
  ELSE
    p = 1d0
  END IF
END FUNCTION ChiSquare1mCDF

!----------------------------------------
! Random number for normal distribution
!   k           Degrees of freedom
FUNCTION RandChiSquare(k) RESULT(x)
  IMPLICIT NONE
  REAL*8 :: x
  INTEGER, INTENT(IN) :: k
  REAL*8 :: u
  u = rand()
  !x = REDUCEDGAMMAINV(0.5d0*k,u)
  x = -1d0
  WRITE(*,*) 'ERROR: RandChiSquare not implemented (REDUCEDGAMMAINV)'
  STOP
END FUNCTION RandChiSquare


!----------POISSON DISTRIBUTION-----------------------------------------
! Poisson distribution has probability:
!   P(N|mu) = EXP(-mu) mu^N / N!
!-----------------------------------------------------------------------

!----------------------------------------
! Probability for Poisson distribution
!   N           Number Poisson distributed
!   mu          Average
ELEMENTAL FUNCTION PoissonPDF(N,mu) RESULT(pdf)
  IMPLICIT NONE
  REAL*8 :: pdf
  INTEGER, INTENT(IN) :: N
  REAL*8, INTENT(IN) :: mu
  ! Bad cases
  IF (N .LT. 0) THEN
    pdf = 0d0
    RETURN
  ELSE IF (mu .LT. 0d0) THEN
    pdf = -HUGE(1d0)
    RETURN
  END IF
  ! Catch special case mu = 0
  IF (mu .EQ. 0d0) THEN
    IF (N .EQ. 0) THEN
      pdf = 1d0
    ELSE
      pdf = 0d0
    END IF
    RETURN
  END IF
  ! General case
  pdf = EXP(-mu + N*LOG(mu) - LOG_GAMMA(N+1d0))
END FUNCTION PoissonPDF

!----------------------------------------
! Log of probability for Poisson distribution
!   N           Number Poisson distributed
!   mu          Average
ELEMENTAL FUNCTION PoissonLogPDF(N,mu) RESULT(lnpdf)
  IMPLICIT NONE
  REAL*8 :: lnpdf
  INTEGER, INTENT(IN) :: N
  REAL*8, INTENT(IN) :: mu
  ! Bad cases
  IF ((N .LT. 0) .OR. (mu .LT. 0d0)) THEN
    lnpdf = -HUGE(1d0)
    RETURN
  END IF
  ! Catch special case mu = 0
  IF (mu .EQ. 0d0) THEN
    IF (N .EQ. 0) THEN
      lnpdf = 0d0
    ELSE
      lnpdf = -HUGE(1d0)
    END IF
    RETURN
  END IF
  ! General case
  lnpdf = -mu + N*LOG(mu) - LOG_GAMMA(N+1d0)
END FUNCTION PoissonLogPDF

!----------------------------------------
! CDF for Poisson distribution
!   N           Number Poisson distributed
!   mu          Average
ELEMENTAL FUNCTION PoissonCDF(N,mu) RESULT(cdf)
  IMPLICIT NONE
  REAL*8 :: cdf
  INTEGER, INTENT(IN) :: N
  REAL*8, INTENT(IN) :: mu
  ! Bad cases
  IF (N .LT. 0) THEN
    cdf = 0d0
    RETURN
  ELSE IF (mu .LT. 0d0) THEN
    cdf = -HUGE(1d0)
    RETURN
  END IF
  ! Catch special case mu = 0
  IF (mu .EQ. 0d0) THEN
    cdf = 1d0
    RETURN
  END IF
  ! General case
  ! Regularized gamma function Q(N+1,mu)
  cdf = GAMMA_Q(N+1d0,mu)
END FUNCTION PoissonCDF

!----------------------------------------
! 1-CDF for Poisson distribution
!   N           Number Poisson distributed
!   mu          Average
ELEMENTAL FUNCTION Poisson1mCDF(N,mu) RESULT(p)
  IMPLICIT NONE
  REAL*8 :: p
  INTEGER, INTENT(IN) :: N
  REAL*8, INTENT(IN) :: mu
  ! Bad cases
  IF (N .LT. 0) THEN
    p = 1d0
    RETURN
  ELSE IF (mu .LT. 0d0) THEN
    p = -HUGE(1d0)
    RETURN
  END IF
  ! Catch special case mu = 0
  IF (mu .EQ. 0d0) THEN
    p = 0d0
    RETURN
  END IF
  ! General case
  ! Regularized gamma function P(N+1,mu)
  p = GAMMA_P(N+1d0,mu)
END FUNCTION Poisson1mCDF

!----------------------------------------
! Random number for Poisson distribution
!   mu          Average
FUNCTION RandPoisson(mu) RESULT(N)
  IMPLICIT NONE
  INTEGER :: N
  REAL*8, INTENT(IN) :: mu
  REAL*8 :: alpha,beta,k,c,u,u2,x,xsum,z
  REAL*8, PARAMETER :: PI = 3.1415926535897932d0
  ! Small mu case: simply iterate over terms.
  IF (mu .LE. 40d0) THEN
    u = rand()
    N = 0
    x = EXP(-mu)
    xsum = x
    DO WHILE (xsum .LT. u)
      N = N+1
      x = mu * x / N
      xsum = xsum + x
      ! Bad case: start over (can occur for u extremely
      ! close to 1).  Note 1-CDF(102|mu=40) = 7e-17
      ! => remaining terms too small.
      IF (N .GT. 102) THEN
        u = rand()
        N = 0
        x = EXP(-mu)
        xsum = x
      END IF
    END DO
    RETURN
  END IF
  ! Larger mu case: use rejection technique.
  ! See Atkinson, Appl. Statist. 28, 29 (1979).
  beta  = PI / SQRT(3*mu)
  alpha = beta*mu
  c = 0.767d0 - 3.36d0/mu
  k = LOG(c/beta) - mu
  DO WHILE (.TRUE.)
    u = rand()
    x = (alpha - LOG((1-u)/u)) / beta
    N = NINT(x)
    IF (N .LT. 0) CYCLE
    u2 = rand()
    z = alpha - beta*x
    IF (z+LOG(u2/(1+EXP(z))**2) .LE. k+N*LOG(mu)-LOG_GAMMA(N+1d0)) THEN
      RETURN
    END IF
  END DO
END FUNCTION RandPoisson


!----------BINOMIAL DISTRIBUTION----------------------------------------
! Binomial distribution has probability:
!   P(k|p,N) = [N!/(k! (N-k)!)] p^k (1-p)^(N-k)
!-----------------------------------------------------------------------

!----------------------------------------
! Probability for binomial distribution
!   k           Number of successful trials
!   p           Probability of success for each trial
!   N           Number of trials
ELEMENTAL FUNCTION BinomialPDF(k,p,N) RESULT(pdf)
  IMPLICIT NONE
  REAL*8 :: pdf
  INTEGER, INTENT(IN) :: k,N
  REAL*8, INTENT(IN) :: p
  INTEGER :: k0
  REAL*8 :: lnCnk
  ! Bad cases
  IF ((k .LT. 0) .OR. (k .GT. N)) THEN
    !pdf = -HUGE(pdf)
    pdf = 0d0
    RETURN
  END IF
  ! Special cases
  IF (p .LE. 0d0) THEN
    IF (k .EQ. 0) THEN
      pdf = 1d0
    ELSE
      pdf = 0d0
    END IF
    RETURN
  ELSE IF (p .GE. 1d0) THEN
    IF (k .EQ. N) THEN
      pdf = 1d0
    ELSE
      pdf = 0d0
    END IF
    RETURN
  END IF
  ! General case
  k0 = MIN(k,N-k)
  IF (N .LE. 30) THEN
    ! quick, exact for given condition
    pdf = BINOMIAL_COEFF(N,k) * p**k * (1-p)**(N-k)
  !ELSE IF (N .LE. 1000) THEN
  !  lnCnk = LOG_GAMMA(N+1d0) - LOG_GAMMA(k+1d0) - LOG_GAMMA((N-k)+1d0)
  !  pdf   = EXP(lnCnk) * p**k * (1-p)**(N-k)
  ! Overflow/underflow issues for N ~ O(10^3)
  ELSE
    lnCnk = LOG_GAMMA(N+1d0) - LOG_GAMMA(k+1d0) - LOG_GAMMA((N-k)+1d0)
    pdf   = EXP(lnCnk + k*LOG(p) + (N-k)*LOG(1-p))
  END IF
END FUNCTION BinomialPDF

!----------------------------------------
! Probability for binomial distribution
!   k           Number of successful trials
!   p           Probability of success for each trial
!   N           Number of trials
ELEMENTAL FUNCTION BinomialLogPDF(k,p,N) RESULT(lnpdf)
  IMPLICIT NONE
  REAL*8 :: lnpdf
  INTEGER, INTENT(IN) :: k,N
  REAL*8, INTENT(IN) :: p
  INTEGER :: k0
  REAL*8 :: lnCnk
  ! Bad cases
  IF ((k .LT. 0) .OR. (k .GT. N)) THEN
    lnpdf = -HUGE(1d0)
    RETURN
  END IF
  ! Special cases
  IF (p .LE. 0d0) THEN
    IF (k .EQ. 0) THEN
      lnpdf = 0d0
    ELSE
      lnpdf = -HUGE(1d0)
    END IF
    RETURN
  ELSE IF (p .GE. 1d0) THEN
    IF (k .EQ. N) THEN
      lnpdf = 0d0
    ELSE
      lnpdf = -HUGE(1d0)
    END IF
    RETURN
  END IF
  ! General case
  k0 = MIN(k,N-k)
  IF (N .LE. 30) THEN
    ! quick, exact for given condition
    lnCnk = LOG_BINOMIAL_COEFF(N,k)
  ELSE
    lnCnk = LOG_GAMMA(N+1d0) - LOG_GAMMA(k+1d0) - LOG_GAMMA((N-k)+1d0)
  END IF
  lnpdf = lnCnk + k*LOG(p) +(N-k)*LOG(1-p)
END FUNCTION BinomialLogPDF

!----------------------------------------
! CDF for binomial distribution
!   k           Number of successful trials
!   p           Probability of success for each trial
!   N           Number of trials
! NOTE: Not optimized for large N (scales as sqrt(N*p*q)).
ELEMENTAL FUNCTION BinomialCDF(k,p,N) RESULT(cdf)
  IMPLICIT NONE
  REAL*8 :: cdf
  INTEGER, INTENT(IN) :: k,N
  REAL*8, INTENT(IN) :: p
  REAL*8 :: ltsum,gtsum,pdf
  ! Note cdf = BetaRegularized[1-p,N-k,k+1]
  CALL BinomialSums(k,p,N,ltsum,gtsum,pdf)
  cdf = ltsum + pdf
END FUNCTION BinomialCDF

!----------------------------------------
! 1-CDF for binomial distribution
!   k           Number of successful trials
!   p           Probability of success for each trial
!   N           Number of trials
! NOTE: Not optimized for large N (scales as sqrt(N*p*q)).
ELEMENTAL FUNCTION Binomial1mCDF(k,p,N) RESULT(sf)
  IMPLICIT NONE
  REAL*8 :: sf
  INTEGER, INTENT(IN) :: k,N
  REAL*8, INTENT(IN) :: p
  REAL*8 :: ltsum,gtsum,pdf
  ! Note cdf = BetaRegularized[1-p,N-k,k+1]
  CALL BinomialSums(k,p,N,ltsum,gtsum,pdf)
  sf = gtsum
END FUNCTION Binomial1mCDF

!----------------------------------------
! Calculates various sums for binomial distribution.
! Utility routine for CDF and 1-CDF calculations.
!   k           Number of successful trials
!   p           Probability of success for each trial
!   N           Number of trials
! Output:
!   ltsum       Sum of PDFs for j=0,..,k-1
!   gtsum       Sum of PDFs for j=k+1,..,N
!   pdf         The PDF for k
! NOTE: Not optimized for large N (scales as sqrt(N*p*q)).
ELEMENTAL SUBROUTINE BinomialSums(k,p,N,ltsum,gtsum,pdf)
  IMPLICIT NONE
  REAL*8 :: cdf
  INTEGER, INTENT(IN) :: k,N
  REAL*8, INTENT(IN) :: p
  REAL*8, INTENT(OUT) :: ltsum,gtsum,pdf
  LOGICAL :: flipped
  INTEGER :: k0,I
  REAL*8 :: p0,q0,mu,xI,tmp
  
  ! Bad cases
  IF (k .LT. 0) THEN
    ltsum = 0d0
    gtsum = 1d0
    pdf   = 0d0
    RETURN
  ELSE IF (k .GT. N) THEN
    ltsum = 1d0
    gtsum = 0d0
    pdf   = 0d0
    RETURN
  END IF
  
  ! Can write CDF(k|p,N) = CDF(N-k|1-p,N) = CDF(k0,p0,N)
  ! Choose smaller of p or 1-p to work with
  IF (p .GT. 0.5d0) THEN
    p0 = 1-p
    q0 = p
    k0 = N-k
    flipped = .TRUE.
  ELSE
    p0 = p
    q0 = 1-p
    k0 = k
    flipped = .FALSE.
  END IF
  mu = N*p0
  
  ! Will find approximately smaller sum first.
  ! Special cases.
  IF (mu .EQ. 0d0) THEN
    IF (k0 .EQ. 0) THEN
      ltsum = 0d0
      gtsum = 0d0
      pdf   = 1d0
    ELSE
      ltsum = 1d0
      gtsum = 0d0
      pdf   = 0d0
    END IF
  ELSE IF (k0 .EQ. 0) THEN
    ltsum = 0d0
    pdf   = q0**N
    gtsum = 1d0 - pdf
  ELSE IF (k0 .EQ. N) THEN
    gtsum = 0d0
    pdf   = p0**N
    ltsum = 1d0 - pdf
  ! Work our way down from k0.
  ELSE IF (k0 .LE. mu) THEN
    I  = k0
    xI = EXP(LOG_BINOMIAL_COEFF(N,I) + I*LOG(p0) + (N-I)*LOG(q0))
    pdf   = xI
    ltsum = 0d0
    DO WHILE ((I .GT. 0) .AND. (xI .GT. EPSILON(1d0)*ltsum))
      I     = I-1
      xI    = xI * ((I+1)*q0) / ((N-I)*p0)
      ltsum = ltsum + xI
    END DO
    gtsum = 1d0 - pdf - ltsum
  ! Work our way up from k0.
  ELSE
    I  = k0
    xI = EXP(LOG_BINOMIAL_COEFF(N,I) + I*LOG(p0) + (N-I)*LOG(q0))
    pdf   = xI
    gtsum = 0d0
    DO WHILE ((I .LT. N) .AND. (xI .GT. EPSILON(1d0)*gtsum))
      I     = I+1
      xI    = xI * ((N-I+1)*p0) / (I*q0)
      gtsum = gtsum + xI
    END DO
    ltsum = 1d0 - pdf - gtsum
  END IF
  ! reverse sums if we reversed p & q, k & N-k
  IF (flipped) THEN
    tmp   = ltsum
    ltsum = gtsum
    gtsum = tmp
  END IF
END SUBROUTINE BinomialSums

!----------------------------------------
! Random number for binomial distribution
!   p           Probability of success for each trial
!   N           Number of trials
FUNCTION RandBinomial(p,N) RESULT(K)
  IMPLICIT NONE
  INTEGER :: K
  INTEGER, INTENT(IN) :: N
  REAL*8, INTENT(IN) :: p
  REAL*8, PARAMETER :: MUMAX = 20d0
  LOGICAL :: flipped
  REAL*8 :: p0,mu
  ! Bad/special cases
  IF ((N .LT. 0) .OR. (p .LT. 0d0) .OR. (p .GT. 1d0)) THEN
    K = -HUGE(K)
    RETURN
  ELSE IF (N .EQ. 0) THEN
    K = 0
    RETURN
  END IF
  ! Can write P(k|p,N) = P(N-k|1-p,N) = P(k0,p0,N)
  ! Choose smaller of p or 1-p to work with
  IF (p .GT. 0.5d0) THEN
    p0 = 1-p
    flipped = .TRUE.
  ELSE
    p0 = p
    flipped = .FALSE.
  END IF
  mu = N*p0
  ! Special case
  IF (mu .EQ. 0d0) THEN
    K = 0
  ! Small mu case: iterative calculation faster
  ELSE IF (mu .LE. MUMAX) THEN
    K = IteratedK(p0,N)
  ! Otherwise, use ~ constant time rejection algorithm
  ELSE
    K = RejectionK(p0,N)
  END IF
  IF (flipped) K = N - K
  
  CONTAINS
  !-----------------------------
  ! Use explicit iteration to find random K.
  ! Assumes 0 < p <= 0.5 and N > 0.
  FUNCTION IteratedK(p,N) RESULT(K)
    IMPLICIT NONE
    INTEGER :: K
    INTEGER, INTENT(IN) :: N
    REAL*8, INTENT(IN) :: p
    REAL*8 :: q,mu,xK,y,sum
    q  = 1-p
    mu = N*p
    y  = rand()
    ! Work our way up from K=0
    IF (y .LE. 0.9d0) THEN
      K = NINT(mu - 9*SQRT(mu*q))
      IF (K .LE. 0) THEN
        K  = 0
        xK = q**N
      ELSE
        xK = EXP(LOG_BINOMIAL_COEFF(N,K) + K*LOG(p) + (N-K)*LOG(q))
      END IF
      sum = xK
      DO WHILE (sum .LT. y)
        K   = K+1
        xK  = xK * ((N-K+1)*p) / (K*q)
        sum = sum + xK
      END DO
    ! Work backwards from tail
    ELSE
      y = 1-y
      K = NINT(mu + 9*SQRT(mu*q))
      IF (K .GE. N) THEN
        K  = N
        xK = p**N
      ELSE
        xK = EXP(LOG_BINOMIAL_COEFF(N,K) + K*LOG(p) + (N-K)*LOG(q))
      END IF
      sum = xk
      DO WHILE (sum .LT. y)
        K   = K-1
        xK  = xK * ((K+1)*q) / ((N-K)*p)
        sum = sum + xK
      END DO
    END IF
  END FUNCTION IteratedK
  
  !-----------------------------
  ! Use rejection technique to find random K, as
  ! described in Hormann, JSCS 46, 101 (1993).
  ! Assumes 0 < p <= 0.5 and N > 0.
  FUNCTION RejectionK(p,N) RESULT(K)
    IMPLICIT NONE
    INTEGER :: K
    INTEGER, INTENT(IN) :: N
    REAL*8, INTENT(IN) :: p
    INTEGER :: M
    REAL*8 :: q,mu,sd,r,a,b,c,alpha,ur,vr,urvr,u,v,lnf,lnf_accept
    q  = 1-p
    mu = N*p
    ! Rejection technique: binomial density is proportional to
    !   f(x) = Pb(floor(x)) / Pb(M) <= 1
    ! where Pb is the binomial probability and  M = floor((N+1)*p0)
    ! is the mode. Use dominating distribution with inverse
    !   G(u) = (2alpha/(0.5-|u|) + b)*u + c
    ! with u in [-0.5,0.5].  Appropriate choices of constants
    ! determined by Hormann, JSCS 46, 101 (1993).
    sd = SQRT(mu*q)
    M = INT((N+1)*p)
    r = p/q
    c = mu + 0.5d0
    b = 1.15d0 + 2.53d0*sd
    a = -0.0873d0 + 0.0248d0*b + 0.01d0*p
    alpha = (2.83d0 + 5.1d0/b) * sd
    !ur = 0.43d0
    vr = 0.92d0 - 4.2d0/b
    urvr = 0.86d0*vr
    ! Will generate point in u in [-0.5,0.5], v in [0,1].
    ! Use a little magic to reduce number of uniform random
    ! numbers needed (on average).
    DO WHILE (.TRUE.)
      v = rand()
      ! Box: u in [-ur,ur], v in [0,vr].
      ! Always accepted.  Uniform in [0,ur*vr] converted to
      ! uniform in [-ur,ur], avoiding second random number.
      IF (v .LE. urvr) THEN
        u = v/vr - 0.43d0
        K = FLOOR((2*a/(0.5d0-ABS(u))+b)*u+c)
        EXIT
      END IF
      ! Outside of above box.  Two regions to consider:
      ! 1) Point in v > vr.  Need random u (can keep v).
      IF (v .GE. vr) THEN
        u = rand() - 0.5d0
      ! 2) Point in v < vr, u in [-0.5,-ur] or [ur,0.5].
      ! Uniform in [2ur*vr,vr] converted to uniform in
      ! [-0.5,-ur]+[ur,0.5].  Need random v in [0,vr].
      ELSE
        u = v/vr - 0.93d0
        u = SIGN(0.5d0,u) - u
        v = vr*rand()
      END IF
      ! Transform U to K = floor(G(U)).
      K = FLOOR((2*a/(0.5d0-ABS(u))+b)*u+c)
      IF ((K .LT. 0) .OR. (K .GT. N)) CYCLE
      ! Acceptance condition: V < f(G(U)) G'(U) / alpha
      ! where f(G(U)) = f(K).  Rewrite as:
      !   alpha*V/G'(U) < f(G(U)) = f(K)
      lnf_accept = LOG(alpha*v / (a/(0.5d0-ABS(u))**2+b))
      ! Calculate f(K).
      ! This step not terribly well optimized: could use recursion
      ! for K ~ M, though LOG_GAMMA is fairly fast as an intrinsic.
      IF (N .LE. 30) THEN
        lnf = LOG_BINOMIAL_COEFF(N,K) - LOG_BINOMIAL_COEFF(N,M)         &
              + (K-M)*LOG(r)
      ELSE
        lnf = LOG_GAMMA(M+1d0) + LOG_GAMMA((N-M)+1d0)                   &
              - LOG_GAMMA(K+1d0) - LOG_GAMMA((N-K)+1d0)                 &
              + (K-M)*LOG(r)
      END IF
      IF (lnf .GE. lnf_accept) EXIT
    END DO
  END FUNCTION RejectionK
  
END FUNCTION RandBinomial


END MODULE
