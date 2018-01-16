MODULE DDNREffectiveTheory

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! DDNREffectiveTheory
!    Routines to calculate reponse functions for the non-relativistic
!    effective theory of dark matter-nucleon interactions
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX


USE DDConstants
USE DDHalo
USE DDCouplings
USE DDTypes

IMPLICIT NONE
PRIVATE


PUBLIC :: NRET_SFunctions, NRET_CreateCoeffList, NRET_SetDMSpin, NRET_SetNRCoefficient

INTERFACE NRET_SFunctions
  MODULE PROCEDURE NRET_SFunctions_fct
END INTERFACE

INTERFACE NRET_CreateCoeffList
  MODULE PROCEDURE NRET_CreateCoeffList_fct
END INTERFACE

INTERFACE NRET_SetDMSpin
  MODULE PROCEDURE NRET_SetDMSpin_fct
END INTERFACE

INTERFACE NRET_SetNRCoefficient
  MODULE PROCEDURE NRET_SetNRCoefficient_fct
END INTERFACE


CONTAINS

! this function creates a list for the coefficients of the non-relativistic operators.
! All coefficients are set to zero, and the DM spin is set to -1.
FUNCTION NRET_CreateCoeffList_fct() &
     RESULT (par)
  IMPLICIT NONE
  REAL*8 :: par(37)
  par(:) = 0.0d0
  par(1) = -1.0d0
END FUNCTION


! this function sets the spin of the DM particle
SUBROUTINE NRET_SetDMSpin_fct(par, DMSpin)
  IMPLICIT NONE
  REAL*8, INTENT(INOUT) :: par(37)
  REAL*8, INTENT(IN) :: DMSpin
  par(1) = DMSpin
END SUBROUTINE


! this function sets a single coefficient to a non-zero value (in units GeV^(-2))
SUBROUTINE NRET_SetNRCoefficient_fct(par, OpString, tau, value)
  IMPLICIT NONE
  REAL*8, INTENT(INOUT) :: par(37)
  CHARACTER(LEN=*), INTENT(IN) :: OpString
  REAL*8, INTENT(IN) :: value
  INTEGER, INTENT(IN) :: tau

  IF ((OpString.EQ.'Op1') .AND. (tau.EQ.0)) THEN
    par(2)=value
  !ELSE IF ((OpString.EQ.'Op1') .AND. (tau.EQ.0)) THEN
  !  par(...
  END IF
END SUBROUTINE

! Auxiliary function providing the response functions S1 and S2 for the non-relativistic effective theory.
! alpha (running from 1 to 8) denotes the index of the corresponding nuclear response function.
! KE and Kiso are the indices of the energy and isotope array.
! The return value is a list [S1_00, S1_01, S1_10, S1_11, S2_00, S2_01, S2_10, S2_11] in units GeV^(-4).
! The notation follows appendix A of 1607.04418 . 
FUNCTION NRET_SFunctions_fct(D, m, p, alpha, KE, Kiso)   &
    RESULT (S1S2)
  IMPLICIT NONE

  REAL*8 :: S1S2(8)
  TYPE(DetectorStruct), INTENT(IN) :: D
  REAL*8, INTENT(IN) :: m
  REAL*8, INTENT(IN) :: p(37)
  INTEGER, INTENT(IN) :: alpha, KE, Kiso
  REAL*8 :: a, b, mp, muT, ttt

  mp = PROTON_MASS
  muT = m*D%Miso(Kiso)/(m+D%Miso(Kiso)) ! DM-nucleus reduced mass
  a = (2*D%Miso(Kiso)*D%E(KE)*1d-6)/(mp**2) ! q^2/mp^2
  b = (2*D%Miso(Kiso)*D%E(KE)*1d-6)/(4.0d0*muT**2) ! q^2/(4*muT^2)

  S1S2 = (/ 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0 /)

  IF ( alpha.EQ.1 ) THEN
    S1S2(1)=(p(2)+a*p(4))**2-(p(1)*(1+p(1))*(a*b*p(12)**2+b*p(18)**2-a*p(24)**2))/3.
    S1S2(2)=(p(2)+a*p(4))*(p(3)+a*p(5))-(p(1)*(1+p(1))*(a*b*p(12)*p(13)+b*p(18)*p(19)-a*p(24)*p(25)))/3.
    S1S2(3)=(p(2)+a*p(4))*(p(3)+a*p(5))-(p(1)*(1+p(1))*(a*b*p(12)*p(13)+b*p(18)*p(19)-a*p(24)*p(25)))/3.
    S1S2(4)=(p(3)+a*p(5))**2-(p(1)*(1+p(1))*(a*b*p(13)**2+b*p(19)**2-a*p(25)**2))/3.
    S1S2(5)=(p(1)*(1+p(1))*(a*p(12)**2+p(18)**2))/3.
    S1S2(6)=(p(1)*(1+p(1))*(a*p(12)*p(13)+p(18)*p(19)))/3.
    S1S2(7)=(p(1)*(1+p(1))*(a*p(12)*p(13)+p(18)*p(19)))/3.
    S1S2(8)=(p(1)*(1+p(1))*(a*p(13)**2+p(19)**2))/3.  

  ELSE IF ( alpha.EQ.2 ) THEN
    S1S2(1)=(a*(3*a*p(6)**2+p(1)*(1+p(1))*(p(26)-a*p(32))**2))/12.
    S1S2(2)=(a*(3*a*p(6)*p(7)+p(1)*(1+p(1))*(p(26)-a*p(32))*(p(27)-a*p(33))))/12.
    S1S2(3)=(a*(3*a*p(6)*p(7)+p(1)*(1+p(1))*(p(26)-a*p(32))*(p(27)-a*p(33))))/12.
    S1S2(4)=(a*(3*a*p(7)**2+p(1)*(1+p(1))*(p(27)-a*p(33))**2))/12.
  
  ELSE IF ( alpha.EQ.3 ) THEN
  
  ELSE IF ( alpha.EQ.4 ) THEN

  ELSE IF ( alpha.EQ.5 ) THEN
  
  ELSE IF ( alpha.EQ.6 ) THEN
  
  ELSE IF ( alpha.EQ.7 ) THEN
  
  ELSE IF ( alpha.EQ.8 ) THEN

  END IF  


END FUNCTION



END MODULE
