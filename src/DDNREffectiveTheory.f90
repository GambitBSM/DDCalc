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


PUBLIC :: NRET_SFunctions, NRET_CreateCoeffList, NRET_SetDMSpin, NRET_GetParamIndex,&
   NRET_SetNRCoefficient, NRET_UpdateNRCoefficients

INTERFACE NRET_SFunctions
  MODULE PROCEDURE NRET_SFunctions_fct
END INTERFACE

INTERFACE NRET_CreateCoeffList
  MODULE PROCEDURE NRET_CreateCoeffList_fct
END INTERFACE

INTERFACE NRET_SetDMSpin
  MODULE PROCEDURE NRET_SetDMSpin_fct
END INTERFACE

INTERFACE NRET_GetParamIndex
  MODULE PROCEDURE NRET_GetParamIndex_fct
END INTERFACE

INTERFACE NRET_SetNRCoefficient
  MODULE PROCEDURE NRET_SetNRCoefficient_fct
END INTERFACE

INTERFACE NRET_UpdateNRCoefficients
  MODULE PROCEDURE NRET_UpdateNRCoefficients_fct
END INTERFACE

CONTAINS

! this function creates a list for the coefficients of the non-relativistic operators.
! All coefficients are set to zero, and the DM spin is set to -1.
FUNCTION NRET_CreateCoeffList_fct() &
     RESULT (par)
  IMPLICIT NONE
  REAL*8 :: par(45)
  par(:) = 0.0d0
  par(1) = -1.0d0
END FUNCTION


! this function sets the spin of the DM particle
SUBROUTINE NRET_SetDMSpin_fct(par, DMSpin)
  IMPLICIT NONE
  REAL*8, INTENT(INOUT) :: par(45)
  REAL*8, INTENT(IN) :: DMSpin
  par(1) = DMSpin
END SUBROUTINE


! This function finds the index of the parameter array corresponding to a given operator index and isospin index
SUBROUTINE NRET_GetParamIndex_fct(OpIndex, tau, i)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: OpIndex, tau
  INTEGER, INTENT(OUT) :: i

  IF (OpIndex.EQ.1) THEN
    i = 2
  ELSE IF (OpIndex.EQ.(-1)) THEN
    i = 4
  ELSE IF (OpIndex.EQ.3) THEN
    i = 6
  ELSE IF (OpIndex.EQ.4) THEN
    i = 8
  ELSE IF (OpIndex.EQ.(-4)) THEN
    i = 10
  ELSE IF (OpIndex.EQ.5) THEN
    i = 12
  ELSE IF (OpIndex.EQ.6) THEN
    i = 14
  ELSE IF (OpIndex.EQ.7) THEN
    i = 16
  ELSE IF (OpIndex.EQ.8) THEN
    i = 18
  ELSE IF (OpIndex.EQ.9) THEN
    i = 20
  ELSE IF (OpIndex.EQ.10) THEN
    i = 22
  ELSE IF (OpIndex.EQ.11) THEN
    i = 24
  ELSE IF (OpIndex.EQ.12) THEN
    i = 26
  ELSE IF (OpIndex.EQ.13) THEN
    i = 28
  ELSE IF (OpIndex.EQ.14) THEN
    i = 30
  ELSE IF (OpIndex.EQ.15) THEN
    i = 32
  ELSE IF (OpIndex.EQ.17) THEN
    i = 34
  ELSE IF (OpIndex.EQ.18) THEN
    i = 36
  ELSE
    WRITE (*,*) 'Error in NRET_SetNRCoefficient_fct: invalid operator index.'
    STOP
  END IF

  IF ((tau.NE.0) .AND. (tau.NE.1)) THEN
    WRITE (*,*) 'Error in NRET_SetNRCoefficient_fct: invalid value of tau.'
    STOP
  END IF

  i=i+tau

END SUBROUTINE


! this function sets a single coefficient to a non-zero value (in units GeV^(-2))
! OpIndex is the index of the non-relativistic operator, e.g. 6 for O_6.
! Additionally, OpIndex = -1 corresponds to q^2*O_1, and OpIndex = -4 corresponds to q^2*O_4
SUBROUTINE NRET_SetNRCoefficient_fct(par, OpIndex, tau, value)
  IMPLICIT NONE
  REAL*8, INTENT(INOUT) :: par(45)
  REAL*8, INTENT(IN) :: value
  INTEGER, INTENT(IN) :: OpIndex, tau
  INTEGER :: i

  CALL NRET_GetParamIndex(OpIndex, tau, i)
  par(i) = value
  CALL NRET_UpdateNRCoefficients(par)

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
  REAL*8, INTENT(IN) :: p(45)
  INTEGER, INTENT(IN) :: alpha, KE, Kiso
  REAL*8 :: a, b, mp, muT, ttt

  mp = PROTON_MASS
  muT = m*D%Miso(Kiso)/(m+D%Miso(Kiso)) ! DM-nucleus reduced mass
  a = (2*D%Miso(Kiso)*D%E(KE)*1d-6)/(mp**2) ! q^2/mp^2
  b = (2*D%Miso(Kiso)*D%E(KE)*1d-6)/(4.0d0*muT**2) ! q^2/(4*muT^2)

  S1S2 = (/ 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0 /)

  IF ( alpha.EQ.1 ) THEN
    S1S2(1)=(p(2)+a*p(4))**2-(p(1)*(1+p(1))*(4*b*p(18)**2+a*(-4*p(24)**2+b*(4*p(12)**2+p(34)**2))))/12.
    S1S2(2)=(p(2)+a*p(4))*(p(3)+a*p(5))-(p(1)*(1+p(1))*(4*b*p(18)*p(19)+a*(4*b*p(12)*p(13)-4*p(24)*p(25)+b*p(34)*p(35))))/12.
    S1S2(3)=(p(2)+a*p(4))*(p(3)+a*p(5))-(p(1)*(1+p(1))*(4*b*p(18)*p(19)+a*(4*b*p(12)*p(13)-4*p(24)*p(25)+b*p(34)*p(35))))/12.
    S1S2(4)=(p(3)+a*p(5))**2-(p(1)*(1+p(1))*(4*b*p(19)**2+a*(-4*p(25)**2+b*(4*p(13)**2+p(35)**2))))/12.
    S1S2(5)=(p(1)*(1+p(1))*(4*p(18)**2+a*(4*p(12)**2+p(34)**2)))/12.
    S1S2(6)=(p(1)*(1+p(1))*(4*a*p(12)*p(13)+4*p(18)*p(19)+a*p(34)*p(35)))/12.
    S1S2(7)=(p(1)*(1+p(1))*(4*a*p(12)*p(13)+4*p(18)*p(19)+a*p(34)*p(35)))/12.
    S1S2(8)=(p(1)*(1+p(1))*(4*p(19)**2+a*(4*p(13)**2+p(35)**2)))/12. 

  ELSE IF ( alpha.EQ.2 ) THEN
    S1S2(1)=(a*(3*a*p(6)**2+p(1)*(1+p(1))*(p(26)-a*p(32))**2))/12.
    S1S2(2)=(a*(3*a*p(6)*p(7)+p(1)*(1+p(1))*(p(26)-a*p(32))*(p(27)-a*p(33))))/12.
    S1S2(3)=(a*(3*a*p(6)*p(7)+p(1)*(1+p(1))*(p(26)-a*p(32))*(p(27)-a*p(33))))/12.
    S1S2(4)=(a*(3*a*p(7)**2+p(1)*(1+p(1))*(p(27)-a*p(33))**2))/12.
  
  ELSE IF ( alpha.EQ.3 ) THEN
    S1S2(1)=a*((p(2)+a*p(4))*p(6)+(p(1)*(1+p(1))*p(24)*(p(26)-a*p(32)))/3.)
    S1S2(2)=a*((p(3)+a*p(5))*p(6)+(p(1)*(1+p(1))*p(25)*(p(26)-a*p(32)))/3.)
    S1S2(3)=a*((p(2)+a*p(4))*p(7)+(p(1)*(1+p(1))*p(24)*(p(27)-a*p(33)))/3.)
    S1S2(4)=a*((p(3)+a*p(5))*p(7)+(p(1)*(1+p(1))*p(25)*(p(27)-a*p(33)))/3.)
  
  ELSE IF ( alpha.EQ.4 ) THEN
    S1S2(1)=(a*p(1)*(1+p(1))*(p(26)**2+a*p(28)**2))/12.
    S1S2(2)=(a*p(1)*(1+p(1))*(p(26)*p(27)+a*p(28)*p(29)))/12.
    S1S2(3)=(a*p(1)*(1+p(1))*(p(26)*p(27)+a*p(28)*p(29)))/12.
    S1S2(4)=(a*p(1)*(1+p(1))*(p(27)**2+a*p(29)**2))/12.

  ELSE IF ( alpha.EQ.5 ) THEN

    S1S2(1)=(p(1)*((1+p(1))*p(8)**2+2*a*(1+p(1))*p(8)*(p(10)+p(14))+a**2*p(10)*((1+p(1))*p(10)+2*p(14))))/12.
    S1S2(1)=S1S2(1)+(a**2*p(1)*p(14)*(p(14)+p(1)*(2*p(10)+p(14)))+3*a*p(22)**2-b*p(1)*p(26)**2)/12.
    S1S2(1)=S1S2(1)+(p(1)*(-2*b*(a*p(28)**2+p(1)*(p(26)**2+a*p(28)**2))+a*(1+p(1))*p(36)**2))/24.

    S1S2(2)=(p(1)*(1+p(1))*(p(8)*(p(9)+a*p(11))+a*(a*p(10)*p(11)+p(9)*(p(10)+p(14)))))/12.
    S1S2(2)=S1S2(2)+(a*p(1)*((1+p(1))*p(8)*p(15)+a*((1+p(1))*p(11)*p(14)+p(10)*p(15))))/12.
    S1S2(2)=S1S2(2)+(a**2*p(1)*(p(14)+p(1)*(p(10)+p(14)))*p(15)+3*a*p(22)*p(23)-b*p(1)*p(26)*p(27))/12.
    S1S2(2)=S1S2(2)+(p(1)*(-2*b*(a*p(28)*p(29)+p(1)*(p(26)*p(27)+a*p(28)*p(29)))+a*(1+p(1))*p(36)*p(37)))/24.

    S1S2(3)=S1S2(2)

    S1S2(4)=(p(1)*((1+p(1))*p(9)**2+2*a*(1+p(1))*p(9)*(p(11)+p(15))+a**2*p(11)*((1+p(1))*p(11)+2*p(15))))/12.
    S1S2(4)=S1S2(4)+(a**2*p(1)*p(15)*(p(15)+p(1)*(2*p(11)+p(15)))+3*a*p(23)**2-b*p(1)*p(27)**2)/12.
    S1S2(4)=S1S2(4)+(p(1)*(-2*b*(a*p(29)**2+p(1)*(p(27)**2+a*p(29)**2))+a*(1+p(1))*p(37)**2))/24.

    S1S2(5)=(p(1)*(1+p(1))*(p(26)**2+a*p(28)**2))/12.
    S1S2(6)=(p(1)*(1+p(1))*(p(26)*p(27)+a*p(28)*p(29)))/12.
    S1S2(7)=(p(1)*(1+p(1))*(p(26)*p(27)+a*p(28)*p(29)))/12.
    S1S2(8)=(p(1)*(1+p(1))*(p(27)**2+a*p(29)**2))/12.
  
  ELSE IF ( alpha.EQ.6 ) THEN

    S1S2(1)=(2*p(1)*(1+p(1))*p(8)**2+2*a**2*p(1)*(1+p(1))*p(10)**2+a*(-3*b*p(6)**2+4*p(1)*(1+p(1))*p(8)*p(10)))/24.
    S1S2(1)=S1S2(1)+(2*a*p(1)*(1+p(1))*p(20)**2-b*(3*p(16)**2+p(1)*(1+p(1))*(p(26)**2+a*p(30)**2)))/24.
    S1S2(1)=S1S2(1)+(a*p(1)*(1+p(1))*(-2*b*p(32)*(-2*p(26)+a*p(32))+p(36)**2))/48.
   
    S1S2(2)=(2*p(1)*(1+p(1))*p(8)*p(9)+a*(-3*b*p(6)*p(7)+2*p(1)*(1+p(1))*p(9)*p(10)))/24.
    S1S2(2)=S1S2(2)+(2*a*p(1)*(1+p(1))*p(8)*p(11)+2*a**2*p(1)*(1+p(1))*p(10)*p(11)-3*b*p(16)*p(17))/24.
    S1S2(2)=S1S2(2)+(-(p(1)*(b*(1+p(1))*p(26)*p(27)+a*(-2*(1+p(1))*p(20)*p(21)+b*p(30)*p(31))))/24.)
    S1S2(2)=S1S2(2)+(a*b*p(1)*(p(27)*p(32)+p(1)*(-(p(30)*p(31))+p(27)*p(32))+p(26)*p(33)))/24.
    S1S2(2)=S1S2(2)+(a*p(1)*(-2*b*(a*p(32)+p(1)*(-p(26)+a*p(32)))*p(33)+(1+p(1))*p(36)*p(37)))/48.

    S1S2(3)=S1S2(2)

    S1S2(4)=(2*p(1)*(1+p(1))*p(9)**2+2*a**2*p(1)*(1+p(1))*p(11)**2+a*(-3*b*p(7)**2+4*p(1)*(1+p(1))*p(9)*p(11)))/24.
    S1S2(4)=S1S2(4)+(2*a*p(1)*(1+p(1))*p(21)**2-b*(3*p(17)**2+p(1)*(1+p(1))*(p(27)**2+a*p(31)**2)))/24.
    S1S2(4)=S1S2(4)+(a*p(1)*(1+p(1))*(-2*b*p(33)*(-2*p(27)+a*p(33))+p(37)**2))/48.

    S1S2(5)=(3*p(16)**2+p(1)*(1+p(1))*p(26)**2+a*(3*p(6)**2+p(1)*(1+p(1))*p(30)**2))/24.
    S1S2(5)=S1S2(5)+(a*p(1)*(1+p(1))*p(32)*(-2*p(26)+a*p(32)))/24.

    S1S2(6)=(3*p(16)*p(17)+p(1)*(1+p(1))*p(26)*p(27)+a*(3*p(6)*p(7)+p(1)*(1+p(1))*(p(30)*p(31)-p(27)*p(32))))/24.
    S1S2(6)=S1S2(6)+(a*p(1)*(1+p(1))*(-p(26)+a*p(32))*p(33))/24.

    S1S2(7)=S1S2(6)

    S1S2(8)=(3*p(17)**2+p(1)*(1+p(1))*p(27)**2+a*(3*p(7)**2+p(1)*(1+p(1))*p(31)**2))/24.
    S1S2(8)=S1S2(8)+(a*p(1)*(1+p(1))*p(33)*(-2*p(27)+a*p(33)))/24.

  ELSE IF ( alpha.EQ.7 ) THEN
    S1S2(1)=(a*p(1)*(1+p(1))*(p(18)**2+a*(p(12)**2+p(34)**2/4.)))/3.
    S1S2(2)=(a*p(1)*(1+p(1))*(4*a*p(12)*p(13)+4*p(18)*p(19)+a*p(34)*p(35)))/12.
    S1S2(3)=(a*p(1)*(1+p(1))*(4*a*p(12)*p(13)+4*p(18)*p(19)+a*p(34)*p(35)))/12.
    S1S2(4)=(a*p(1)*(1+p(1))*(p(19)**2+a*(p(13)**2+p(35)**2/4.)))/3.  

  ELSE IF ( alpha.EQ.8 ) THEN
    S1S2(1)=(a*p(1)*(1+p(1))*((p(8)+a*p(10))*p(12)-p(18)*p(20)))/3.
    S1S2(2)=(a*p(1)*(1+p(1))*((p(9)+a*p(11))*p(12)-p(18)*p(21)))/3.
    S1S2(3)=(a*p(1)*(1+p(1))*((p(8)+a*p(10))*p(13)-p(19)*p(20)))/3.
    S1S2(4)=(a*p(1)*(1+p(1))*((p(9)+a*p(11))*p(13)-p(19)*p(21)))/3.

  END IF  

  !WRITE (*,*) 'a = ', a, ', b = ', b, ', S1S2 = ', S1S2


END FUNCTION


SUBROUTINE NRET_UpdateNRCoefficients_fct(par)
  IMPLICIT NONE
  REAL*8, INTENT(INOUT) :: par(45)
  INTEGER :: alpha, i
  INTEGER :: ptoalphaMatrix(288) =                                                                                     &
          (/1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0,&
            0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0,&
            1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0,&
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0,&
            0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1,&
            0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 1, 1, 0, 0, 1, 1,&
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0,&
            0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0/)

  DO alpha = 1,8
    par(37 + alpha) = 0
    DO i = 1,36
      par(37 + alpha) = par(37 + alpha) + ptoalphaMatrix(i+(alpha-1)*36) * par(i+1)**2
    END DO
  END DO

END SUBROUTINE

END MODULE
