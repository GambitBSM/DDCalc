!#######################################################################
! DIRECT DETECTION RATES AND LIKELIHOODS
! Program to calculate the dark matter direct detection event rates
! and corresponding likelihoods/exclusion levels.  See DDCalc_run.help
! for instructions.
! 
! To see usage:
!     ./DDCalc_run --help
! 
!   Created by Chris Savage
!   University of Utah   (2013 - 2014)
!   Nordita              (2014 - 2015)
! 
!#######################################################################

PROGRAM DDCalc_run
  USE DDCALC
  IMPLICIT NONE
  CALL DDCalc_Main()
END PROGRAM


!#######################################################################
! END OF FILE
!#######################################################################

 
