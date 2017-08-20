MODULE DDOutput

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! DDOutput
!    DDCalc routines for export of data.
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

USE DDConstants
USE DDTypes
USE DDInput
USE DDStats
USE DDUtils
USE DDRates
USE DDCouplings
USE DDCommandLine

IMPLICIT NONE
PRIVATE

PUBLIC :: WriteCommandHeader, WriteHaloHeader, WriteDetectorHeader, &
          WriteLimitsSDData, WriteLimitsSDHeader, WriteLimitsSDColumnHeader, &
          WriteLimitsSIData, WriteLimitsSIHeader, WriteLimitsSIColumnHeader, &
          WriteConstraintsSDData, WriteConstraintsSDHeader, WriteConstraintsSDColumnHeader, &
          WriteConstraintsSIData, WriteConstraintsSIHeader, WriteConstraintsSIColumnHeader, &
          WriteWIMPHeader, WriteEventsByMassHeader, WriteEventsByMassColumnHeader, &
          WriteEventsByMassData, WriteSpectrumHeader, WriteSpectrumData, &
          WriteSpectrumColumnHeader, WriteInteractiveHeader, WriteEventsAndLikelihoodsData, &
          WriteEventsAndLikelihoodsColumnHeader, WriteEventsAndLikelihoodsHeader, &
          WriteLogLikelihoodHeader, WriteLogPValueHeader
       
CONTAINS


!-----------------------------------------------------------------------
! Prints the given number of empty comment lines (just the comment
! prefix).  Utility function.
! 
SUBROUTINE WriteEmptyCommentLines(N)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: N
  INTEGER :: I
  DO I = 1,N
    WRITE(*,'(A)') COMMENT_PREFIX
  END DO
END SUBROUTINE


!-----------------------------------------------------------------------
! Write command header.
! Outputs the command used and date.
! 
! Optional input arguments:
!   extra_lines     Blank lines (albeit with prefix) to add after
!                   output
! 
SUBROUTINE WriteCommandHeader(extra_lines)
  IMPLICIT NONE
  INTEGER, INTENT(IN), OPTIONAL :: extra_lines
  CHARACTER*1024 :: cmd
  CHARACTER*10 :: date,time
  CHARACTER*100 :: datetime
  
  ! Command
  CALL GetFullCommand(cmd)
  WRITE(*,'(A)') COMMENT_PREFIX &
      //  'Command: ' // TRIM(cmd)
  
  ! Date & time
  CALL DATE_AND_TIME(date,time)
  datetime = 'Generated on '                                            &
      // date(1:4) // '-' // date(5:6) // '-' // date(7:8)  // ' at '   &
      // time(1:2) // ':' // time(3:4) // ':' // time(5:10)
  WRITE(*,'(A)') COMMENT_PREFIX &
      // TRIM(datetime)
  
  ! Version
  !WRITE(*,'(A)') COMMENT_PREFIX &
  !    // 'DDCalc version: ' // TRIM(VERSION_STRING)
  !WRITE(*,'(A)') COMMENT_PREFIX &
  !    // 'See/cite C. Savage et al., arxiv:15MM.NNNNN.'
  WRITE(*,'(A)') COMMENT_PREFIX &
      // 'DDCalc version ' // TRIM(VERSION_STRING) // '.  See/cite C. Savage et al., arxiv:15MM.NNNNN.'
  
  WRITE(*,'(A)') COMMENT_PREFIX
  
  IF (PRESENT(extra_lines)) CALL WriteEmptyCommentLines(extra_lines)
  
END SUBROUTINE


!-----------------------------------------------------------------------
! Prints to standard output information regarding the WIMP.
! 
! Optional input arguments:
!   extra_lines     Blank lines (albeit with prefix) to add after
!                   output
! 
SUBROUTINE WriteWIMPHeader(WIMP,extra_lines)
  IMPLICIT NONE
  TYPE(WIMPStruct), INTENT(IN) :: WIMP
  INTEGER, INTENT(IN), OPTIONAL :: extra_lines
  INTEGER :: paramID
  
  !WRITE(*,'(A)') COMMENT_LINE
  
  ! Mass
  WRITE(*,'(A,F11.4)') COMMENT_PREFIX &
      // 'WIMP mass [GeV]     =',WIMP%m
  WRITE(*,'(A)') COMMENT_PREFIX

  ! Type
  WRITE(*,'(A,A)') COMMENT_PREFIX &
      // 'WIMP type           =',WIMP%DMtype
  WRITE(*,'(A)') COMMENT_PREFIX
  
  ! Parameters
  DO paramID = 1,WIMP%Nparams

    WRITE(*,'(A,1(X,I2),A,1(2X,1PG12.4))') COMMENT_PREFIX &
        // ' WIMP parameter',paramID,':',WIMP%params(paramID)
    WRITE(*,'(A)') COMMENT_PREFIX

  END DO
  
  IF (PRESENT(extra_lines)) CALL WriteEmptyCommentLines(extra_lines)
  
END SUBROUTINE


!-----------------------------------------------------------------------
! Prints to standard output information regarding the dark matter halo,
! notably the density and velocity distribution.
! 
! Optional input arguments:
!   extra_lines     Blank lines (albeit with prefix) to add after
!                   output
! 
SUBROUTINE WriteHaloHeader(Halo, extra_lines)
  IMPLICIT NONE
  TYPE(HaloStruct), INTENT(IN) :: Halo
  INTEGER, INTENT(IN), OPTIONAL :: extra_lines
  
  !WRITE(*,'(A)') COMMENT_LINE
  
  ! Solar motion
  WRITE(*,'(A)') COMMENT_PREFIX &
      // 'Sun''s velocity in the Galactic rest frame [km/s] in Galactic'
  WRITE(*,'(A)') COMMENT_PREFIX &
      // 'coordinates (U,V,W), where U is anti-radial (towards the'
  WRITE(*,'(A)') COMMENT_PREFIX &
      // 'Galactic center), V is in the direction of the disk rotation,'
  WRITE(*,'(A)') COMMENT_PREFIX &
      // 'and W is vertical (out of the disk plane):'
  WRITE(*,'(A,3(1X,F8.2),3X,A)') COMMENT_PREFIX &
      // '  ',Halo%vsun
  WRITE(*,'(A)') COMMENT_PREFIX
  
  ! Density
  WRITE(*,'(A,F11.4)') COMMENT_PREFIX &
      // 'Local dark matter density [GeV/cm^3]  =',Halo%rho
  WRITE(*,'(A)') COMMENT_PREFIX
  
  ! Halo velocity distribution
  IF (.NOT. Halo%tabulated) THEN
    ! Maxwell-Boltzmann parameters
    WRITE(*,'(A)') COMMENT_PREFIX &
        // 'SHM-like velocity distribution (Maxwell-Boltzmann with finite cutoff):'
    WRITE(*,'(A,3(1X,F8.2),3X,A)') COMMENT_PREFIX &
        // '  Bulk motion (galactic frame) [km/s] =',Halo%vbulk
    WRITE(*,'(A,F9.2)') COMMENT_PREFIX &
        // '  Most probable speed (v0) [km/s]     =',Halo%v0
    IF (Halo%vesc .GE. 1e6) THEN
      WRITE(*,'(A)') COMMENT_PREFIX &
        // '  Escape speed (vesc) [km/s]          = (infinite)'
    ELSE
      WRITE(*,'(A,F9.2)') COMMENT_PREFIX &
        // '  Escape speed (vesc) [km/s]          =',Halo%vesc
    END IF
  ELSE IF (Halo%eta_file .NE. '') THEN
    ! Tabulated from file
    !WRITE(*,'(A)') COMMENT_PREFIX &
    !    // 'The mean inverse speed was provided in tabulated form in the following'
    !WRITE(*,'(A)') COMMENT_PREFIX &
    !    // 'file:'
    !WRITE(*,'(A)') COMMENT_PREFIX &
    !    // '    ' // TRIM(Halo%eta_file)
    WRITE(*,'(A,A)') COMMENT_PREFIX &
        // 'Mean inverse speed tabulation file    = ',TRIM(Halo%eta_file)
  ELSE
    ! User provided tabulation
    ! This really shouldn't occur here because there is no way
    ! to provide the tabulation this way in program modes!
    WRITE(*,'(A)') COMMENT_PREFIX &
        // 'The mean inverse speed was provided in tabulated form by the user.'
  END IF
  WRITE(*,'(A)') COMMENT_PREFIX
  
  IF (PRESENT(extra_lines)) CALL WriteEmptyCommentLines(extra_lines)
  
END SUBROUTINE


!-----------------------------------------------------------------------
! Prints to standard output information regarding the detector.
! 
! Optional input arguments:
!   extra_lines     Blank lines (albeit with prefix) to add after
!                   output
! 
SUBROUTINE WriteDetectorHeader(Detector,extra_lines)
  IMPLICIT NONE
  TYPE(DetectorStruct), INTENT(IN) :: Detector
  INTEGER, INTENT(IN), OPTIONAL :: extra_lines
  
  !WRITE(*,'(A)') COMMENT_LINE
  
  ! Exposure, observed events, expected background events
  WRITE(*,'(A,1PG12.4)') COMMENT_PREFIX &
      // 'Detector exposure [kg day]            =',Detector%exposure
  IF (Detector%Nevents(0) .GE. 0) THEN
    WRITE(*,'(A,I6)') COMMENT_PREFIX &
      // 'Observed events                       =',Detector%Nevents(0)
  END IF
  WRITE(*,'(A,F11.4)') COMMENT_PREFIX &
      // 'Average expected background events    =',Detector%Backgr(0)
  WRITE(*,'(A)') COMMENT_PREFIX
  
  ! Detector isotopes
  WRITE(*,'(A,I6)') COMMENT_PREFIX &
      // 'Isotopes                              =',Detector%Niso
  WRITE(*,'(A,99(3X,I4,2X))') COMMENT_PREFIX &
      // '  Atomic number Z       ',Detector%Ziso
  WRITE(*,'(A,99(3X,I4,2X))') COMMENT_PREFIX &
      // '  Atomic mass number A  ',Detector%Aiso
  WRITE(*,'(A,99(1X,F8.5))') COMMENT_PREFIX &
      // '  Mass fraction         ',Detector%fiso
  WRITE(*,'(A)') COMMENT_PREFIX
  
  ! Efficiencies and intervals/bins
  WRITE(*,'(A,A)') COMMENT_PREFIX &
      // 'Efficiency file                       = ',TRIM(Detector%eff_file)
  !IF (Detector%Nbins .GT. 0) THEN
  IF (Detector%intervals .AND. (Detector%Nbins .GT. 0)) THEN
    WRITE(*,'(A,I6)') COMMENT_PREFIX &
      // 'Number of bins/sub-intervals          =',Detector%Nbins
  END IF
  WRITE(*,'(A)') COMMENT_PREFIX
  
  IF (PRESENT(extra_lines)) CALL WriteEmptyCommentLines(extra_lines)
  
END SUBROUTINE


!-----------------------------------------------------------------------
! Write the interactive-mode instruction header.
! 
! Optional input arguments:
!   extra_lines     Blank lines (albeit with prefix) to add after
!                   output
! 
SUBROUTINE WriteInteractiveHeader(extra_lines)
  IMPLICIT NONE
  INTEGER, INTENT(IN), OPTIONAL :: extra_lines
  
  WRITE(*,'(A)') COMMENT_PREFIX &
      // 'Enter WIMP parameters in one of the following formats:'
  WRITE(*,'(A)') COMMENT_PREFIX &
      // '  $>  m'
  WRITE(*,'(A)') COMMENT_PREFIX &
      // '  $>  m sigmaSI'
  WRITE(*,'(A)') COMMENT_PREFIX &
      // '  $>  m sigmaSI sigmaSD'
  WRITE(*,'(A)') COMMENT_PREFIX &
      // '  $>  m sigmaSI sigmapSD sigmanSD'
  WRITE(*,'(A)') COMMENT_PREFIX &
      // '  $>  m sigmapSI sigmanSI sigmapSD sigmanSD'
  WRITE(*,'(A)') COMMENT_PREFIX &
      // 'where m is the WIMP mass [GeV] and sigma is a WIMP-nucleon cross-section'
  WRITE(*,'(A)') COMMENT_PREFIX &
      // '[pb].  In the first case, all cross-sections are set to 1 pb.  In the'
  WRITE(*,'(A)') COMMENT_PREFIX &
      // 'second case, SD couplings are set to zero.  Negative cross-sections may'
  WRITE(*,'(A)') COMMENT_PREFIX &
      // 'be given to indicate the corresponding coupling is negative.  A blank'
  WRITE(*,'(A)') COMMENT_PREFIX &
      // 'line will terminate the program (as will an invalid format).'
  !WRITE(*,'(A)') COMMENT_PREFIX &
  !    // '  $>  m GpSI GnSI GpSD GnSD'
  !WRITE(*,'(A)') COMMENT_PREFIX &
  !    // 'where m is the WIMP mass [GeV], sigma is a WIMP-nucleon cross-section'
  !WRITE(*,'(A)') COMMENT_PREFIX &
  !    // '[pb], and G is a WIMP-nucleon coupling [GeV^-2].  In the first case,'
  !WRITE(*,'(A)') COMMENT_PREFIX &
  !    // 'all cross-sections are set to 1 pb.  In the second case, SD couplings'
  !WRITE(*,'(A)') COMMENT_PREFIX &
  !    // 'are set to zero.  Negative cross-sections may be given to indicate'
  !WRITE(*,'(A)') COMMENT_PREFIX &
  !    // 'the corresponding coupling is negative.  A blank line will terminate'
  !WRITE(*,'(A)') COMMENT_PREFIX &
  !    // 'the program (as will an invalid format).'
  !WRITE(*,'(A)') COMMENT_PREFIX &
  !    // ''
  
  IF (PRESENT(extra_lines)) CALL WriteEmptyCommentLines(extra_lines)
  
END SUBROUTINE


!-----------------------------------------------------------------------
! Write log-likelihood header.
! 
! Optional input arguments:
!   extra_lines     Blank lines (albeit with prefix) to add after
!                   output
! 
SUBROUTINE WriteLogLikelihoodHeader(extra_lines)
  IMPLICIT NONE
  INTEGER, INTENT(IN), OPTIONAL :: extra_lines
  
  WRITE(*,'(A)') COMMENT_PREFIX &
      // 'The log-likelihood for the given parameters is given below,'
  WRITE(*,'(A)') COMMENT_PREFIX &
      // 'defined as L = P(N|s+b) where P is the Poisson distribution of'
  WRITE(*,'(A)') COMMENT_PREFIX &
      // 'N observed events given an expected signal s and background b.'
  
  IF (PRESENT(extra_lines)) CALL WriteEmptyCommentLines(extra_lines)
  
END SUBROUTINE


!-----------------------------------------------------------------------
! Write log p-value header.
! 
! Optional input arguments:
!   extra_lines     Blank lines (albeit with prefix) to add after
!                   output
! 
SUBROUTINE WriteLogPValueHeader(Detector,extra_lines)
  IMPLICIT NONE
  TYPE(DetectorStruct), INTENT(IN) :: Detector
  INTEGER, INTENT(IN), OPTIONAL :: extra_lines
  
  IF (Detector%intervals .AND. (Detector%Nbins .EQ. Detector%Nevents(0)+1)) THEN
    WRITE(*,'(A)') COMMENT_PREFIX &
        // 'The log of the p-value for the given parameters is given below,'
    WRITE(*,'(A)') COMMENT_PREFIX &
        // 'calculated using the maximum gap method.  See:'
    WRITE(*,'(A)') COMMENT_PREFIX &
        // '  S. Yellin, PRD 66, 032005 (2002) [physics/0203002]'
  ELSE
    WRITE(*,'(A)') COMMENT_PREFIX &
        // 'The log of the p-value for the given parameters are given below,'
    WRITE(*,'(A)') COMMENT_PREFIX &
        // 'calculated using a Poisson distribution on the number of events.'
  END IF
  
  IF (PRESENT(extra_lines)) CALL WriteEmptyCommentLines(extra_lines)
  
END SUBROUTINE


!-----------------------------------------------------------------------
! Prints to standard output a header containing information regarding
! the events and likelihoods data to follow.
! 
! Optional input arguments:
!   extra_lines     Blank lines (albeit with prefix) to add after
!                   output
! 
SUBROUTINE WriteEventsAndLikelihoodsHeader(Detector,extra_lines)
  IMPLICIT NONE
  TYPE(DetectorStruct), INTENT(IN) :: Detector
  INTEGER, INTENT(IN), OPTIONAL :: extra_lines
  
  !WRITE(*,'(A)') COMMENT_LINE
  
  IF (VerbosityLevel .GE. 2) THEN
    WRITE(*,'(A)') COMMENT_PREFIX &
        // 'The columns below contain the following data:'
  END IF
  
  IF (VerbosityLevel .GE. 3) THEN
    WRITE(*,'(A)') COMMENT_PREFIX &
        // '  m            WIMP mass [GeV].'
    WRITE(*,'(A)') COMMENT_PREFIX &
        // '  GpSI,GnSI,GpSD,GnSD'
    WRITE(*,'(A)') COMMENT_PREFIX &
        // '               WIMP-nucleon spin-independent and spin-dependent couplings [GeV^2].'
    WRITE(*,'(A)') COMMENT_PREFIX &
        // '  sigmapSI,sigmanSI,sigmapSD,sigmanSD'
    WRITE(*,'(A)') COMMENT_PREFIX &
        // '               WIMP-nucleon SI and SD scattering cross-sections [pb].'
  END IF
  
  IF (VerbosityLevel .GE. 2) THEN
    WRITE(*,'(A)') COMMENT_PREFIX &
        // '  observed     The observed number of events.'
    WRITE(*,'(A)') COMMENT_PREFIX &
        // '  background   Average expected background events.'
    WRITE(*,'(A)') COMMENT_PREFIX &
        // '  signal(SI)   Average expected spin-independent signal events.'
    WRITE(*,'(A)') COMMENT_PREFIX &
        // '  signal(SD)   Average expected spin-dependent signal events.'
    WRITE(*,'(A)') COMMENT_PREFIX &
        // '  log(L)       Log-likelihood using the Poisson distribution (signal+background).'
    IF (Detector%intervals .AND. (Detector%Nbins .EQ. Detector%Nevents(0)+1)) THEN
      WRITE(*,'(A)') COMMENT_PREFIX &
        // '  log(p)       Log of the p-value determined using the maximum gap method;'
      WRITE(*,'(A)') COMMENT_PREFIX &
        // '               see Yellin, Phys. Rev. D 66, 032005 (2002) [physics/0203002].'
    ELSE
      WRITE(*,'(A)') COMMENT_PREFIX &
        // '  log(p)       Log of the p-value determined using the Poisson distribution'
      WRITE(*,'(A)') COMMENT_PREFIX &
        // '               (signal only: no background subtraction).'
    END IF
  END IF
  
  IF (VerbosityLevel .GE. 2) THEN
    WRITE(*,'(A)') COMMENT_PREFIX
  END IF
  
  IF (PRESENT(extra_lines)) CALL WriteEmptyCommentLines(extra_lines)
  
END SUBROUTINE


!-----------------------------------------------------------------------
! Prints to standard output a column header for the events and
! likelihoods data to follow.
! 
SUBROUTINE WriteEventsAndLikelihoodsColumnHeader()
  IMPLICIT NONE
  
  SELECT CASE (VerbosityLevel)
  CASE (:0)
    CONTINUE
  CASE (1:2)
    WRITE(*,'(A,1(1X,A8),3(1X,A11),2(1X,A11))') COMMENT_PREFIX,         &
        'observed','background ',           &
        '  log(L)   ','  log(p)   '
  CASE (3:)
    WRITE(*,'(A,1(1X,A10),4(1X,A10),4(1X,A10),1(1X,A8),3(1X,A11),2(1X,A11))') &
        COMMENT_PREFIX,'    mass  ',                                    &
        '   GpSI   ','   GnSI   ','   GpSD   ','   GnSD   ',            &
        ' sigmapSI ',' sigmanSI ',' sigmapSD ',' sigmanSD ',            &
        'observed','background ','signal(SI) ','signal(SD) ',           &
        '  log(L)   ','  log(p)   '
  END SELECT
  
END SUBROUTINE


!-----------------------------------------------------------------------
! Prints to standard output the events and likelihoods data for the
! current WIMP mass and couplings.
! 
SUBROUTINE WriteEventsAndLikelihoodsData(Detector,WIMP)
  IMPLICIT NONE
  TYPE(DetectorStruct), INTENT(IN) :: Detector
  TYPE(WIMPStruct), INTENT(IN) :: WIMP
  REAL*8 :: lnLike
  REAL*8 :: lnp
  
  ! Get log-likelihood and p-value
  lnLike = DDCalc_LogLikelihood(Detector)
  lnp    = DDCalc_LogPValue(Detector)
  
  ! Columns to print depend on the verbosity level.
  ! For data, only absolute value of verbosity is used.
 
  IF (ABS(VerbosityLevel) .GE. 2) THEN
    WRITE(*,'(A,1(2X,I5,2X),3(1X,1PG11.4),2(1X,1PG11.4))')              &
        DATA_PREFIX,                                                    &
        Detector%Nevents(0),Detector%Backgr(0),                         &
        lnLike,lnp
  END IF
 
  
END SUBROUTINE


!-----------------------------------------------------------------------
! Prints to standard output a header containing information regarding
! the events and likelihoods data to follow.
! 
! Optional input arguments:
!   extra_lines     Blank lines (albeit with prefix) to add after
!                   output
! 
SUBROUTINE WriteSpectrumHeader(extra_lines)
  IMPLICIT NONE
  INTEGER, INTENT(IN), OPTIONAL :: extra_lines
  
  !WRITE(*,'(A)') COMMENT_LINE
  
  IF (VerbosityLevel .GE. 4) THEN
    WRITE(*,'(A)') COMMENT_PREFIX &
        // 'The table below gives differential rate components at reference'
    WRITE(*,'(A)') COMMENT_PREFIX &
        // 'couplings that yield WIMP-nucleon scattering cross-sections of 1 pb.'
    WRITE(*,'(A)') COMMENT_PREFIX &
        // 'Separate columns are given for the spectrum contribution arising from'
    WRITE(*,'(A)') COMMENT_PREFIX &
        // 'the proton-proton [dRdEpp0], neutron-neutron [dRdEnn0], and proton-'
    WRITE(*,'(A)') COMMENT_PREFIX &
        // 'neutron [dRdEpn0] components of the cross-section/coupling formula.'
    WRITE(*,'(A)') COMMENT_PREFIX &
        // 'This allows rates for arbitrary couplings to be constructed as follows:'
    WRITE(*,'(A)') COMMENT_PREFIX &
        // '    dRdE = (sigmap/[pb])*dRdEpp0 + (sigman/[pb])*dRdEnn0'
    WRITE(*,'(A)') COMMENT_PREFIX &
        // '              +/- (sqrt{sigmap*sigman}/[pb])*dRdEpn0'
    WRITE(*,'(A)') COMMENT_PREFIX &
        // 'where sigmap and sigman are the WIMP-nucleon scattering cross-sections'
    WRITE(*,'(A)') COMMENT_PREFIX &
        // 'in pb.  The sign of the cross-term should be the same as the sign of'
    WRITE(*,'(A)') COMMENT_PREFIX &
        // 'the product of couplings Gp*Gn.'
    !WRITE(*,'(A)') COMMENT_PREFIX &
    !    // ''
    WRITE(*,'(A)') COMMENT_PREFIX
  END IF
  
  IF (VerbosityLevel .GE. 2) THEN
    WRITE(*,'(A)') COMMENT_PREFIX &
        // 'The columns below contain the following data:'
  END IF
  
  IF (VerbosityLevel .GE. 2) THEN
    WRITE(*,'(A)') COMMENT_PREFIX &
        // '  E            Recoil energy [keV].'
  END IF
  
  SELECT CASE (VerbosityLevel)
  CASE (2)
    WRITE(*,'(A)') COMMENT_PREFIX &
        // '  dRdE         Differential recoil spectrum [cpd/kg/keV].'
  CASE (3)
    WRITE(*,'(A)') COMMENT_PREFIX &
        // '  dRdE         Differential recoil spectrum [cpd/kg/keV].'
    WRITE(*,'(A)') COMMENT_PREFIX &
        // '               One column for each of the total, SI, and SD spectra.'
  CASE (4:)
    WRITE(*,'(A)') COMMENT_PREFIX &
        // '  dRdEpp0,dRdEpn0,dRdEnn0'
    WRITE(*,'(A)') COMMENT_PREFIX &
        // '               Differential recoil spectrum components [cpd/kg/keV],'
    WRITE(*,'(A)') COMMENT_PREFIX &
        // '               given separately for SI and SD interactions.'
  END SELECT
  
  IF (VerbosityLevel .GE. 2) THEN
    WRITE(*,'(A)') COMMENT_PREFIX
  END IF
  
  IF (PRESENT(extra_lines)) CALL WriteEmptyCommentLines(extra_lines)
  
END SUBROUTINE


!-----------------------------------------------------------------------
! Prints to standard output a column header for the tabulated
! differential rate spectrum to follow.
! 
SUBROUTINE WriteSpectrumColumnHeader()
  IMPLICIT NONE
  
  ! Columns to print depend on the verbosity level.
  ! For data, only absolute value of verbosity is used.
  
  SELECT CASE (ABS(VerbosityLevel))
  CASE (1)
    WRITE(*,'(A,1(1X,A10))',ADVANCE='NO') COMMENT_PREFIX,'  E [keV] '
    WRITE(*,'(1(1X,A12))',ADVANCE='NO')                                 &
        ' dR/dE [dru]'
  CASE (2)
    ! Differential rate
    WRITE(*,'(A,1(1X,A10))',ADVANCE='NO') COMMENT_PREFIX,'    E     '
    WRITE(*,'(1(1X,A12))',ADVANCE='NO')                                 &
        '   dR/dE    '
  CASE (3)
    ! Combined, si, sd
    WRITE(*,'(A,1(1X,A10))',ADVANCE='NO') COMMENT_PREFIX,'    E     '
    WRITE(*,'(3(1X,A12))',ADVANCE='NO')                                 &
        '   dR/dE    ',' dR/dE(SI)  ',' dR/dE(SD)  '
  CASE (4:)
    ! Reference rate components
    WRITE(*,'(A,1(1X,A10))',ADVANCE='NO') COMMENT_PREFIX,'    E     '
    WRITE(*,'(6(1X,A12))',ADVANCE='NO')                                 &
        'dRdEpp0(SI) ','dRdEpn0(SI) ','dRdEnn0(SI) ',                   &
        'dRdEpp0(SD) ','dRdEpn0(SD) ','dRdEnn0(SD) '
  END SELECT
  IF (VerbosityLevel .GE. 1) WRITE(*,'(A)') ''
  
END SUBROUTINE


!-----------------------------------------------------------------------
! Prints to standard output the tabulated differential rate spectrum
! for the current WIMP mass and couplings.
! 
SUBROUTINE WriteSpectrumData(Detector)
  IMPLICIT NONE
  TYPE(DetectorStruct), INTENT(IN) :: Detector
  INTEGER :: KE
  
  ! Columns to print depend on the verbosity level.
  ! For data, only absolute value of verbosity is used.
  
  DO KE = 1,Detector%NE
    WRITE(*,'(A,1(1X,F10.4))',ADVANCE='NO') DATA_PREFIX,                &
        Detector%E(KE)
    SELECT CASE (ABS(VerbosityLevel))
    CASE (:2)
      ! Differential rate
      WRITE(*,'(1(1X,1PG12.5))',ADVANCE='NO')
    CASE (3)
      ! Combined, si, sd
      WRITE(*,'(3(1X,1PG12.5))',ADVANCE='NO')
    CASE (4:)
      ! Reference rate components
      WRITE(*,'(6(1X,1PG12.5))',ADVANCE='NO')
    END SELECT
    WRITE(*,'(A)') ''
  END DO
  
END SUBROUTINE


!-----------------------------------------------------------------------
! Prints to standard output a header containing information regarding
! the events table (tabulated by mass) to follow.
! 
! Optional input arguments:
!   extra_lines     Blank lines (albeit with prefix) to add after
!                   output
! 
SUBROUTINE WriteEventsByMassHeader(WIMP,extra_lines)
  IMPLICIT NONE
  TYPE(WIMPStruct), INTENT(IN) :: WIMP
  INTEGER, INTENT(IN), OPTIONAL :: extra_lines
  REAL*8 :: sigmapSI,sigmanSI,sigmapSD,sigmanSD

  IF ( WIMP%DMtype .EQ. 'SISD' ) THEN
  
    !WRITE(*,'(A)') COMMENT_LINE
    
    ! Will write out the fixed WIMP-nucleon cross-sections.
    ! Negative for negative couplings.
    sigmapSI = FpToSigmapSI(WIMP%m,WIMP%params(1))
    sigmanSI = FnToSigmanSI(WIMP%m,WIMP%params(2))
    sigmapSD = ApToSigmapSD(WIMP%m,WIMP%params(3))
    sigmanSD = AnToSigmanSD(WIMP%m,WIMP%params(4))
    
    ! Description and fixed cross-sections
    IF (VerbosityLevel .GE. 2) THEN
      WRITE(*,'(A)') COMMENT_PREFIX &
          // 'The table below gives the expected spin-independent (SI) and spin-'
      WRITE(*,'(A)') COMMENT_PREFIX &
          // 'dependent (SD) interaction events, tabulated by WIMP mass, for fixed'
      WRITE(*,'(A)') COMMENT_PREFIX &
          // 'WIMP-nucleon cross-sections.  The fixed cross-sections are [pb]:'
      WRITE(*,'(A,1(2X,1PG12.4))') COMMENT_PREFIX &
          // '  sigmapSI (proton SI)  =',CoerceExponent(sigmapSI,2,4)
      WRITE(*,'(A,1(2X,1PG12.4))') COMMENT_PREFIX &
          // '  sigmanSI (neutron SI) =',CoerceExponent(sigmanSI,2,4)
      WRITE(*,'(A,1(2X,1PG12.4))') COMMENT_PREFIX &
          // '  sigmapSD (proton SD)  =',CoerceExponent(sigmapSD,2,4)
      WRITE(*,'(A,1(2X,1PG12.4))') COMMENT_PREFIX &
          // '  sigmanSD (neutron SD) =',CoerceExponent(sigmanSD,2,4)
      WRITE(*,'(A)') COMMENT_PREFIX &
          // 'A negative cross-section means the corresponding WIMP-nucleon coupling'
      WRITE(*,'(A)') COMMENT_PREFIX &
          // 'will be taken to be negative.'
      !WRITE(*,'(A)') COMMENT_PREFIX &
      !    // ''
      WRITE(*,'(A)') COMMENT_PREFIX
    END IF
    
    IF (VerbosityLevel .GE. 2) THEN
      WRITE(*,'(A)') COMMENT_PREFIX &
          // 'The columns below contain the following data:'
    END IF
    
    IF (VerbosityLevel .GE. 2) THEN
      WRITE(*,'(A)') COMMENT_PREFIX &
          // '  mass         WIMP mass [GeV].'
    END IF
    
    IF (VerbosityLevel .GE. 3) THEN
      WRITE(*,'(A)') COMMENT_PREFIX &
          // '  G            WIMP-nucleon couplings for spin-independent (SI) and'
      WRITE(*,'(A)') COMMENT_PREFIX &
          // '               spin-dependent (SD) interactions [GeV^-2].'
    END IF
    
    IF (VerbosityLevel .GE. 2) THEN
      WRITE(*,'(A)') COMMENT_PREFIX &
          // '  events(SI)   Average expected spin-independent events.'
      WRITE(*,'(A)') COMMENT_PREFIX &
          // '  events(SD)   Average expected spin-dependent events.'
    END IF
    
    IF (VerbosityLevel .GE. 2) THEN
      WRITE(*,'(A)') COMMENT_PREFIX
    END IF
    
    IF (PRESENT(extra_lines)) CALL WriteEmptyCommentLines(extra_lines)
  
  ELSE

    WRITE(0,*) 'WARNING: WriteEventsByMassHeader called with WIMP type other than SISD.'

  END IF

END SUBROUTINE


!-----------------------------------------------------------------------
! Prints to standard output a column header for the table of expected
! events tabulated by mass.
! 
SUBROUTINE WriteEventsByMassColumnHeader(Detector)
  IMPLICIT NONE
  TYPE(DetectorStruct), INTENT(IN) :: Detector
  INTEGER :: Keff
  
  ! Columns to print depend on the verbosity level.
  ! For data, only absolute value of verbosity is used.
  
  ! Multiple lines in some cases
  IF (VerbosityLevel .GE. 3) THEN
    ! Mass column
    WRITE(*,'(A,1(1X,A10))',ADVANCE='NO') COMMENT_PREFIX,''
    ! Coupling columns
    IF (VerbosityLevel .GE. 4) THEN
      WRITE(*,'(4(1X,A10))',ADVANCE='NO') '','','',''
    END IF
    ! Events for full range
    WRITE(*,'(1X,A1)',ADVANCE='NO') ' '
    WRITE(*,'(1(1X,A23))',ADVANCE='NO')                                 &
        '----- full range ------'
    ! Events for sub-intervals
    IF (Detector%intervals) THEN
      DO Keff = 1,Detector%Nbins
        WRITE(*,'(1X,A1)',ADVANCE='NO') '|'
        WRITE(*,'(1(1X,A14,I3,A6))',ADVANCE='NO')                       &
            '----- interval',Keff,' -----'
      END DO
    END IF
    WRITE(*,'(A)') ''
  END IF
  
  ! Main column header line below
  
  ! Mass column
  IF (VerbosityLevel .EQ. 1) THEN
    WRITE(*,'(A,1(1X,A10))',ADVANCE='NO') COMMENT_PREFIX,'mass [GeV]'
  ELSE IF (VerbosityLevel .GE. 2) THEN
    WRITE(*,'(A,1(1X,A10))',ADVANCE='NO') COMMENT_PREFIX,'    mass  '
  END IF
  
  ! Coupling columns
  IF (VerbosityLevel .GE. 4) THEN
    WRITE(*,'(4(1X,A10))',ADVANCE='NO')                                 &
        '   GpSI   ','   GnSI   ','   GpSD   ','   GnSD   '
  END IF
  
  ! Events for full range
  IF (VerbosityLevel .GE. 3) WRITE(*,'(1X,A1)',ADVANCE='NO') ' '
  IF (VerbosityLevel .GE. 1) THEN
    WRITE(*,'(2(1X,A11))',ADVANCE='NO') ' events(SI)',' events(SD)'
  END IF
  
  ! Events for sub-intervals
  IF (VerbosityLevel .GE. 3) THEN
    IF (Detector%intervals) THEN
      DO Keff = 1,Detector%Nbins
        WRITE(*,'(1X,A1)',ADVANCE='NO') ' '
        WRITE(*,'(2(1X,A11))',ADVANCE='NO') ' events(SI)',' events(SD)'
      END DO
    END IF
  END IF
  
  IF (VerbosityLevel .GE. 1) WRITE(*,'(A)') ''
  
END SUBROUTINE


!-----------------------------------------------------------------------
! Prints to standard output the mass and expected events for the
! current WIMP.  Used for tabulation of events by WIMP mass.
! 
SUBROUTINE WriteEventsByMassData(Detector,WIMP)
  IMPLICIT NONE
  TYPE(DetectorStruct), INTENT(IN) :: Detector
  TYPE(WIMPStruct), INTENT(IN) :: WIMP
  INTEGER :: Keff
  
  ! Columns to print depend on the verbosity level.
  ! For data, only absolute value of verbosity is used.
  
  ! Mass
  WRITE(*,'(A,1(1X,F10.4))',ADVANCE='NO') DATA_PREFIX,                  &
      CoerceNumber(WIMP%m,10,4)
  
  ! ---Outdated---
  ! Couplings
  !IF (ABS(VerbosityLevel) .GE. 4) THEN
  !  WRITE(*,'(4(1X,1PG10.3))',ADVANCE='NO')                             &
  !      CoerceExponent(WIMP%GpSI,2,3),CoerceExponent(WIMP%GnSI,2,3),    &
  !      CoerceExponent(WIMP%GpSD,2,3),CoerceExponent(WIMP%GnSD,2,3)
  !END IF
  
  
  
  WRITE(*,'(A)') ''
  
END SUBROUTINE


!-----------------------------------------------------------------------
! Prints to standard output a header containing information regarding
! the SI cross-section constraints table (tabulated by mass) to follow.
! 
! Required input arguments:
!   lnp             The logarithm of the p-value to use for the
!                   constraint CL.
!   thetaG          The angle of (Gp,Gn) that is fixed for each
!                   constraint determination.
!   s1,s2           The range of allowed signal expectation values.
! Optional input arguments:
!   extra_lines     Blank lines (albeit with prefix) to add after
!                   output
! 
SUBROUTINE WriteConstraintsSIHeader(lnp,thetaG,s1,s2,extra_lines)
  IMPLICIT NONE
  REAL*8, INTENT(IN) :: lnp,thetaG,s1,s2
  INTEGER, INTENT(IN), OPTIONAL :: extra_lines
  
  !WRITE(*,'(A)') COMMENT_LINE
  
  ! Description
  IF (VerbosityLevel .GE. 2) THEN
    WRITE(*,'(A)') COMMENT_PREFIX &
        // 'The table below gives the range of allowed spin-independent (SI) cross-'
    WRITE(*,'(A)') COMMENT_PREFIX &
        // 'sections.  The allowed cross-sections are those that predict a mean'
    WRITE(*,'(A)') COMMENT_PREFIX &
        // 'signal compatible with the observed number of events and estimated'
    WRITE(*,'(A)') COMMENT_PREFIX &
        // 'background.  The confidence interval (CI) for the mean signal at the'
    WRITE(*,'(A)') COMMENT_PREFIX &
        // 'given confidence level CL is determined using a Poisson distribution'
    WRITE(*,'(A)') COMMENT_PREFIX &
        // 'with Feldman-Cousins ordering; see Feldman & Cousins, Phys. Rev. D 57,'
    WRITE(*,'(A)') COMMENT_PREFIX &
        // '3873 (1998) [physics/9711021].'
    WRITE(*,'(A,1(2X,1PG12.4))') COMMENT_PREFIX &
        // '  1-CL                  =',CoerceExponent(EXP(lnp),2,4)
    WRITE(*,'(A,1(2X,1PG12.4),2X,A,1(2X,1PG12.4))') COMMENT_PREFIX &
        // '  signal events CI      =',CoerceExponent(s1,2,4),'-',CoerceExponent(s2,2,4)
    WRITE(*,'(A)') COMMENT_PREFIX
    
    WRITE(*,'(A)') COMMENT_PREFIX &
        // 'The constraints are determined using a fixed ratio of the two WIMP-nucleon'
    WRITE(*,'(A)') COMMENT_PREFIX &
        // 'couplings Gp & Gn, described by the angle theta s.t. tan(theta) = Gn/Gp.'
    WRITE(*,'(A)') COMMENT_PREFIX &
        // 'The conventional isospin-invariant case Gn=Gp corresponds to theta=PI/4.'
    WRITE(*,'(A,1(2X,F9.5))') COMMENT_PREFIX &
        // '  theta/PI              =',CoerceNumber(thetaG/PI,10,5)
    WRITE(*,'(A)') COMMENT_PREFIX
  END IF
  
  IF (VerbosityLevel .GE. 2) THEN
    WRITE(*,'(A)') COMMENT_PREFIX &
        // 'The columns below contain the following data:'
  END IF
  
  IF (VerbosityLevel .GE. 2) THEN
    WRITE(*,'(A)') COMMENT_PREFIX &
        // '  mass         WIMP mass [GeV].'
  END IF
  
  IF (VerbosityLevel .GE. 2) THEN
    WRITE(*,'(A)') COMMENT_PREFIX &
        // '  sigmapSI     The WIMP-proton spin-independent cross-section lower and'
    WRITE(*,'(A)') COMMENT_PREFIX &
        // '               upper limits [pb].'
  END IF
  
  IF (VerbosityLevel .GE. 3) THEN
    WRITE(*,'(A)') COMMENT_PREFIX &
        // '  sigmanSI     The WIMP-neutron spin-independent cross-section lower and'
    WRITE(*,'(A)') COMMENT_PREFIX &
        // '               upper limits [pb].'
  END IF
  
  IF (VerbosityLevel .GE. 2) THEN
    WRITE(*,'(A)') COMMENT_PREFIX
  END IF
  
  IF (PRESENT(extra_lines)) CALL WriteEmptyCommentLines(extra_lines)
  
END SUBROUTINE


!-----------------------------------------------------------------------
! Prints to standard output a column header for the table of SI cross-
! section constraints tabulated by mass.
! 
SUBROUTINE WriteConstraintsSIColumnHeader()
  IMPLICIT NONE
  
  ! Columns to print depend on the verbosity level.
  ! For data, only absolute value of verbosity is used.
  
  ! Mass column
  IF (VerbosityLevel .EQ. 1) THEN
    WRITE(*,'(A,1(1X,A10))',ADVANCE='NO') COMMENT_PREFIX,'mass [GeV]'
  ELSE IF (VerbosityLevel .GE. 2) THEN
    WRITE(*,'(A,1(1X,A10))',ADVANCE='NO') COMMENT_PREFIX,'    mass  '
  END IF
  
  ! WIMP-proton cross-section columns
  IF (VerbosityLevel .EQ. 1) THEN
    WRITE(*,'(1(1X,A25))',ADVANCE='NO') ' --sigmapSI range [pb]-- '
  ELSE IF (VerbosityLevel .GE. 2) THEN
    WRITE(*,'(1(1X,A25))',ADVANCE='NO') ' --- sigmapSI range ---- '
  END IF
  
  ! WIMP-neutron cross-section column
  IF (VerbosityLevel .GE. 3) THEN
    WRITE(*,'(1(1X,A25))',ADVANCE='NO') ' --- sigmanSI range ---- '
  END IF
  
  IF (VerbosityLevel .GE. 1) WRITE(*,'(A)') ''
  
END SUBROUTINE


!-----------------------------------------------------------------------
! Prints to standard output the mass and SI cross-section limit for the
! current WIMP.  Used for tabulation of limits by WIMP mass.
! 
! Input argument:
!   s1,s2       The range of allowed signal expectation values.
! 
SUBROUTINE WriteConstraintsSIData(s1,s2,Detector,WIMP)
  IMPLICIT NONE
  TYPE(DetectorStruct), INTENT(IN) :: Detector
  TYPE(WIMPStruct), INTENT(IN) :: WIMP
  REAL*8, INTENT(IN) :: s1,s2
  REAL*8 :: mu,x1,x2,sigmapSI,sigmanSI
  
  IF ( WIMP%DMtype .EQ. 'SIonly' ) THEN

    ! Need scale factors x s.t. sigma -> x*sigma gives desired
    ! number of events.
    CALL DDCalc_GetRates(Detector,signal_si=mu)
    ! Empty set case
    IF (s2 .EQ. 0d0) THEN
      x1 = 0d0
      x2 = 0d0
    ! General case
    ELSE IF (mu .GT. 0d0) THEN
      x1 = s1/mu
      x2 = s2/mu
    ! No events case (at any scale)
    ELSE
      x1 = 0d0
      x2 = HUGE(1d0)
    END IF
  
    ! Cross-sections (multiply by x for limit)
    sigmapSI = FpToSigmapSI(WIMP%m,WIMP%params(1))
    sigmanSI = FnToSigmanSI(WIMP%m,WIMP%params(2))
  
    ! Columns to print depend on the verbosity level.
    ! For data, only absolute value of verbosity is used.
  
    ! Mass
    WRITE(*,'(A,1(1X,F10.4))',ADVANCE='NO') DATA_PREFIX,                  &
        CoerceNumber(WIMP%m,10,4)
  
    ! WIMP-proton cross-section
    WRITE(*,'(2(1X,1PG12.5))',ADVANCE='NO')                               &
        CoerceExponent(x1*sigmapSI,2,5),CoerceExponent(x2*sigmapSI,2,5)
  
    ! WIMP-neutron cross-section
    IF (ABS(VerbosityLevel) .GE. 3) THEN
      WRITE(*,'(2(1X,1PG12.5))',ADVANCE='NO')                             &
          CoerceExponent(x1*sigmanSI,2,5),CoerceExponent(x2*sigmanSI,2,5)
    END IF
  
    WRITE(*,'(A)') ''

  ELSE

    WRITE(0,*) 'WARNING: WriteConstraintsSIData called with WIMP type other than SIonly.'

  END IF

  
END SUBROUTINE


!-----------------------------------------------------------------------
! Prints to standard output a header containing information regarding
! the SD cross-section constraints table (tabulated by mass) to follow.
! 
! Required input arguments:
!   lnp             The logarithm of the p-value to use for the
!                   constraint CL.
!   thetaG          The angle of (Gp,Gn) that is fixed for each
!                   constraint determination.
!   s1,s2           The range of allowed signal expectation values.
! Optional input arguments:
!   extra_lines     Blank lines (albeit with prefix) to add after
!                   output
! 
SUBROUTINE WriteConstraintsSDHeader(lnp,thetaG,s1,s2,extra_lines)
  IMPLICIT NONE
  REAL*8, INTENT(IN) :: lnp,thetaG,s1,s2
  INTEGER, INTENT(IN), OPTIONAL :: extra_lines
  
  !WRITE(*,'(A)') COMMENT_LINE
  
  ! Description
  IF (VerbosityLevel .GE. 2) THEN
    WRITE(*,'(A)') COMMENT_PREFIX &
        // 'The table below gives the range of allowed spin-dependent (SD) cross-'
    WRITE(*,'(A)') COMMENT_PREFIX &
        // 'sections.  The allowed cross-sections are those that predict a mean'
    WRITE(*,'(A)') COMMENT_PREFIX &
        // 'signal compatible with the observed number of events and estimated'
    WRITE(*,'(A)') COMMENT_PREFIX &
        // 'background.  The confidence interval (CI) for the mean signal at the'
    WRITE(*,'(A)') COMMENT_PREFIX &
        // 'given confidence level CL is determined using a Poisson distribution'
    WRITE(*,'(A)') COMMENT_PREFIX &
        // 'with Feldman-Cousins ordering; see Feldman & Cousins, Phys. Rev. D 57,'
    WRITE(*,'(A)') COMMENT_PREFIX &
        // '3873 (1998) [physics/9711021].'
    WRITE(*,'(A,1(2X,1PG12.4))') COMMENT_PREFIX &
        // '  1-CL                  =',CoerceExponent(EXP(lnp),2,4)
    WRITE(*,'(A,1(2X,1PG12.4),2X,A,1(2X,1PG12.4))') COMMENT_PREFIX &
        // '  signal events CI      =',CoerceExponent(s1,2,4),'-',CoerceExponent(s2,2,4)
    WRITE(*,'(A)') COMMENT_PREFIX
    
    WRITE(*,'(A)') COMMENT_PREFIX &
        // 'The constraints are determined using a fixed ratio of the two WIMP-nucleon'
    WRITE(*,'(A)') COMMENT_PREFIX &
        // 'couplings Gp & Gn, described by the angle theta s.t. tan(theta) = Gn/Gp.'
    WRITE(*,'(A)') COMMENT_PREFIX &
        // 'The isospin-invariant case Gn=Gp corresponds to theta=PI/4, proton-only'
    WRITE(*,'(A)') COMMENT_PREFIX &
        // 'coupling is theta=0, and neutron-only coupling is theta=PI/2.'
    WRITE(*,'(A,1(2X,F9.5))') COMMENT_PREFIX &
        // '  theta/PI              =',CoerceNumber(thetaG/PI,10,5)
    WRITE(*,'(A)') COMMENT_PREFIX
  END IF
  
  IF (VerbosityLevel .GE. 2) THEN
    WRITE(*,'(A)') COMMENT_PREFIX &
        // 'The columns below contain the following data:'
  END IF
  
  IF (VerbosityLevel .GE. 2) THEN
    WRITE(*,'(A)') COMMENT_PREFIX &
        // '  mass         WIMP mass [GeV].'
  END IF
  
  IF (VerbosityLevel .GE. 2) THEN
    WRITE(*,'(A)') COMMENT_PREFIX &
        // '  sigmapSD     The WIMP-proton spin-dependent cross-section lower and'
    WRITE(*,'(A)') COMMENT_PREFIX &
        // '               upper limits [pb].'
  END IF
  
  IF (VerbosityLevel .GE. 3) THEN
    WRITE(*,'(A)') COMMENT_PREFIX &
        // '  sigmanSD     The WIMP-neutron spin-dependent cross-section lower and'
    WRITE(*,'(A)') COMMENT_PREFIX &
        // '               upper limits [pb].'
  END IF
  
  IF (VerbosityLevel .GE. 2) THEN
    WRITE(*,'(A)') COMMENT_PREFIX
  END IF
  
  IF (PRESENT(extra_lines)) CALL WriteEmptyCommentLines(extra_lines)
  
END SUBROUTINE


!-----------------------------------------------------------------------
! Prints to standard output a column header for the table of SD cross-
! section constraints tabulated by mass.
! 
SUBROUTINE WriteConstraintsSDColumnHeader()
  IMPLICIT NONE
  
  ! Columns to print depend on the verbosity level.
  ! For data, only absolute value of verbosity is used.
  
  ! Mass column
  IF (VerbosityLevel .EQ. 1) THEN
    WRITE(*,'(A,1(1X,A10))',ADVANCE='NO') COMMENT_PREFIX,'mass [GeV]'
  ELSE IF (VerbosityLevel .GE. 2) THEN
    WRITE(*,'(A,1(1X,A10))',ADVANCE='NO') COMMENT_PREFIX,'    mass  '
  END IF
  
  ! WIMP-proton cross-section columns
  IF (VerbosityLevel .EQ. 1) THEN
    WRITE(*,'(1(1X,A25))',ADVANCE='NO') ' --sigmapSD range [pb]-- '
  ELSE IF (VerbosityLevel .GE. 2) THEN
    WRITE(*,'(1(1X,A25))',ADVANCE='NO') ' --- sigmapSD range ---- '
  END IF
  
  ! WIMP-neutron cross-section column
  IF (VerbosityLevel .GE. 3) THEN
    WRITE(*,'(1(1X,A25))',ADVANCE='NO') ' --- sigmanSD range ---- '
  END IF
  
  IF (VerbosityLevel .GE. 1) WRITE(*,'(A)') ''
  
END SUBROUTINE


!-----------------------------------------------------------------------
! Prints to standard output the mass and SD cross-section limit for the
! current WIMP.  Used for tabulation of limits by WIMP mass.
! 
! Input argument:
!   s1,s2       The range of allowed signal expectation values.
! 
SUBROUTINE WriteConstraintsSDData(s1,s2,Detector,WIMP)
  IMPLICIT NONE
  TYPE(DetectorStruct), INTENT(IN) :: Detector
  TYPE(WIMPStruct), INTENT(IN) :: WIMP
  REAL*8, INTENT(IN) :: s1,s2
  REAL*8 :: mu,x1,x2,sigmapSD,sigmanSD
  
  IF ( WIMP%DMtype .EQ. 'SDonly' ) THEN

    ! Need scale factors x s.t. sigma -> x*sigma gives desired
    ! number of events.
    CALL DDCalc_GetRates(Detector,signal_sd=mu)
    ! Empty set case
    IF (s2 .EQ. 0d0) THEN
      x1 = 0d0
      x2 = 0d0
    ! General case
    ELSE IF (mu .GT. 0d0) THEN
      x1 = s1/mu
      x2 = s2/mu
    ! No events case (at any scale)
    ELSE
      x1 = 0d0
      x2 = HUGE(1d0)
    END IF
  
    ! Cross-sections (multiply by x for limit)
    sigmapSD = ApToSigmapSD(WIMP%m,WIMP%params(1))
    sigmanSD = AnToSigmanSD(WIMP%m,WIMP%params(2))
 
    ! Columns to print depend on the verbosity level.
    ! For data, only absolute value of verbosity is used.
  
    ! Mass
    WRITE(*,'(A,1(1X,F10.4))',ADVANCE='NO') DATA_PREFIX,                  &
        CoerceNumber(WIMP%m,10,4)
  
    ! WIMP-proton cross-section
    WRITE(*,'(2(1X,1PG12.5))',ADVANCE='NO')                               &
        CoerceExponent(x1*sigmapSD,2,5),CoerceExponent(x2*sigmapSD,2,5)
  
    ! WIMP-neutron cross-section
    IF (ABS(VerbosityLevel) .GE. 3) THEN
      WRITE(*,'(2(1X,1PG12.5))',ADVANCE='NO')                             &
          CoerceExponent(x1*sigmanSD,2,5),CoerceExponent(x2*sigmanSD,2,5)
    END IF
  
    WRITE(*,'(A)') ''
  
  ELSE

    WRITE(0,*) 'WARNING: WriteConstraintsSDData called with WIMP type other than SDonly.'

  END IF

END SUBROUTINE


!-----------------------------------------------------------------------
! Prints to standard output a header containing information regarding
! the SI cross-section limits table (tabulated by mass) to follow.
! 
! Required input arguments:
!   lnp             The logarithm of the p-value to use for the limit
!                   CL.
!   thetaG          The angle of (Gp,Gn) that is fixed for each limit
!                   determination.
! Optional input arguments:
!   extra_lines     Blank lines (albeit with prefix) to add after
!                   output
! 
SUBROUTINE WriteLimitsSIHeader(lnp,thetaG,Detector,extra_lines)
  IMPLICIT NONE
  TYPE(DetectorStruct), INTENT(IN) :: Detector
  REAL*8, INTENT(IN) :: lnp,thetaG
  INTEGER, INTENT(IN), OPTIONAL :: extra_lines
  
  !WRITE(*,'(A)') COMMENT_LINE
  
  ! Description
  IF (VerbosityLevel .GE. 2) THEN
    WRITE(*,'(A)') COMMENT_PREFIX &
        // 'The table below gives the upper limit on spin-independent (SI) cross-'
    WRITE(*,'(A)') COMMENT_PREFIX &
        // 'section(s) that are not excluded.  Cross-sections are excluded if their'
    WRITE(*,'(A)') COMMENT_PREFIX &
        // 'p-value is smaller than the given p-value, where the p-value is'
    IF (Detector%intervals .AND. (Detector%Nbins .EQ. Detector%Nevents(0)+1)) THEN
      WRITE(*,'(A)') COMMENT_PREFIX &
        // 'determined using the maximum gap method; see Yellin, Phys. Rev. D 66,'
      WRITE(*,'(A)') COMMENT_PREFIX &
        // '032005 (2002) [physics/0203002].'
    ELSE
      WRITE(*,'(A)') COMMENT_PREFIX &
        // 'determined using the Poisson distribution (signal only: no background).'
    END IF
    WRITE(*,'(A,1(2X,1PG12.4))') COMMENT_PREFIX &
        // '  p-value               =',CoerceExponent(EXP(lnp),2,4)
    WRITE(*,'(A)') COMMENT_PREFIX
    
    WRITE(*,'(A)') COMMENT_PREFIX &
        // 'The limit is determined using a fixed ratio of the two WIMP-nucleon'
    WRITE(*,'(A)') COMMENT_PREFIX &
        // 'couplings Gp & Gn, described by the angle theta s.t. tan(theta) = Gn/Gp.'
    WRITE(*,'(A)') COMMENT_PREFIX &
        // 'The conventional isospin-invariant case Gn=Gp corresponds to theta=PI/4.'
    WRITE(*,'(A,1(2X,F9.5))') COMMENT_PREFIX &
        // '  theta/PI              =',CoerceNumber(thetaG/PI,10,5)
    WRITE(*,'(A)') COMMENT_PREFIX
  END IF
  
  IF (VerbosityLevel .GE. 2) THEN
    WRITE(*,'(A)') COMMENT_PREFIX &
        // 'The columns below contain the following data:'
  END IF
  
  IF (VerbosityLevel .GE. 2) THEN
    WRITE(*,'(A)') COMMENT_PREFIX &
        // '  mass         WIMP mass [GeV].'
  END IF
  
  IF (VerbosityLevel .GE. 2) THEN
    WRITE(*,'(A)') COMMENT_PREFIX &
        // '  sigmapSI     The WIMP-proton spin-independent cross-section upper limit [pb].'
  END IF
  
  IF (VerbosityLevel .GE. 3) THEN
    WRITE(*,'(A)') COMMENT_PREFIX &
        // '  sigmanSI     The WIMP-neutron spin-independent cross-section upper limit [pb].'
  END IF
  
  IF (VerbosityLevel .GE. 2) THEN
    WRITE(*,'(A)') COMMENT_PREFIX
  END IF
  
  IF (PRESENT(extra_lines)) CALL WriteEmptyCommentLines(extra_lines)
  
END SUBROUTINE


!-----------------------------------------------------------------------
! Prints to standard output a column header for the table of SI cross-
! section limits tabulated by mass.
! 
SUBROUTINE WriteLimitsSIColumnHeader()
  IMPLICIT NONE
  
  ! Columns to print depend on the verbosity level.
  ! For data, only absolute value of verbosity is used.
  
  ! Mass column
  IF (VerbosityLevel .EQ. 1) THEN
    WRITE(*,'(A,1(1X,A10))',ADVANCE='NO') COMMENT_PREFIX,'mass [GeV]'
  ELSE IF (VerbosityLevel .GE. 2) THEN
    WRITE(*,'(A,1(1X,A10))',ADVANCE='NO') COMMENT_PREFIX,'    mass  '
  END IF
  
  ! WIMP-proton cross-section column
  IF (VerbosityLevel .EQ. 1) THEN
    WRITE(*,'(1(1X,A12))',ADVANCE='NO') 'sigmapSI[pb]'
  ELSE IF (VerbosityLevel .GE. 2) THEN
    WRITE(*,'(1(1X,A12))',ADVANCE='NO') '  sigmapSI  '
  END IF
  
  ! WIMP-neutron cross-section column
  IF (VerbosityLevel .GE. 3) THEN
    WRITE(*,'(1(1X,A12))',ADVANCE='NO') '  sigmanSI  '
  END IF
  
  IF (VerbosityLevel .GE. 1) WRITE(*,'(A)') ''
  
END SUBROUTINE


!-----------------------------------------------------------------------
! Prints to standard output the mass and SI cross-section limit for the
! passed WIMP and detector.  Used for tabulation of limits by WIMP mass.
! 
! Input argument:
!   lnp         Logarithm of the p-value for the limit CL
! 
SUBROUTINE WriteLimitsSIData(lnp,Detector,WIMP)
  IMPLICIT NONE
  TYPE(DetectorStruct), INTENT(IN) :: Detector
  TYPE(WIMPStruct), INTENT(IN) :: WIMP
  REAL*8, INTENT(IN) :: lnp
  REAL*8 :: x,sigmapSI,sigmanSI
  
  IF ( WIMP%DMtype .EQ. 'SIonly' ) THEN

    ! Need scale factor
    x = DDCalc_ScaleToPValue(Detector,lnp=lnp)
  
    ! Cross-sections (multiply by x for limit)
    sigmapSI = FpToSigmapSI(WIMP%m,WIMP%params(1))
    sigmanSI = FnToSigmanSI(WIMP%m,WIMP%params(1))
  
    ! Columns to print depend on the verbosity level.
    ! For data, only absolute value of verbosity is used.
  
    ! Mass
    WRITE(*,'(A,1(1X,F10.4))',ADVANCE='NO') DATA_PREFIX,                  &
        CoerceNumber(WIMP%m,10,4)
  
    ! WIMP-proton cross-section
    WRITE(*,'(1(1X,1PG12.5))',ADVANCE='NO')                               &
        CoerceExponent(x*sigmapSI,2,5)
  
    ! WIMP-neutron cross-section
    IF (ABS(VerbosityLevel) .GE. 3) THEN
      WRITE(*,'(1(1X,1PG12.5))',ADVANCE='NO')                             &
          CoerceExponent(x*sigmanSI,2,5)
    END IF
  
    WRITE(*,'(A)') ''
  
  ELSE

    WRITE(0,*) 'WARNING: WriteLimitsSIData called with WIMP type other than SIonly.'

  END IF

END SUBROUTINE


!-----------------------------------------------------------------------
! Prints to standard output a header containing information regarding
! the SD cross-section limits table (tabulated by mass) to follow.
! 
! Required input arguments:
!   lnp             The logarithm of the p-value to use for the limit
!                   CL.
!   thetaG          The angle of (Gp,Gn) that is fixed for each limit
!                   determination.
! Optional input arguments:
!   extra_lines     Blank lines (albeit with prefix) to add after
!                   output
! 
SUBROUTINE WriteLimitsSDHeader(lnp,thetaG,Detector,extra_lines)
  IMPLICIT NONE
  TYPE(DetectorStruct), INTENT(IN) :: Detector
  REAL*8, INTENT(IN) :: lnp,thetaG
  INTEGER, INTENT(IN), OPTIONAL :: extra_lines
  
  !WRITE(*,'(A)') COMMENT_LINE
  
  ! Description
  IF (VerbosityLevel .GE. 2) THEN
    WRITE(*,'(A)') COMMENT_PREFIX &
        // 'The table below gives the upper limit on spin-dependent (SD) cross-'
    WRITE(*,'(A)') COMMENT_PREFIX &
        // 'section(s) that are not excluded.  Cross-sections are excluded if their'
    WRITE(*,'(A)') COMMENT_PREFIX &
        // 'p-value is smaller than the given p-value, where the p-value is'
    IF (Detector%intervals .AND. (Detector%Nbins .EQ. Detector%Nevents(0)+1)) THEN
      WRITE(*,'(A)') COMMENT_PREFIX &
        // 'determined using the maximum gap method; see Yellin, Phys. Rev. D 66,'
      WRITE(*,'(A)') COMMENT_PREFIX &
        // '032005 (2002) [physics/0203002].'
    ELSE
      WRITE(*,'(A)') COMMENT_PREFIX &
        // 'determined using the Poisson distribution (signal only: no background'
      WRITE(*,'(A)') COMMENT_PREFIX &
        // ').'
    END IF
    WRITE(*,'(A,1(2X,1PG12.4))') COMMENT_PREFIX &
        // '  p-value               =',CoerceExponent(EXP(lnp),2,4)
    WRITE(*,'(A)') COMMENT_PREFIX
    
    WRITE(*,'(A)') COMMENT_PREFIX &
        // 'The limit is determined using a fixed ratio of the two WIMP-nucleon'
    WRITE(*,'(A)') COMMENT_PREFIX &
        // 'couplings Gp & Gn, described by the angle theta s.t. tan(theta) = Gn/Gp.'
    WRITE(*,'(A)') COMMENT_PREFIX &
        // 'The isospin-invariant case Gn=Gp corresponds to theta=PI/4, proton-only'
    WRITE(*,'(A)') COMMENT_PREFIX &
        // 'coupling is theta=0, and neutron-only coupling is theta=PI/2.'
    WRITE(*,'(A,1(2X,F9.5))') COMMENT_PREFIX &
        // '  theta/PI              =',CoerceNumber(thetaG/PI,10,5)
    WRITE(*,'(A)') COMMENT_PREFIX
    
    !WRITE(*,'(A)') COMMENT_PREFIX &
    !    // 'NOTE: SPIN-DEPENDENT FORM FACTORS NOT IMPLEMENTED.'
    !WRITE(*,'(A)') COMMENT_PREFIX
  END IF
  
  IF (VerbosityLevel .GE. 2) THEN
    WRITE(*,'(A)') COMMENT_PREFIX &
        // 'The columns below contain the following data:'
  END IF
  
  IF (VerbosityLevel .GE. 2) THEN
    WRITE(*,'(A)') COMMENT_PREFIX &
        // '  mass         WIMP mass [GeV].'
  END IF
  
  IF (VerbosityLevel .GE. 2) THEN
    WRITE(*,'(A)') COMMENT_PREFIX &
        // '  sigmapSD     The WIMP-proton spin-dependent cross-section upper limit [pb].'
  END IF
  
  IF (VerbosityLevel .GE. 3) THEN
    WRITE(*,'(A)') COMMENT_PREFIX &
        // '  sigmanSD     The WIMP-neutron spin-dependent cross-section upper limit [pb].'
  END IF
  
  IF (VerbosityLevel .GE. 2) THEN
    WRITE(*,'(A)') COMMENT_PREFIX
  END IF
  
  IF (PRESENT(extra_lines)) CALL WriteEmptyCommentLines(extra_lines)
  
END SUBROUTINE


!-----------------------------------------------------------------------
! Prints to standard output a column header for the table of SD cross-
! section limits tabulated by mass.
! 
SUBROUTINE WriteLimitsSDColumnHeader()
  IMPLICIT NONE
  
  ! Columns to print depend on the verbosity level.
  ! For data, only absolute value of verbosity is used.
  
  ! Mass column
  IF (VerbosityLevel .EQ. 1) THEN
    WRITE(*,'(A,1(1X,A10))',ADVANCE='NO') COMMENT_PREFIX,'mass [GeV]'
  ELSE IF (VerbosityLevel .GE. 2) THEN
    WRITE(*,'(A,1(1X,A10))',ADVANCE='NO') COMMENT_PREFIX,'    mass  '
  END IF
  
  ! WIMP-proton cross-section column
  IF (VerbosityLevel .EQ. 1) THEN
    WRITE(*,'(1(1X,A12))',ADVANCE='NO') 'sigmapSD[pb]'
  ELSE IF (VerbosityLevel .GE. 2) THEN
    WRITE(*,'(1(1X,A12))',ADVANCE='NO') '  sigmapSD  '
  END IF
  
  ! WIMP-neutron cross-section column
  IF (VerbosityLevel .GE. 3) THEN
    WRITE(*,'(1(1X,A12))',ADVANCE='NO') '  sigmanSD  '
  END IF
  
  IF (VerbosityLevel .GE. 1) WRITE(*,'(A)') ''
  
END SUBROUTINE


!-----------------------------------------------------------------------
! Prints to standard output the mass and SD cross-section limit for the
! passed WIMP and detector.  Used for tabulation of limits by WIMP mass.
! 
! Input argument:
!   lnp         Logarithm of the p-value for the limit CL
! 
SUBROUTINE WriteLimitsSDData(lnp,Detector,WIMP)
  IMPLICIT NONE
  TYPE(DetectorStruct), INTENT(IN) :: Detector
  TYPE(WIMPStruct), INTENT(IN) :: WIMP
  REAL*8, INTENT(IN) :: lnp
  REAL*8 :: x,sigmapSD,sigmanSD
  
  IF ( WIMP%DMtype .EQ. 'SDonly' ) THEN

    ! Need scale factor
    x = DDCalc_ScaleToPValue(Detector,lnp=lnp)
  
    ! Cross-sections (multiply by x for limit)
    sigmapSD = ApToSigmapSD(WIMP%m,WIMP%params(1))
    sigmanSD = AnToSigmanSD(WIMP%m,WIMP%params(2))
  
    ! Columns to print depend on the verbosity level.
    ! For data, only absolute value of verbosity is used.
  
    ! Mass
    WRITE(*,'(A,1(1X,F10.4))',ADVANCE='NO') DATA_PREFIX,                  &
        CoerceNumber(WIMP%m,10,4)
  
    ! WIMP-proton cross-section
    WRITE(*,'(1(1X,1PG12.5))',ADVANCE='NO')                               &
        CoerceExponent(x*sigmapSD,2,5)
  
    ! WIMP-neutron cross-section
    IF (ABS(VerbosityLevel) .GE. 3) THEN
      WRITE(*,'(1(1X,1PG12.5))',ADVANCE='NO')                             &
          CoerceExponent(x*sigmanSD,2,5)
    END IF
  
    WRITE(*,'(A)') ''

  ELSE

    WRITE(0,*) 'WARNING: WriteLimitsSDData called with WIMP type other than SDonly.'

  END IF

  
END SUBROUTINE


ENDMODULE
