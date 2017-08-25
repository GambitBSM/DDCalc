DDCalc
======

The DDCalc software package is a set of routines and a frontend for
doing various dark matter direct detection calculations for several
experimental results.  Most notably, these include likelihood
constraints: either Poisson-based with Feldman-Cousins ordering
(Feldman & Cousins 1998), or exclusion limits  based on the maximum
gap method (Yellin 2002).

A full description of this package and the physics framework behind
it can be found in the [GAMBIT DarkBit paper](DarkBitPaper.pdf):

T Bringmann, J Conrad, JM Cornell, LA Dal, J Edsjö, B Farmer,
F Kahlhoefer, A Kvellestad, A Putze, C Savage, P Scott, C Weniger,
M White & S Wild, DarkBit: A GAMBIT module for computing dark matter
observables and likelihoods, EPJC submitted (2017), arXiv:1705.07920.

If you write a paper that uses DDCalc, please cite this paper.


Experiments
--

The experimental results implemented in DDCalc at this time are:

  * Xenon1T 2017: 1042*34.2 kg-days (Xe)  
    E. Aprile et al., [arxiv:1705.06655]

  * LUX 2016: 33500 kg-days (Xe)  
    D.S. Akerib et al., PRL 118, 021303 (2017) [arxiv:1608.07648]

  * PandaX 2016: 34*224.6 kg-days (Xe)  
    A. Tan et al., PRL 117, 121303 (2016) [arxiv:1607.07400]

  * LUX 2015: 118*85.3 kg-days (Xe)  
    D.S. Akerib et al., PRL 116, 161301 (2016) [arXiv:1512.03506]

  * LUX 2013: 118*85.3 kg-days (Xe)  
    D.S. Akerib et al., PRL 112, 091303 (2014) [arxiv:1310.8214]

  * XENON100 2012: 34*224.6 kg-days (Xe)  
    E. Aprile et al., PRL 109, 181301 (2012) [arxiv:1207.5988]

  * PICO: 60 (CF3), 2L (C3F8) and 60 (C3F8)
    C. Amole et al., PRD 93, 052014 (2016) [arXiv:1510.07754]  
    C. Amole et al., PRD 93, 061101 (2016) [arXiv:1601.03729] 
    C. Amole et al., arXiv:1702.07666 

  * SuperCDMS 2014: 577 kg-days (Ge)  
    R. Agnese et al., PRL 112, 241302 (2014) [arxiv:1402.7137]

  * SIMPLE 2014: 18.24 kg-days (C2ClF5)  
    M. Felizardo et al., PRD 89, 072013 (2014) [arXiv:1404.4309]

  * DARWIN projections (argon & xenon)  
    J. Conrad et al., in preparation.


Citation Policy
--

When using DDCalc, the following publications must be cited; they can be
found in the 'DDCalc.bib' file.

The DDCalc software itself is presented in the GAMBIT DarkBit paper:

T Bringmann, J Conrad, JM Cornell, LA Dal, J Edsjö, B Farmer,
F Kahlhoefer, A Kvellestad, A Putze, C Savage, P Scott, C Weniger,
M White & S Wild, DarkBit: A GAMBIT module for computing dark matter
observables and likelihoods, EPJC submitted (2017), arXiv:1705.07920.

Constraints are generated using the Feldman-Cousins method for a
likelihood-based analysis with background subtraction, or using the
Yellin maximum gap method for a conservative exclusion limit with no
background subtraction.  If you use one of those constraints,
please cite the associated paper:

  G. Feldman & R. Cousins, PRD 57, 3873 (1998) [physics/9711021]  
    http://inspirehep.net/record/454197  
  S. Yellin, PRD 66, 032005 (2002) [physics/0203002]  
    http://inspirehep.net/record/583791  

Efficiencies for liquid noble gas time-projection chamber (TPC)
experiments can be generated using the TPCMC software, which itself uses
the NEST model.  If using any TPC experimental results that rely on
TPCMC-generated efficiencies, these papers should be cited:

  NEST:  
    M. Szydagis et al., JINST 6, P10002 (2011) [arxiv:1106.1613]  
      http://inspirehep.net/record/913031  
    M. Szydagis et al., JINST 8, C10003 (2013) [arxiv:1307.6601]  
      http://inspirehep.net/record/1244479


If using DDCalc for calculations involving an  experimental result,
please cite the appropriate one of the experimental papers listed above.


Compiling
---------

DDCalc is written as a Fortran 95 module.  It supports both the  Intel
and GNU compilers.

To build, simply do

  make all

This will generate the following:

  * DDCalc_run
    The main program for performing various rate and likelihood
    calculations.

  * DDCalc.o, DDCalc.a, DDCalc.so
    Object, static library, and shared object files for linking
    with other software.

  * DDCalc_exampleF, DDCalc_exampleC
    Simple example programs with source code showing how to use
    the DDCalc routines externally from either Fortran or C++.


Command-line usage
---------------

The program 'DDCalc_run' is a frontend to the DDCalc routines.
It can perform a variety of calculations.  Full details of its
usage can be found by running:

  ./DDCalc_run --help

or by looking at the 'doc/DDCalc_run.help' file.

As a quick example, try:

  ./DDCalc_run --Xenon1T-2017 --limits-SI --verbosity=3

This will use the Xenon1T results (--Xenon1T-2017) to generate
upper limits on the spin-independent WIMP-nucleon cross-section
using the maximum gap method (--limits-SI). The last flag
(--verbosity=3) gives more detailed output.  Likelihood constraints
can be generated using '--constraints-SI' instead.  Replace 'SI'
with 'SD' for spin-dependent limits.


Usage as a library
---------------

The routines are designed to be accessible to other software packages.
The example source file 'examples/DDCalc_exampleF.f90' gives details
on the available routines and how they must be called from Fortran.

For ease of use, a C++ interface is provided in the 'include/DDCalc.hpp'
header file.  This defines C++ functions of the same names and
signatures as in Fortran.  The 'examples/DDCalc_exampleC.cpp'
source file provides an example of how to call/use DDCalc routines
from C++.


Off the beaten track: momentum-dependent cross-sections
-------------------------------------------------------

A variety of operators exist to describe WIMP-fermion scattering, and many
lead to a non-trivial velocity- or momentum-dependence of the scattering
cross-section. This requires direct search limits to be re-evaluated from
first principles, since integrals over velocity or momentum must now take
account of the changed functional form. Since DDCalc derives dark matter
constraints from scratch rather than using digitised limits from experimental
results derived using simple operators, it is in principle capable of
supporting non-trivial velocity and momentum dependence. Although the
current release does not support this behaviour by default, we here provide
some instructions for the interested for adapting the code, using an example
of a spin-independent cross-section with an extra q^2/m_DM factor which would
arise, to take one example, from a pseudoscalar WIMP-nucleon interaction
with a fermionic WIMP (see Beniwal et al, PRD 2016, arXiv:1512.06458).

The extra factor can be included by making the following modifications to
the DDcalc source code:

1. Create an array of size [1:NE,1:Niso] in the DetectorStruct definition
in DDTypes.f90 that can hold momentum transfer values in units of GeV:

  ! Array of size [1:NE,1:Niso] containing the momenta transfer
  ! corresponding to the tabulated energies and target isotopes [GeV].
  REAL*8, ALLOCATABLE :: q(:,:)

2. Initialise the momentum array in the SetDetector subroutine of
DDDetectors.f90, by adding the following lines immediately before the
definition of the weighted form factors:

  IF ((iso_change .OR. E_change) .AND. (DP%Niso .GE. 0) .AND. (DP%NE .GE. 0)) THEN
    IF (ALLOCATED(DP%q)) DEALLOCATE(DP%q)
    ALLOCATE(DP%q(DP%NE,DP%Niso))
    DO Kiso = 1,DP%Niso
      DP%q(:,Kiso) = EToQ(DP%E,DP%Miso(Kiso))
    END DO
  END IF

3. In the weighted form factor definitions that follow the previous modification,
use DP\%q(:,Kiso) in place of EToQ when calling the CalcWSI and CalcWSD routines.

4. Finally, modify the cross-section for the scattering process to include the
extra q^2/m_DM factor. The relevant code lives in the CalcRates subroutine of
DDRates.f90. To modify the spin independent cross-section, add an extra factor:

  * (DP%q(:,Kiso) / (2 * WIMP%m))**2

to the definitions of DP\%dRdEsi0(+1,:), DP\%dRdEsi0( 0,:) and DP\%dRdEsi0(-1,:).



Contact
-------

Any questions, comments, or suggestions should be directed to the GAMBIT
Dark Matter Workgroup, via ddcalc@projects.hepforge.org.


Contributors
------------

The following have contributed to the development and
testing of the DDCalc package:

  * Felix Kahlhoefer, DESY  
    LUX 2015/16, PICO, Panda X, Xenon1T

  * Lauren Hsu, FermiLab  
    SuperCDMS.

  * Miguel Pato, Stockholm University  
    Calibration/tuning of LUX parameters.

  * Chris Savage  
    Most fortran procedures.

  * Andre Scaffidi, Univ of Adelaide  
    SD form factors, C++ interface and example, testing.

  * Pat Scott, Imperial College London  
    C++ interface and example, testing, fortran structure, bug fixes, writeup.

  * Martin White, Univ of Adelaide  
    C++ interface and example, testing, writeup.
