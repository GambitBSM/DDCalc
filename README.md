DDCalc
======

The DDCalc software package is a set of routines and a frontend for
doing various dark matter direct detection calculations for several
experimental results.  Most notably, these include Poisson likelihoods
(binned and unbinned) and p-values based on the maximum gap method (Yellin 2002).

A full description of this package and the physics framework behind
it can be found in the [GAMBIT DarkBit paper](DarkBitPaper.pdf):

T Bringmann, J Conrad, JM Cornell, LA Dal, J Edsjö, B Farmer,
F Kahlhoefer, A Kvellestad, A Putze, C Savage, P Scott, C Weniger,
M White & S Wild, DarkBit: A GAMBIT module for computing dark matter
observables and likelihoods, EPJC submitted (2017), arXiv:1705.07920.

The new features in DDCalc_v2 (in particular momentum and velocity
dependent interactions and binned likelihoods) are described in:

XXX

If you write a paper that uses DDCalc, please cite both papers. 


Experiments
--

The experimental results implemented in DDCalc at this time are:

  * Xenon1T 2017: 1042*34.2 kg-days (Xe)  
    E. Aprile et al., [arxiv:1705.06655]

  * LUX 2016: 33500 kg-days (Xe)  
    D.S. Akerib et al., PRL 118, 021303 (2017) [arxiv:1608.07648]

  * PandaX 2016: 334.3*98.7 kg-days (Xe)  
    A. Tan et al., PRL 117, 121303 (2016) [arxiv:1607.07400]

  * PandaX 2017: 361.5*77.1 kg-days (Xe)  
    X. Cui et al., PRL 119, 181302 (2017) [arxiv:1708.06917]

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

  * CDMSlite: 70.1 kg-days (Ge)  
    R. Agnese et al., PRL 116, 071301 (2016) [arxiv:1509.02448]

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

The new features in DDCalc_v2 are described in:

XXX

If you use the Yellin maximum gap method, please cite:

  S. Yellin, PRD 66, 032005 (2002) [physics/0203002]  
    http://inspirehep.net/record/583791  

If using DDCalc for calculations involving experimental results,
please cite the appropriate experimental papers listed above.


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
    the DDCalc routines externally from either Fortran or C++ to 
    calculate observables and likelihoods.

  * DDCalc_exclusionF, DDCalc_exclusionC
    Simple example programs with source code showing how to use
    the DDCalc routines externally from either Fortran or C++ to
    calculate exclusion limits.

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


Contact
-------

Any questions, comments, or suggestions should be directed to the GAMBIT
Dark Matter Workgroup, via ddcalc@projects.hepforge.org.


Contributors
------------

The following have contributed to the development and
testing of the DDCalc package:

  * Chris Savage  
    Most fortran procedures.

  * Felix Kahlhoefer, RWTH Aachen  
    DDCalc_v2.

  * Sebastian Wild, DESY  
    DDCalc_v2.

  * Pat Scott, Imperial College London  
    C++ interface and example, testing, fortran structure, bug fixes, writeup.

  * Martin White, Univ of Adelaide  
    C++ interface and example, testing, writeup.

  * Andre Scaffidi, Univ of Adelaide  
    SD form factors, C++ interface and example, testing.

  * Lauren Hsu, FermiLab  
    SuperCDMS.

  * Miguel Pato
    Calibration/tuning of LUX parameters.

