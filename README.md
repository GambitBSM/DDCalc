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
observables and likelihoods, Eur.Phys.J. C77 (2017) no.12, 831, arXiv:1705.07920.

The new features in DDCalc_v2 (in particular momentum and velocity
dependent interactions and binned likelihoods) are described in the 
[Higgs Portal paper](HiggsPortalPaper.pdf)

Peter Athron, Csaba Balázs, Ankit Beniwal, Sanjay Bloor, José Eliel Camargo-Molina, 
Jonathan M. Cornell, Ben Farmer, Andrew Fowlie, Tomás E. Gonzalo, Felix Kahlhoefer, 
Anders Kvellestad, Gregory D. Martinez, Pat Scott, Aaron C. Vincent, Sebastian Wild, 
Martin White, Anthony G. Williams, Global analyses of Higgs portal singlet dark 
matter models using GAMBIT, EPJC submitted (2018), arXiv:1808.10465.

If you write a paper that uses DDCalc, please cite both papers. 


Experiments
--

The experimental results implemented in DDCalc at this time are (in alphabetical order):

  * CDMSlite: 70.1 kg-days (Ge)  
    R. Agnese et al., PRL 116, 071301 (2016) [arxiv:1509.02448]

  * CRESST-II: 52 kg-days (CaWO4)  
    G. Angloher et al., Eur.Phys.J. C76, 25 (2016) [arXiv:1509.01515]

  * CRESST-III: 3.64 kg-days (CaWO4)  
    G. Angloher et al., PRD 100, 102002 (2019) [arXiv:1904.00498]

  * DarkSide-50: S2-only (Ar)  
    P. Agnes et al., PRL 121, 081307 (2018) [arXiv:1802.06994]

  * DarkSide-50: 19600 kg-days (Ar)  
    P. Agnes et al., PRD 98, 102006 (2018) [arXiv:1802.07198]

  * LUX 2013: 118*85.3 kg-days (Xe)  
    D.S. Akerib et al., PRL 112, 091303 (2014) [arxiv:1310.8214]

  * LUX 2015: 118*85.3 kg-days (Xe)  
    D.S. Akerib et al., PRL 116, 161301 (2016) [arXiv:1512.03506]

  * LUX 2016: 33500 kg-days (Xe)  
    D.S. Akerib et al., PRL 118, 021303 (2017) [arxiv:1608.07648]

  * PandaX 2016: 334.3*98.7 kg-days (Xe)  
    A. Tan et al., PRL 117, 121303 (2016) [arxiv:1607.07400]

  * PandaX 2017: 361.5*77.1 kg-days (Xe)  
    X. Cui et al., PRL 119, 181302 (2017) [arxiv:1708.06917]

  * PICO 2L: 104 kg-days (C3F8)
    C. Amole et al., PRD 93, 061101 (2016) [arXiv:1601.03729] 

  * PICO 60: 741 kg-days (CF3)
    C. Amole et al., PRD 93, 052014 (2016) [arXiv:1510.07754]  

  * PICO 60 2017: (C3F8)
    C. Amole et al., PRL 118, 251301 (2017) [arXiv:1702.07666]

  * PICO 60 2019: (C3F8)
    C. Amole et al., PRD 100, 022001 (2019) [arXiv:1902.04031]

  * SIMPLE 2014: 18.24 kg-days (C2ClF5)  
    M. Felizardo et al., PRD 89, 072013 (2014) [arXiv:1404.4309]

  * SuperCDMS 2014: 577 kg-days (Ge)  
    R. Agnese et al., PRL 112, 241302 (2014) [arxiv:1402.7137]

  * XENON100 2012: 34*224.6 kg-days (Xe)  
    E. Aprile et al., PRL 109, 181301 (2012) [arxiv:1207.5988]

  * XENON1T 2017: 1042*34.2 kg-days (Xe)  
    E. Aprile et al., PRL 119 181301 (2017) [arxiv:1705.06655]

  * XENON1T 2018: 1300*279 kg-days (Xe)  
    E. Aprile et al., PRL 121, 111302 (2018) [arxiv:1805.12562]

Note that the following data sets are not statistically independent and should therefore not be combined:
  * LUX 2013, LUX 2015 and LUX 2016
  * XENON1T 2017 and XENON1T 2018
In contrast, the following data sets are independent and can be combined:
  * PICO 60, PICO 60 2017, PICO 60 2019
  * PandaX 2016 and PandaX 2017 (these correspond to Run 9 and Run 10, respectively, as described in arXiv:1708.06917)

Moreover, projections for a number of future facilities have been implemented.

Near future:

* LZ (Xe)
  D.S. Akerib et al., arXiv:1509.02910

* PICO-500 (F)
  S. Fallows, https://indico.cern.ch/event/606690/contributions/2623446/attachments/1497228/2330240/Fallows_2017_07_24__TAUP__PICO-40L_v1.2.pdf

Far future:

* DARWIN (Xe)
  J. Aalbers et al., JCAP 1611 (2016) 017 [arXiv:1606.07001]
  M. Schumann et al., JCAP 1510 (2015) 016 [arXiv:1506.08309]

* DarkSide-20k (Ar)
  C.E. Aalseth et al., Eur.Phys.J.Plus 133 (2018) 131 [arXiv:1707.08145]

Citation Policy
--

When using DDCalc, the following publications must be cited; they can be
found in the 'DDCalc.bib' file.

The DDCalc software itself is presented in the GAMBIT DarkBit paper:

T Bringmann, J Conrad, JM Cornell, LA Dal, J Edsjö, B Farmer,
F Kahlhoefer, A Kvellestad, A Putze, C Savage, P Scott, C Weniger,
M White & S Wild, DarkBit: A GAMBIT module for computing dark matter
observables and likelihoods, Eur.Phys.J. C77 (2017) no.12, 831, arXiv:1705.07920.

The new features in DDCalc_v2 are described in:

Peter Athron, Csaba Balázs, Ankit Beniwal, Sanjay Bloor, José Eliel Camargo-Molina, 
Jonathan M. Cornell, Ben Farmer, Andrew Fowlie, Tomás E. Gonzalo, Felix Kahlhoefer, 
Anders Kvellestad, Gregory D. Martinez, Pat Scott, Aaron C. Vincent, Sebastian Wild, 
Martin White, Anthony G. Williams, Global analyses of Higgs portal singlet dark 
matter models using GAMBIT, EPJC submitted (2018), arXiv:1808.10465.

If you use the Yellin maximum gap method, please cite:

  S. Yellin, PRD 66, 032005 (2002) [physics/0203002]  

If you use the tabulated velocity integrals, please cite:

  M. Kuhlen et al., JCAP 1002 (2010) 030 [arXiv:0912.2358]

If you use the interface with DirectDM, please cite:

  Bishara, Brod, Grinstein, and Zupan, JCAP 1702 (2017) 009 [arXiv:1611.00368]
  Bishara, Brod, Grinstein, and Zupan, JHEP 1711 (2017) 059 [arXiv:1707.06998]
  Bishara, Brod, Grinstein, and Zupan [arXiv:1708.02678]
  Brod, Grinstein, and Zupan [arXiv:1710.10218]
  Brod, Gootjes-Dreesbach, Tammaro, and Zupan, JHEP 1802 (2018) 174 [arXiv:1801.04240]
  Bishara, Brod, Grinstein, and Zupan [arXiv:1809.03506]

If using DDCalc for calculations involving experimental results,
please cite the appropriate experimental papers listed above.


Compiling
---------

DDCalc is written as a Fortran 95 module.  It supports both the  Intel
and GNU compilers.

To build, simply do

  make all

This will generate the following:

  * DDCalc_test
    Tests whether DDCalc has been correctly compiled. This program simply serves as a placeholder for more complex main programs.

  * DDCalc.o, DDCalc.a, DDCalc.so
    Object, static library, and shared object files for linking
    with other software.

  * DDCalc_exampleF, DDCalc_exampleC, DDCalc_examplePython.py
    Simple example programs with source code showing how to use
    the DDCalc routines externally from either Fortran, C++ or
    Python to calculate observables and likelihoods.

  * DDCalc_exclusionF, DDCalc_exclusionC, DDCalc_exclusionPython.py
    Simple example programs with source code showing how to use
    the DDCalc routines externally from either Fortran, C++ or
    Python to calculate exclusion limits.

  * DDCalc_DirectDM_InterfaceExample.py
    Simple example programs with source code showing how to use
    the DDCalc-DirectDM interface in Python to calculate observables
    based on effective operators defined at some high scale.

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

  * Gonzalo Herrera
    CRESST-III.
