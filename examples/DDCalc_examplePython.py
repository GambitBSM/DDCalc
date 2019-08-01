#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
/**********************************************************************
 * DDCALC EXAMPLE PROGRAM (python2 or python3)
 * This program shows how to use the DDCalc module from python, making
 * use of the interface defined in the DDCalc.py wrapper file.
 * 
 * For various types of DM-nucleon interactions, as well as for various 
 * direct detection experiments, it computes the expected signal rates,
 * any by comparing to the observed number of events, also the 
 * log(likelihood) value for each of the examples.
 * 
 * 
 * Run:
 *         python2 DDCalc_examplePython.py 
 *   or    python3 DDCalc_examplePython.py
 * 
 *
 *       F. Kahlhoefer   RWTH Aachen   2018
 *       S. Wild     	  DESY          2018
 *       ddcalc@projects.hepforge.org
 * 
 **********************************************************************/
"""

import DDCalcInclude
import sys
import os
sys.path.append(DDCalcInclude.include_dir)
import DDCalc



# Initialise a DM Halo object to default values.
halo = DDCalc.InitHalo()

# Initialise a WIMP objects to default values.
wimp = DDCalc.InitWIMP()

# Optionally set the Standard Halo Model parameters to values different from the default choices:
#       rho     Local dark matter density [GeV/cm^3]
#       vrot    Local disk rotation speed [km/s]
#       v0      Maxwell-Boltzmann most probable speed [km/s]
#       vesc    Galactic escape speed [km/s]
# DDCalc.SetSHM(halo, 0.3, 235.0, 235.0, 550.0)

# Optionally tell DDCalc to load the velocity integral from file.
# The file must be named halo.txt and must be located in the DDCalc main folder
# The first column of the file contains the values of vmin (in km/s) for which the velocity integral is provided
# The argument Nvmin specifies the number of rows in the file
# The column from which the actual velocity integral should be read can be provided via the argument gcol 
# Finally, the assumed local DM density can be provided as a separate argument
DDCalc.HaloFromFile(halo, 0.3, 1000, 8)



# ****************************************************************************************************************
# Example 1: XENON1T (2017) analysis, with standard SI/SD interactions specified by WIMP-nucleon cross sections.
detector = DDCalc.InitExperiment('Xenon1T_2017')    # Initalize the XENON1T_2017 detector.
mDM = 100.0	 												# DM Mass in GeV                  
sigmap_SI = 4.0e-9											# SI WIMP-proton cross section in pb.
sigman_SI = -0.3e-9											# SI WIMP-neutron cross section in pb.
																#   The negative value indicates that the corresponding WIMP-nucleon coupling is negative.
sigmap_SD = 2.0e-5											# SD WIMP-proton cross section in pb.
sigman_SD = 8.0e-5											# SD WIMP-neutron cross section in pb.

DDCalc.SetWIMP_msigma(wimp, mDM, sigmap_SI, sigman_SI, sigmap_SD, sigman_SD)
DDCalc.CalcRates(detector, wimp, halo) 	             # This performs the actual calculation of the rates.


print("********************************************************************************")
print("Example 1: Mixed SI and SD interactions at XENON1T (2017),")
print("           specified by the WIMP-nucleon cross sections.")
print("mDM       = %.5e GeV" % mDM)
print("sigmap_SI = %.5e pb" % sigmap_SI)
print("sigman_SI = %.5e pb" % sigman_SI)
print("sigmap_SD = %.5e pb" % sigmap_SD)
print("sigman_SD = %.5e pb" % sigman_SD)
print("")
print("Expected number of signal events:      %.5e" % DDCalc.Signal(detector))
print("Observed number of events:             %d" % DDCalc.Events(detector))
print("Expected number of background events:  %.5e" % DDCalc.Background(detector))
print("Log(likelihood):                       %.5e" % DDCalc.LogLikelihood(detector))
print("********************************************************************************")
# ****************************************************************************************************************




# ****************************************************************************************************************
# Example 2: PICO-60 (2017) analysis, with standard SI and SD interactions
#                specified by the couplings fp, fn, ap, an.          
#                fp and fn are the coefficients of the DM DM N N term in the Lagrangian, assuming a Majorana DM particle.
#                The definition of ap and an follow the convention in Jungman&Kamionkowski, 'SUPERSYMMETRIC DARK MATTER'.
detector = DDCalc.InitExperiment('PICO_60_2017')    # Initalize the PICO_60_2017 detector.
mDM = 30.0	 												# DM Mass in GeV                  
fp = 1.0e-8					# SI WIMP-proton coupling fp, in units [GeV^(-2)].
fn = 3.0e-8					# SI WIMP-neutron coupling fn, in units [GeV^(-2)].
ap = 1.0e-2					# SD WIMP-proton coupling ap, unitless.
an = 0.0						# SD WIMP-neutron coupling an, unitless.
	
DDCalc.SetWIMP_mfa(wimp, mDM, fp, fn, ap, an)
DDCalc.CalcRates(detector, wimp, halo) 	             # This performs the actual calculation of the rates.


print("********************************************************************************")
print("Example 2: Mixed SI and SD interactions at PICO-60 (2017),")
print("           specified by the WIMP-nucleon couplings.")
print("mDM       = %.5e GeV" % mDM)
print("fp        = %.5e GeV^(-2)" % fp)
print("fn        = %.5e GeV^(-2)" % fn)
print("ap        = %.5e" % ap)
print("an        = %.5e" % an)
print("")
print("Expected number of signal events:      %.5e" % DDCalc.Signal(detector))
print("Observed number of events:             %d" % DDCalc.Events(detector))
print("Expected number of background events:  %.5e" % DDCalc.Background(detector))
print("Log(likelihood):                       %.5e" % DDCalc.LogLikelihood(detector))
print("********************************************************************************")
# ****************************************************************************************************************





# ****************************************************************************************************************
# Example 3: CRESST (2017) analysis, with a couple of non-relativistic scattering operators 
#            set to a non-zero value.           
#	          In particular, the analysis of this experiment involves several bins in energy.
detector = DDCalc.InitExperiment('CRESST_II')    # Initalize the CRESST_II detector.
mDM = 3.0	 					# DM Mass in GeV                  
DM_spin = 0.5					# DM spin 
	
DDCalc.SetWIMP_NREffectiveTheory(wimp, mDM, DM_spin)
	# This defines a WIMP within the non-relativistic effective theory of DM-nucleon interactions, with a given mass in [GeV] and spin.
	# Initially, all the couplings corresponding to the various operators are set to zero.


# The operator coefficients are set by DDCalc::SetNRCoefficient(WIMP, OpIndex, tau, value).
# OpIndex is an integer specifing the operator, ranging from 1 to 18 (excluding 2 and 16). Additional possible values are -1 and -4,
# corresponding to the momentum suppressed operators O_1 * (q^2/mp^2) and O_4 * (q^2/mp^2), respectively.
# tau is an integer, with 0 corresponding to the isoscalar part of the operator, and 1 to the isovector part.
# value specifies the operator coefficient corresponding to OpIndex and tau, in units [GeV^(-2)].
DDCalc.SetNRCoefficient(wimp, 11, 0, 5.0e-4)	       # This sets the isoscalar part of O_11 to 5.0e-4 GeV^(-2).
DDCalc.SetNRCoefficient(wimp, 11, 1, -2.0e-4)       # This sets the isovector part of O_11 to -2.0e-4 GeV^(-2).
DDCalc.CalcRates(detector,wimp,halo)               #  This performs the actual calculation of the rates.


print("********************************************************************************")
print("Example 3 : Non-standard scattering operators at CRESST-II.")
print("            This is also an example for an experiment with several bins.")
print("mDM       = %.5e GeV" % mDM)
print("DM_spin   = %.5e" % DM_spin)
print("")
print("\t\tExpected signal\t\tObserved\tExpected background")
print("All bins\t%.5e\t\t%d\t\t%.5e" % \
				  (DDCalc.Signal(detector), DDCalc.Events(detector), DDCalc.Background(detector)))
for i_bin in range(1,DDCalc.Bins(detector)+1):
	print("Bin %d\t\t%.5e\t\t%d\t\t%.5e" % \
			   (i_bin, DDCalc.BinSignal(detector, i_bin), DDCalc.BinEvents(detector, i_bin), DDCalc.BinBackground(detector, i_bin)))
print("Log(likelihood):\t\t%.5e" % DDCalc.LogLikelihood(detector))
print("********************************************************************************")
# ****************************************************************************************************************

DDCalc.FreeAll() # Clean up all the objects

