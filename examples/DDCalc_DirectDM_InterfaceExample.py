#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
/**********************************************************************
 * EXAMPLE PROGRAM FOR INTERFACING DirectDM WITH DDCalc 
 * (python3)
 * This program shows how to interface DirectDM with DDCalc.
 * Only python3 is supported for this interface.
 *
 * This requires DirectDM to be installed, e.g. by
 *    pip3 install directdm
 * Alternatively, one can manually put all the directdm source files in a
 * directory called 'directdm' within the main DDCalc directory
 * 
 * Starting from a relativistic Lagrangian, DirectDM performs the running
 * of the Wilson coefficients down to a scale where they can be matched
 * to non-relativistic operators. The Wilson coefficients of these
 * non-relativistic operators are then fed into DDCalc, which calculates
 * event rates, the likelihood, etc. for experiments.
 * 
 * If you use the DirectDM-DDCalc interface, please cite both
 * the relevant DDCalc papers (in particular 1705.07920 and 1808.10465), 
 * as well as the relevant DirectDM papers (in particular 1708.02678
 * , 1710.10218 and 1809.03506)
 * 
 * Run:    python3 DDCalc_DirectDM_InterfaceExample.py
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
import DDCalc_DirectDM_Interface




halo = DDCalc.InitHalo()	   # Initialise a DM Halo object to default values.

# Optionally set the Standard Halo Model parameters to values different from the default choices:
#       rho     Local dark matter density [GeV/cm^3]
#       vrot    Local disk rotation speed [km/s]
#       v0      Maxwell-Boltzmann most probable speed [km/s]
#       vesc    Galactic escape speed [km/s]
DDCalc.SetSHM(halo, 0.3, 235.0, 235.0, 550.0)
detector = DDCalc.InitExperiment('Xenon1T_2018')    # Initalize the XENON1T_2018 detector.

###############################################################################
### Example 1 : Specifying the Lagrangian above the electroweak scale, 
###             i.e. in terms of its couplings to SM quarks, gauge bosons,
###	             the Higgs boson etc.

mDM = 50.0	 		# DM Mass in GeV 
DM_type = "D"       # DM type. Only Dirac DM ("D") is currently supported
						# by DirectDM for the case of specifying the Lagrangian
						# above the EW scale.
DM_hypercharge = 0 # DM hypercharge
dim_SU2_repr = 1	   # dimension of the SU(2) representation of the WIMP
Lambda = 100.0     # Scale in GeV at which the Lagrangian is specified

# Definition of the non-zero couplings. This follows the convention of 
# relativistic operators given in 1809.03506 and the manual of DirectDM.
# In this example, only the operator Q^(5)_3 is set to a non-zero value,
# corresponding to a coupling of the dark matter particle to the
# SM Higgs doublet via a dimension five operator.
coupling_dict = {'C53' : 5.0e-5}

# This function defines a WIMP object with the properties given as arguments,
# runs the couplings down to the hadronic scale, matches onto non-relativistic
# operators, and defines a corresponding WIMP in DDCalc.
# It returns the WIMP index to be used later on for rate calculations.
wimp = DDCalc_DirectDM_Interface.SetWIMP_EW(coupling_dict, mDM, \
								DM_type, DM_hypercharge, dim_SU2_repr, Lambda)

# This performs the actual calculation of the rates.
DDCalc.CalcRates(detector, wimp, halo)

print("**********************************************************************")
print("Example 1: Specify the Lagrangian above the EW scale:")
print("           Singlet Dirac DM particle coupled to the SM Higgs doublet")
print("           via an effective dimension-five operator.")
print("           mDM =", mDM, "GeV")
print("           Experiment: Xenon1T_2018")
print("Expected number of signal events:      %.5e" % DDCalc.Signal(detector))
print("Observed number of events:             %d" % DDCalc.Events(detector))
print("Expected number of background events:  %.5e" \
				  % DDCalc.Background(detector))
print("Log(likelihood):                       %.5e" \
				  % DDCalc.LogLikelihood(detector))
print("**********************************************************************")
###############################################################################






###############################################################################
### Example 2 : Specifying the Lagrangian below the electroweak scale, 
###             i.e. in terms of its couplings to the photon, gluon,
###	             and to a given number of quarks

mDM = 20.0	 		# DM Mass in GeV 
DM_type = "C"       # DM type. Can be "D" for Dirac, "M" for Majorana,
						# "R" for a real scalar and "C" for a complex scalar
Nf = 5					# Flavor scheme, specifying which quark flavors
						# are already considered to be integrated out.
						# This can be 3,4 or 5.
						# See 1708.02678 and 1809.03506 for more details.

# Definition of the non-zero couplings. This follows the convention of 
# relativistic operators given in 1809.03506 and the manual of DirectDM.
# In this example, we use the five-flavor scheme, and consider the operator
# Q^(6)_4, corresponding to a CP-violating coupling of a complex scalar 
# to the pseudoscalar quark current, together with Q^(6)_6 for the
# effective coupling to gluons.
coupling_value = 3.0e-5
coupling_dict = {'C64u' : coupling_value, 'C64d' : coupling_value, \
				    'C64s' : coupling_value, 'C64c' : coupling_value, \
					 'C64b' : coupling_value, 'C66' : 0.5*coupling_value}

# This function defines a WIMP object with the properties given as arguments,
# runs the couplings down to the hadronic scale, matches onto non-relativistic
# operators, and defines a corresponding WIMP in DDCalc.
# It returns the WIMP index to be used later on for rate calculations.
wimp = DDCalc_DirectDM_Interface.SetWIMP_Nf(coupling_dict, \
									mDM, DM_type, Nf = Nf)

# This performs the actual calculation of the rates.
DDCalc.CalcRates(detector, wimp, halo)

print("**********************************************************************")
print("Example 2: Specify the Lagrangian within the five-flavor scheme:")
print("           Complex scalar with CP-violating couplings to quarks")
print("           and gluons.")
print("           mDM =", mDM, "GeV")
print("           Experiment: Xenon1T_2018")
print("Expected number of signal events:      %.5e" % DDCalc.Signal(detector))
print("Observed number of events:             %d" % DDCalc.Events(detector))
print("Expected number of background events:  %.5e" \
				  % DDCalc.Background(detector))
print("Log(likelihood):                       %.5e" \
				  % DDCalc.LogLikelihood(detector))
print("**********************************************************************")
###############################################################################




DDCalc.FreeAll() # Clean up all the objects
