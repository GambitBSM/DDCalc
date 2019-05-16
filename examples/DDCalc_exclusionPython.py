#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
/**********************************************************************
 * DDCALC EXAMPLE PROGRAM (python2 or python3)
 * This program shows how to use the DDCalc module from python, making
 * use of the interface defined in the DDCalc.py wrapper file.
 * 
 * For standard SI/SD scattering, this program calculates upper limits
 * on the scattering cross section for null search experiments,
 * as well as best-fit regions for an experiment seeing an excess.
 * 
 * 
 * Run:
 *         python2 DDCalc_exclusionPython.py 
 *   or    python3 DDCalc_exclusionPython.py
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
import math
import numpy as np


###############################################################################
########################## PREPARATIONS #######################################
###############################################################################

# Define the smallest and largest DM mass to be considered,
# as well as the number of intermediate values.
mDMmin = 1.0
mDMmax = 1.0e4
mDMsteps = 80

sigmin = 1.0e-11
sigmax = 1.0e-7
sigsteps = 80

# Define the (arbitrary) reference cross section (in pb)
# for which event rates will be calculated
sigmatest = 1.0e-10

#Initialize WIMP and DM Halo object with default values
Halo = DDCalc.InitHalo()
WIMP = DDCalc.InitWIMP()	

# Optionally set the Standard Halo Model parameters to values different from the default choices:
#       rho     Local dark matter density [GeV/cm^3]
#       vrot    Local disk rotation speed [km/s]
#       v0      Maxwell-Boltzmann most probable speed [km/s]
#       vesc    Galactic escape speed [km/s]
DDCalc.SetSHM(Halo, 0.3, 220.0, 220.0, 550.0)

# Initialize the Xenon1T_2018 experiment.
Detector = DDCalc.InitExperiment('Xenon1T_2018') 


###############################################################################
#################### Routines for calculating exclusion limits ################
###############################################################################

# Here we assume that the experiment sees no excess, such that the best-fit
# point is identical to the background-only hypothesis.
# See below for the case of experiments seeing an excess.
z
# Step 6: 2d likelihood scan
# To include spectral information, we need to perform a full scan over the
# 2d parameter space. We can then determine the combination of mass and
# cross section that give the best fit to the data and construct
# confidence regions around this point.
logLbest = 1.0*BGlogL
logL = np.zeros((mDMsteps+1, sigsteps+1))
for mDMi in range(mDMsteps+1):
	mDM = mDMmin * (mDMmax/mDMmin)**(mDMi/(1.0*mDMsteps))
	for sigi in range(sigsteps+1):
		sig = sigmin * (sigmax/sigmin)**(sigi/(1.0*sigsteps))
		DDCalc.SetWIMP_msigma(WIMP, mDM, sig, sig, 0., 0.)
		DDCalc.CalcRates(Detector,WIMP,Halo)
		logL[mDMi,sigi] = DDCalc.LogLikelihood(Detector)
		if(logL[mDMi,sigi] > logLbest):
			logLbest = 1.0*logL[mDMi,sigi]
			mDMbest = 1.0*mDM
			sigbest = 1.0*sig
	 
print("")
print("**************************** 2d likelihood scan ****************************")
print("")
print("Assuming spin-independent scattering with equal couplings to protons and neutrons.")
print("")
print("Best fit point is m_DM = %7.1f GeV, sigma_p = %10.3e pb with log L = %7.3f" % \
			  (mDM,sigbest,logLbest))
print("")
print("    log_10(m_DM/GeV)  log_10(sigma_p/pb)               log L      -2 Delta log L")
print("")
for mDMi in range(mDMsteps+1):
	mDM = mDMmin * (mDMmax/mDMmin)**(mDMi/(1.0*mDMsteps))
	for sigi in range(sigsteps+1):
		sig = sigmin * (sigmax/sigmin)**(sigi/(1.0*sigsteps))
		print("          %10.4f          %10.4f          %10.4f          %10.4f" % \
			(math.log(mDM)/math.log(10.),\
			math.log(sig)/math.log(10.),\
			logL[mDMi,sigi],\
			-2.0*(logL[mDMi,sigi]-logLbest)))
    
DDCalc.FreeAll() # Clean up all the objects
