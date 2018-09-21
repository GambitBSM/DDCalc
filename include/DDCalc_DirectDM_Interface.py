#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Routines for interfacing DirectDM with DDCalc.
Only works with python3.

This requires DirectDM to be installed, e.g. by
   pip3 install directdm
Alternatively, one can specify the location of the directdm folder below

       F. Kahlhoefer   RWTH Aachen   2018
       S. Wild     	  DESY          2018
       ddcalc@projects.hepforge.org
"""

import DDCalcInclude
import sys
import os
directdm_path = DDCalcInclude.ddcalc_dir + 'directdm'
sys.path.append(directdm_path)
import DDCalc
import directdm as ddm

def DM_Spin(DM_type):
	if ((DM_type == "D") or (DM_type == "M")):
		return 0.5
	elif ((DM_type == "C") or (DM_type == "R")):
		return 0.0
	else:
		print("Error in DM_Spin: wrong DM_type.")
		sys.exit()

# This function sets the NR coefficients in DDCalc (WIMPType 'NREFT_CPT') according to
# the dictionary cNR_coeffs, which is the output of DirectDM
def Set_NREFT_CPT_couplings(wimp, cNR_coeffs):
	for key in list(cNR_coeffs.keys()):
		if(not((key[:3] == 'cNR') and (key[-1] in ['p', 'n']) and (int(key[3:-1]) in [i for i in range(1,24)] + [100,104]))):
			print("Set_NREFT_CPT_couplings: Invalid key in cNR_coeffs. Stop.")
			sys.exit()
		OpIndex = int(key[3:-1])
		if key[-1] == 'p':
			tau = 0
		elif key[-1] == 'n':
			tau = 1
		else:
			print("Set_NREFT_CPT_couplings: Invalid key in cNR_coeffs. Stop.")
			sys.exit()
		DDCalc.SetNRCoefficient(wimp, OpIndex, tau, cNR_coeffs[key])
	return


# Define a WIMP by specifying its interactions above the EW scale, at the scale Lambda.
# coupling_dict has to be a dictionary containing the non-zero coefficients
# of the relativistic operators above the EW scale, following 1809.03506.
def SetWIMP_EW(coupling_dict, mDM, DM_type, DM_hypercharge, dim_SU2_repr, Lambda):
	wc_ew = ddm.WC_EW(coupling_dict, DM_hypercharge, dim_SU2_repr, DM_type=DM_type)
	cNR_coeffs = wc_ew._my_cNR(mDM, Lambda)
	wimp = DDCalc.InitWIMP()		# Initialise a WIMP objects to default values.
	DDCalc.SetWIMP_NREFT_CPT(wimp, mDM, DM_Spin(DM_type)) # Set the WIMP to be of type 'NREFT_CPT', with given mass and spin
	Set_NREFT_CPT_couplings(wimp, cNR_coeffs)    # Set all non-relativistic coefficients according to cNR_coeffs
	return wimp


# Define a WIMP by specifying its interactions within the 3-, 4- or 5-flavor scheme.
# coupling_dict has to be a dictionary containing the non-zero coefficients
# of the corresponding relativistic operators, following 1809.03506.
def SetWIMP_Nf(coupling_dict, mDM, DM_type, Nf = 5):
	if(Nf == 5):
		wc = ddm.WC_5f(coupling_dict, DM_type=DM_type)
	elif(Nf == 4):
		wc = ddm.WC_4f(coupling_dict, DM_type=DM_type)
	elif(Nf == 3):
		wc = ddm.WC_3f(coupling_dict, DM_type=DM_type)
	else:
		print("Error in SetWIMP_Nf: Invalid value of Nf. Stop.")
		sys.exit()
	cNR_coeffs = wc._my_cNR(mDM)
	wimp = DDCalc.InitWIMP()		# Initialise a WIMP objects to default values.
	DDCalc.SetWIMP_NREFT_CPT(wimp, mDM, DM_Spin(DM_type)) # Set the WIMP to be of type 'NREFT_CPT', with given mass and spin
	Set_NREFT_CPT_couplings(wimp, cNR_coeffs)    # Set all non-relativistic coefficients according to cNR_coeffs
	return wimp

