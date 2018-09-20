#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 20 12:56:12 2018

Python interface for DDCalc.
This works by wrapping the functions provided by the C interface for DDCalc
into appropriate python functions.

Created by Sebastian Wild, Sep 2018
"""

DDCalc_library_path = 'lib/libDDCalc.so'

import ctypes
import sys


##### Load the DDCalc library #################################################
try:
	ddcalc_lib = ctypes.CDLL(DDCalc_library_path) 
	print("Successfully loaded DDCalc library " + DDCalc_library_path)
except OSError:
	print("DDCalc library " + DDCalc_library_path + " does not exist. Stop.")
	sys.exit()
###############################################################################




##### Auxiliary functions #####################################################	
int_byref = lambda x : ctypes.byref(ctypes.c_int(x))
double_byref = lambda x : ctypes.byref(ctypes.c_double(x))
###############################################################################



##### Initializing experiments ################################################
def InitExperiment(ExperimentName):
	experiments = {
'Xenon100_2012' : ddcalc_lib.C_DDCalc_xenon100_2012_init(),
'LUX_2013' : ddcalc_lib.C_DDCalc_lux_2013_init(),
'LUX_2016' : ddcalc_lib.C_DDCalc_lux_2016_init(),
'PandaX_2016' : ddcalc_lib.C_DDCalc_pandax_2016_init(),
'PandaX_2017' : ddcalc_lib.C_DDCalc_pandax_2017_init(),
'Xenon1T_2017' : ddcalc_lib.C_DDCalc_xenon1t_2017_init(),
'Xenon1T_2018' : ddcalc_lib.C_DDCalc_xenon1t_2018_init(),
'LUX_2015' : ddcalc_lib.C_DDCalc_lux_2015_init(),
'PICO_2L' : ddcalc_lib.C_DDCalc_pico_2l_init(),
'PIC0_60' : ddcalc_lib.C_DDCalc_pico_60_init(),
'PICO_60_2017' : ddcalc_lib.C_DDCalc_pico_60_2017_init(),
'SuperCDMS_2014' : ddcalc_lib.C_DDCalc_supercdms_2014_init(),
'CDMSlite' : ddcalc_lib.C_DDCalc_cdmslite_init(),
'Simple_2014' : ddcalc_lib.C_DDCalc_simple_2014_init(),
'CRESST_II' : ddcalc_lib.C_DDCalc_cresst_ii_init(),
'LZ' : ddcalc_lib.C_DDCalc_lz_init(),
'Darwin' : ddcalc_lib.C_DDCalc_darwin_init(),
'Darkside_20k' : ddcalc_lib.C_DDCalc_darkside_20k_init(),
'Darkside_50' : ddcalc_lib.C_DDCalc_darkside_50_init(),
'PICO_500' : ddcalc_lib.C_DDCalc_pico_500_init()
			}
	if not ExperimentName in list(experiments.keys()):
		print("InitExperiment: Invalid experiment name " + ExperimentName + \
			". Available experiments are:")
		print(list(experiments.keys()))
		sys.exit()
	return experiments[ExperimentName]
###############################################################################


##### Initializing WIMP/Halo ##################################################
def InitWIMP():
	return ddcalc_lib.C_DDWIMP_ddcalc_initwimp()

def InitHalo():
	return ddcalc_lib.C_DDHalo_ddcalc_inithalo()
###############################################################################




##### Setting WIMP properties #################################################
def SetWIMP_mfa(wimp,m,fp,fn,ap,an):
	ddcalc_lib.C_DDCalc_ddcalc_setwimp_mfa(int_byref(wimp),\
								double_byref(m),double_byref(fp),\
								double_byref(fn),double_byref(ap),\
								double_byref(an))

def SetWIMP_mG(wimp,m,GpSI,GnSI,GpSD,GnSD):
	ddcalc_lib.C_DDCalc_ddcalc_setwimp_mg(int_byref(wimp),\
								double_byref(m),double_byref(GpSI),\
								double_byref(GnSI),double_byref(GpSD),\
								double_byref(GnSD))

def SetWIMP_msigma(wimp, mDM, sigmap_SI, sigman_SI, sigmap_SD, sigman_SD):
	ddcalc_lib.C_DDCalc_ddcalc_setwimp_msigma(int_byref(wimp),\
				   double_byref(mDM),double_byref(sigmap_SI),\
				   double_byref(sigman_SI),double_byref(sigmap_SD),\
				   double_byref(sigman_SD))
###############################################################################

def SetSHM(halo, rho, vrot, v0, vesc):
	ddcalc_lib.C_DDCalc_ddcalc_setshm(int_byref(halo), \
				   double_byref(rho), double_byref(vrot), \
				   double_byref(v0), double_byref(vesc))
	
def SetWIMP_msigma(wimp, mDM, sigmap_SI, sigman_SI, sigmap_SD, sigman_SD):
	ddcalc_lib.C_DDCalc_ddcalc_setwimp_msigma(int_byref(wimp),\
				   double_byref(mDM),double_byref(sigmap_SI),\
				   double_byref(sigman_SI),double_byref(sigmap_SD),\
				   double_byref(sigman_SD))
	
def CalcRates(detector, wimp, halo):
	ddcalc_lib.C_DDRates_ddcalc_calcrates(int_byref(detector),\
								   int_byref(wimp), int_byref(halo))
	
def Signal(detector):
	ddcalc_lib.C_DDRates_ddcalc_signal.restype = ctypes.c_double
	return ddcalc_lib.C_DDRates_ddcalc_signal(int_byref(detector))
	


