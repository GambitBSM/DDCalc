#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 20 12:56:12 2018

Python interface for DDCalc.
This works by wrapping the functions provided by the C interface for DDCalc
into appropriate python functions.
See DDCalc.hpp for more explanations regarding the meaning of each function.

Created by Sebastian Wild, Sep 2018
"""

import DDCalcInclude
import os
import ctypes
import sys

DDCalc_library = DDCalcInclude.lib_dir + '/libDDCalc.so'

##### Load the DDCalc library #################################################
try:
	ddcalc_lib = ctypes.CDLL(DDCalc_library) 
	#print("Successfully loaded DDCalc library " + DDCalc_library)
except OSError:
	print("DDCalc library " + DDCalc_library + " does not exist. Stop.")
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

def SetWIMP_msigma(wimp, m, sigmap_SI, sigman_SI, sigmap_SD, sigman_SD):
	ddcalc_lib.C_DDCalc_ddcalc_setwimp_msigma(int_byref(wimp),\
				   double_byref(m),double_byref(sigmap_SI),\
				   double_byref(sigman_SI),double_byref(sigmap_SD),\
				   double_byref(sigman_SD))

def SetWIMP_longrange(wimp, m, gp, gn, mmed):
	ddcalc_lib.C_DDCalc_ddcalc_setwimp_longrange(int_byref(wimp),\
				   double_byref(m),double_byref(gp),\
				   double_byref(gn),double_byref(mmed))
	
def SetWIMP_NREffectiveTheory(wimp, m, spin):
	ddcalc_lib.C_DDCalc_ddcalc_setwimp_nreffectivetheory(int_byref(wimp), \
					double_byref(m), double_byref(spin))
	
def SetWIMP_NREFT_CPT(wimp, m, spin):
	ddcalc_lib.C_DDCalc_ddcalc_setwimp_nreft_cpt(int_byref(wimp),\
					  double_byref(m), double_byref(spin))
	
def SetNRCoefficient(wimp, OpIndex, tau, value):
	ddcalc_lib.C_DDCalc_ddcalc_setnrcoefficient(int_byref(wimp), \
					 int_byref(OpIndex), int_byref(tau), double_byref(value))
###############################################################################



##### Set Halo properties #####################################################
def SetSHM(halo, rho, vrot, v0, vesc):
	ddcalc_lib.C_DDCalc_ddcalc_setshm(int_byref(halo), \
				   double_byref(rho), double_byref(vrot), \
				   double_byref(v0), double_byref(vesc))
###############################################################################


##### Set or retrieve detector properties #####################################
def SetDetectorEmin(detector, Emin):
	ddcalc_lib.C_DDCalc_ddcalc_setdetectoremin(int_byref(detector),\
					double_byref(Emin))

def Bins(detector):
	return ddcalc_lib.C_DDRates_ddcalc_bins(int_byref(detector))

def BinEvents(detector, BinIndex):
	return ddcalc_lib.C_DDRates_ddcalc_binevents(int_byref(detector), \
						  int_byref(BinIndex))
	
def BinBackground(detector, BinIndex):
	ddcalc_lib.C_DDRates_ddcalc_binbackground.restype = ctypes.c_double
	return ddcalc_lib.C_DDRates_ddcalc_binbackground(int_byref(detector), \
						  int_byref(BinIndex))
	
def BinSignal(detector, BinIndex):
	ddcalc_lib.C_DDRates_ddcalc_binsignal.restype = ctypes.c_double
	return ddcalc_lib.C_DDRates_ddcalc_binsignal(int_byref(detector), \
						  int_byref(BinIndex))

	
###############################################################################
	
	
	
##### Calculate rates #########################################################	
def CalcRates(detector, wimp, halo):
	ddcalc_lib.C_DDRates_ddcalc_calcrates(int_byref(detector),\
								   int_byref(wimp), int_byref(halo))
###############################################################################	


##### Inspect results of calculations #########################################
def Events(detector):
	return ddcalc_lib.C_DDRates_ddcalc_events(int_byref(detector))		
	
def Background(detector):
	ddcalc_lib.C_DDRates_ddcalc_background.restype = ctypes.c_double
	return ddcalc_lib.C_DDRates_ddcalc_background(int_byref(detector))	
	
def Signal(detector):
	ddcalc_lib.C_DDRates_ddcalc_signal.restype = ctypes.c_double
	return ddcalc_lib.C_DDRates_ddcalc_signal(int_byref(detector))

def LogLikelihood(detector):
	ddcalc_lib.C_DDStats_ddcalc_loglikelihood.restype = ctypes.c_double
	return ddcalc_lib.C_DDStats_ddcalc_loglikelihood(int_byref(detector))

def ScaleToPValue(detector, logp=-2.302585):
	ddcalc_lib.C_DDStats_ddcalc_scaletopvalue.restype = ctypes.c_double
	return ddcalc_lib.C_DDStats_ddcalc_scaletopvalue(int_byref(detector), \
					  double_byref(logp))
	
def FeldmanCousinsUpper(LnP, N, B):
	ddcalc_lib.C_DDStats_ddcalc_feldmancousinsupper.restype = ctypes.c_double
	return ddcalc_lib.C_DDStats_ddcalc_feldmancousinsupper(double_byref(LnP), \
				int_byref(N), double_byref(B))
	
def FeldmanCousinsLower(LnP, N, B):
	ddcalc_lib.C_DDStats_ddcalc_feldmancousinslower.restype = ctypes.c_double
	return ddcalc_lib.C_DDStats_ddcalc_feldmancousinslower(double_byref(LnP), \
				int_byref(N), double_byref(B))
###############################################################################	



##### Memory cleanup ##########################################################
def FreeWIMPs():
	ddcalc_lib.C_DDUtils_ddcalc_freewimps()
	
def FreeHalos():
	ddcalc_lib.C_DDUtils_ddcalc_freehalos()
	
def FreeDetectors():
	ddcalc_lib.C_DDUtils_ddcalc_freedetectors()
	
def FreeAll():
	ddcalc_lib.C_DDUtils_ddcalc_freeall()
###############################################################################	
