#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 20 12:04:29 2018

@author: sebwild
"""

import DDCalc

wimp = DDCalc.InitWIMP()
halo = DDCalc.InitHalo()

detector = DDCalc.InitExperiment('Xenon1T_2017')



mDM = 100.0;                   
sigmap_SI = 4.0e-9
sigman_SI = -0.3e-9
sigmap_SD = 2.0e-5
sigman_SD = 8.0e-5
	

DDCalc.SetSHM(halo, 0.3, 235.0, 235.0, 550.0)
DDCalc.SetWIMP_msigma(wimp, mDM, sigmap_SI, sigman_SI, sigmap_SD, sigman_SD)
DDCalc.CalcRates(detector, wimp, halo)
print(DDCalc.Signal(detector))
