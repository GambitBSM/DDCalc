#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plot exclusion contours from experiments calculated by DDCalc
"""

import DDCalcInclude
import sys
sys.path.append(DDCalcInclude.include_dir)
import DDCalc
import numpy as np
import matplotlib.pyplot as plt
import argparse


# Parse user arguments

parser = argparse.ArgumentParser(description='Plot exclusion limits.')
parser.add_argument('--z_limit', type=float, help='Exclusion limit, e.g. 1.64', default=1.64)

parser.add_argument('--si', dest='si', action='store_true', help='Plot SI scattering cross sections')
parser.add_argument('--no-si', dest='si', action='store_false', help='Do not plot SI scattering cross sections')
parser.set_defaults(si=True)

parser.add_argument('--sd', dest='sd', action='store_true', help='Plot sd scattering cross sections')
parser.add_argument('--no-sd', dest='sd', action='store_false', help='Do not plot sd scattering cross sections')
parser.set_defaults(sd=True)

parser.add_argument('experiments', type=str, default=["PICO_60_2017"], nargs='*', help='List of experiments to plot')

args = parser.parse_args()

# Define the (arbitrary) reference cross section (in pb)
# for which event rates will be calculated
sigma_reference = 1.0e-10

#Initialize WIMP and DM Halo object with default values
Halo = DDCalc.InitHalo()
WIMP = DDCalc.InitWIMP()

# DM masses
mDM = np.logspace(0, 4, 100)

for experiment in args.experiments:

    # Initialize the experiment etc
    Detector = DDCalc.InitExperiment(experiment)
    DDCalc.SetWIMP_msigma(WIMP, mDM[0], 0., 0., 0., 0.)
    DDCalc.CalcRates(Detector, WIMP, Halo)
    BGlogL = DDCalc.LogLikelihood(Detector)

    # For finding contour
    logLlimit = BGlogL - args.z_limit / 2.

    si_limit = []
    sd_limit = []

    for m in mDM:

        DDCalc.SetWIMP_msigma(WIMP, m, 0., 0., sigma_reference, sigma_reference)
        DDCalc.CalcRates(Detector, WIMP, Halo)
        sd_limit.append(DDCalc.ScaleToPValue(Detector, logLlimit) * sigma_reference)

        DDCalc.SetWIMP_msigma(WIMP, m, sigma_reference, sigma_reference, 0., 0.)
        DDCalc.CalcRates(Detector, WIMP, Halo)
        si_limit.append(DDCalc.ScaleToPValue(Detector, logLlimit) * sigma_reference)

    if args.si:
        l = plt.loglog(mDM, si_limit, ls="-", label="SI {}".format(experiment))
        c = l[0].get_color()
        ls = "--"
    else:
        c = None
        ls = None

    if args.sd:
        plt.loglog(mDM, sd_limit, ls=ls, label="SD {}".format(experiment), color=c)

plt.ylim(1e-10, 1e2)
plt.xlabel("$m$ (GeV)")
plt.ylabel(r"$\sigma$ (pb)")
plt.legend(title=r"${}\sigma$ limit with $n = p$".format(args.z_limit))
plt.show()

DDCalc.FreeAll()
