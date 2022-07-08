import sys
import os
from numpy import sqrt
import numpy as np
from scipy.special import erf
from scipy.interpolate import interp1d
import math
import matplotlib.pyplot as plt

import DDCalcInclude
sys.path.append(DDCalcInclude.include_dir)
import DDCalc

# Define the smallest and largest DM mass to be considered,
# as well as the number of intermediate values.
mDMmin = 1
mDMmax = 10000
mDMsteps = 80

sigmin = 1.0e-12
sigmax = 1.0e-6
sigsteps = 80

# Define the (arbitrary) reference cross section (in pb)
# for which event rates will be calculated
sigmatest = 1.0

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
Detector = DDCalc.InitExperiment('LZ_2022') 

DDCalc.SetWIMP_msigma(WIMP, mDMmin, 0., 0., 0., 0.)
DDCalc.CalcRates(Detector,WIMP,Halo)
BGlogL = DDCalc.LogLikelihood(Detector)

limit = 4.6
logLlimit = BGlogL - limit/2.0
print("logLlimit = %.5e" % logLlimit)

x = np.zeros(mDMsteps+1)
y = np.zeros(mDMsteps+1)
for mDMi in range(mDMsteps+1):
	mDM = mDMmin * (mDMmax/mDMmin)**(mDMi/(1.0*mDMsteps))
	DDCalc.SetWIMP_msigma(WIMP, mDM, sigmatest, sigmatest, 0., 0.)
	DDCalc.CalcRates(Detector,WIMP,Halo)
	print("          %.5f          %.5f" % \
	   (math.log(mDM)/math.log(10.0),\
	   math.log(DDCalc.ScaleToPValue(Detector,logLlimit)*sigmatest)/math.log(10.0)))
	x[mDMi]=mDM
	y[mDMi]=DDCalc.ScaleToPValue(Detector,logLlimit)*sigmatest

data = np.loadtxt('LZ_limit.dat')

print(data[:,0])
print(data[:,1]*1e36)
print(x,y)
plt.xlim([mDMmin, mDMmax])
plt.ylim([sigmin,sigmax])
plt.xscale('log')
plt.yscale('log')
plt.xlabel(r'$m_\mathrm{DM} [\mathrm{GeV}]$')
plt.ylabel(r'$\sigma_p [\mathrm{pb}]$')
plt.tight_layout()
plt.plot(x,y)
plt.plot(data[:,0],data[:,1]*1e36)
plt.savefig('LZ.pdf')  
plt.show()

exit()
