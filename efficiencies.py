from __future__ import division

from numpy import sqrt
import numpy as np
from scipy.special import erf
from scipy.interpolate import interp1d

GeneralAcceptanceFunction = lambda ER, E_lower, E_upper, sigmaE : 0.5 * (erf((E_upper - ER)/(np.sqrt(2)*sigmaE)) - erf((E_lower - ER)/(np.sqrt(2)*sigmaE)))
	# acceptance function to have a detected recoil within [E_lower, E_upper], for a true recoil ER, and a constant energy resolution sigmaE
	# this does not take into account experimental-specific additional efficiencies

# LZ

bins = [i for i in range(6,31,4)]

print bins

energies = [0.2*j + 2 for j in range(241)]

sigmaE = lambda ER: 0.065 * ER + 0.24 * np.sqrt(ER)

efficiencies = np.array([[sum([0.5 * GeneralAcceptanceFunction(energies[j], bins[i], bins[i+1], sigmaE(energies[j]))for i in range(len(bins) - 1)])] + [0.5 * GeneralAcceptanceFunction(energies[j], bins[i], bins[i+1], sigmaE(energies[j]))for i in range(len(bins) - 1)] for j in range(len(energies))])

np.savetxt('data/LZ/efficiencies.dat',efficiencies)

np.savetxt('data/LZ/energies.dat',np.array(energies))

# DarkSide

acceptance = np.loadtxt('data/DarkSide/DarkSide.dat')

acceptancefunc = interp1d(acceptance[:,0], acceptance[:,1])

print acceptancefunc(40)

bins = [i for i in range(30,81,5)]

energies = [0.4*j + 20 for j in range(201)]

sigmaE = lambda ER: 0.05 * ER + 2

efficiencies = np.array([ [sum([acceptancefunc(energies[j]) * GeneralAcceptanceFunction(energies[j], bins[i], bins[i+1], sigmaE(energies[j]))for i in range(len(bins) - 1)])] + [acceptancefunc(energies[j]) * GeneralAcceptanceFunction(energies[j], bins[i], bins[i+1], sigmaE(energies[j]))for i in range(len(bins) - 1)] for j in range(len(energies))])

np.savetxt('data/DarkSide/efficiencies.dat',efficiencies)

np.savetxt('data/DarkSide/energies.dat',np.array(energies))

# PICO

bins = [3.2, 10, 50]

energies = [0.1*j + 3.2 for j in range(68)] + [0.4*j + 10 for j in range(101)] 

efficiency = [[0.5 for j in range(68)] + [1 for j in range(101)], [0.5 for j in range(68)] + [0 for j in range(101)], [0 for j in range(68)] + [1 for j in range(101)]]

np.savetxt('data/PICO-500/efficiencies.dat',np.transpose(np.array(efficiency)))

np.savetxt('data/PICO-500/energies.dat',np.array(energies))

