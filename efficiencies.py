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

# DARWIN

bins = [i for i in range(5,21,3)]

print bins

energies = [0.1*j + 2 for j in range(241)]

sigmaE = lambda ER: 0.05 * ER + np.sqrt(0.05 * ER)

efficiencies = np.array([[sum([0.3 * GeneralAcceptanceFunction(energies[j], bins[i], bins[i+1], sigmaE(energies[j]))for i in range(len(bins) - 1)])] + [0.3 * GeneralAcceptanceFunction(energies[j], bins[i], bins[i+1], sigmaE(energies[j]))for i in range(len(bins) - 1)] for j in range(len(energies))])

np.savetxt('data/DARWIN/efficiencies.dat',efficiencies)

np.savetxt('data/DARWIN/energies.dat',np.array(energies))

# DarkSide

acceptance = np.loadtxt('data/DarkSide/DarkSide_f200.dat')

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

energies_old = [3.2,     3.39712, 3.49709, 3.6,     3.70595, 3.81501, 3.92728, 4.04286, 4.16184, 4.28431, 4.4104,  4.54019, 4.6738,  4.81135, 4.95294, 5.0987,  5.24875, 5.40322, 5.56223, 5.72592, 5.89442, 6.06789, 6.24646, 6.43029, 6.61953, 6.81433, 7.01487, 7.22131, 7.43383, 7.6526,  7.8778,  8.10964, 8.3483,  8.59398, 8.84689, 9.10724, 9.37526, 9.65117, 9.93519, 10.2276, 10.5286, 10.8384, 11.1574, 11.4857, 11.8237, 12.1717, 12.5299, 12.8986, 13.2782, 13.669,  14.0712, 14.4853, 14.9116, 15.3505, 15.8022, 16.2673, 16.746,  17.2388, 17.7461, 18.2684, 18.806,  19.3594, 19.9292, 20.5157, 21.1194, 21.7409, 22.3807, 23.0394, 23.7174, 24.4154, 25.1339, 25.8736, 26.635,  27.4188, 28.2258, 29.0564, 29.9115, 30.7918, 31.6979, 32.6308, 33.5911, 34.5796, 35.5973, 36.6448, 37.7233, 38.8334, 39.9762, 41.1527, 42.3638, 43.6105, 44.8939, 46.2151, 47.5752, 48.9752, 50.4165, 51.9002, 53.4276, 54.9999, 56.6185, 58.2847, 60.]

efficiencies_old = [0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.01260, 0.18145, 0.32697, 0.44480, 0.57556, 0.75532, 0.90435, 1.00000, 1.00000, 1.00000, 1.00000, 1.00000, 1.00000, 1.00000, 1.00000, 1.00000, 1.00000, 1.00000, 1.00000, 1.00000, 1.00000, 1.00000, 1.00000, 1.00000, 1.00000, 1.00000, 1.00000, 1.00000, 1.00000, 1.00000, 1.00000, 1.00000, 1.00000, 1.00000, 1.00000, 1.00000, 1.00000, 1.00000, 1.00000, 1.00000, 1.00000, 1.00000, 1.00000, 1.00000, 1.00000, 1.00000, 1.00000, 1.00000, 1.00000, 1.00000, 1.00000, 1.00000, 1.00000, 1.00000, 1.00000, 1.00000, 1.00000, 1.00000, 1.00000, 1.00000, 1.00000, 1.00000, 1.00000, 1.00000, 1.00000, 1.00000, 1.00000, 1.00000, 1.00000, 1.00000, 1.00000, 1.00000, 1.00000, 1.00000, 1.00000, 1.00000, 1.00000, 1.00000, 1.00000, 1.00000, 1.00000, 1.00000, 1.00000]

acceptancefunc = interp1d(energies_old, efficiencies_old)

energies = [0.1*j + 3.2 for j in range(68)] + [0.4*j + 10 for j in range(101)] 

efficiency = [[0.5 * acceptancefunc(energies[j]) for j in range(68)] + [1. for j in range(101)], [0.5  * acceptancefunc(energies[j]) for j in range(68)] + [0 for j in range(101)], [0 for j in range(68)] + [1.0 for j in range(101)]]

np.savetxt('data/PICO-500/efficiencies.dat',np.transpose(np.array(efficiency)))

np.savetxt('data/PICO-500/energies.dat',np.array(energies))

