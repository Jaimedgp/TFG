"""
    Simulation of the rate equations used to discribe the dynamics of the
    carrier density (N(t)), the photon density (S(t)) and optical phase (phi(t))
    in the transient interval.

    The simulation is done without taking into account the inyection terms

"""

__author__ = 'Jaime Diez G-P'
__version__ = '1.0.0'
__email__ = "jaimediezgp@gmail.com"
__date__ = "Jun 13, 2019"

import matplotlib
import numpy as np
import matplotlib.pyplot as plt
import os.path

import sys
sys.path.insert(0, '../')

from simulation import Simulation
from Constants import nTr

font = {'family' : 'serif',
        'weight' : 'normal',
        'size'   : 15}
matplotlib.rc('font', **font)

colors = ['g', '#1f77b4', '#ff7f0e', 'r']

iBias = 30  # bias current [mA] / must be in [C ns^-1] by multiplying *10**-12
vRF = 0 #RMS voltage value of the signal generator [V]
fR = 5.0

existData = False

fig, axs = plt.subplots(4, 1, sharex=True, figsize=(10, 10))

nameFile = "Data/RateEquations_%imA_%imV_%iGHZ.npz" %(iBias, vRF *10**(12), fR)

if os.path.isfile(nameFile) and existData:

    dataPSD = np.load(nameFile)
    print "Opening file " + nameFile

    time = dataPSD['time']
    I = dataPSD['I']
    S = dataPSD['S']
    dPhi = dataPSD['dPhi']
    N = dataPSD['N']

else:
    laser = Simulation(iBias, vRF, fR)
    laser.rateEquations()

    time, I, S, dPhi, N = laser.time, laser.I, laser.S, laser.dPhi, laser.N/nTr

indexes = np.where(time < 1.2)

time = time[indexes]
I = I[indexes]
S = S[indexes]
dPhi = dPhi[indexes]
N = N[indexes]

axs[0].plot(time, I, colors[0])
axs[0].axhline(y=13.2, linestyle=":", color='k', linewidth=3)
axs[0].grid(linestyle='-.')

axs[1].plot(time, S, colors[1])
axs[1].set_yscale('log')
axs[1].grid(linestyle='-.')

axs[2].plot(time, N, colors[2])
axs[2].grid(linestyle='-.')

axs[3].plot(time, dPhi,colors[3], label="N(t)")
axs[3].grid(linestyle='-.')
axs[3].set_ylim([-40, 20])

axs[0].set_ylabel("I(t) [$mA$]")
axs[1].set_ylabel("S(t) [$m^{-3}$]")
axs[2].set_ylabel("$N(t) / N_{Tr}$")
axs[3].set_ylabel("Chirp [GHz]")

axs[3].set_xlabel("t [ns]")

plt.tight_layout()
plt.show()
