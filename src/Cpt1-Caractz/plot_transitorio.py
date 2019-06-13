"""
    Simulation of the rate equations used to discribe the dynamics of the
    carrier density (N(t)), the photon density (S(t)) and optical phase (phi(t))

    The simulation is done without taking into account the inyection terms

    author: JaimeDGP
    latest version: 19 April 2019
"""

import matplotlib
import numpy as np
import matplotlib.pyplot as plt
import os.path

import sys
sys.path.insert(0, '../')

from Constantes import nTr
from simulacion import Simulacion

font = {'family' : 'serif',
        'weight' : 'normal',
        'size'   : 15}
matplotlib.rc('font', **font)

iBias = 30  # bias current [mA] / must be in [C ns^-1] by multiplying *10**-12
vRF = 0 #RMS voltage value of the signal generator [V]
fR = 5.0

existData = False

fig, axs = plt.subplots(4, 1, sharex=True, figsize=(10, 10))

nameFile = "Data/RateEquations_%imA_%imV_%iGHZ.npz" %(iBias, vRF *10**(12), fR)

if os.path.isfile(nameFile) and existData:

    dataPSD = np.load(nameFile)
    print "Opening file " + nameFile

    time, I, S, dPhi, N = dataPSD['time'], dataPSD['I'], dataPSD['S'], dataPSD['dPhi'], dataPSD['N']

else:
    laser = Simulacion(iBias, vRF, fR)
    laser.rateEquations()

    time, I, S, dPhi, N = laser.time, laser.I, laser.S, laser.dPhi, laser.N

indexes = np.where(time < 1.2)

time = time[indexes]
I = I[indexes]
S = S[indexes]
dPhi = dPhi[indexes]
N = N[indexes]

axs[0].plot(time, I, 'r')
axs[0].axhline(y=13.2, linestyle=":", color='k', linewidth=3)
axs[0].grid(linestyle='-.')

axs[1].plot(time, S, 'b')
axs[1].set_yscale('log')
axs[1].grid(linestyle='-.')

axs[2].plot(time, N/nTr, 'r')
axs[2].grid(linestyle='-.')

axs[3].plot(time, dPhi, 'r', label="N(t)")
axs[3].grid(linestyle='-.')
axs[3].set_ylim([-40, 20])

axs[0].set_ylabel("I(t) [$mA$]", fontsize=15)
axs[1].set_ylabel("S(t) [$m^{-3}$]", fontsize=15)
axs[2].set_ylabel("$N(t) / N_{Tr}$", fontsize=15)
axs[3].set_ylabel("Chirp [GHz]", fontsize=15)

axs[3].set_xlabel("t [ns]", fontsize=15)

plt.tight_layout()
plt.show()
