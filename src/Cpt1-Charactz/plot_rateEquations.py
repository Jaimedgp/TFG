"""
    Simulation of the rate equations used to discribe the dynamics of the
    carrier density (N(t)), the photon density (S(t)) and optical phase (phi(t))

    The simulation is done without taking into account the inyection terms

    author: JaimeDGP
    latest version: 29 April 2019
"""
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import os.path

import sys
sys.path.insert(0, '../')

from simulation import Simulation
from Constants import nTr

font = {
    'family' : 'serif',
    'weight' : 'normal',
    'size'   : 20
}

matplotlib.rc('font', **font)

colors = ['g', '#1f77b4', '#ff7f0e']
graphLabel = [
    ['(a)', '(b)', '(c)'],
    ['(d)', '(e)', '(f)'],
    ['(g)', '(h)', '(i)'],
    ['(j)', '(k)', '(l)']
]

iBias = 30  # bias current [mA] / must be in [C ns^-1] by multiplying *10**-12
vRF = [0, 0.05 *10**(-9), 1 *10**(-9)]#, 1.5 * 10**(-9)] #RMS voltage value of the signal generator [V]
fR = 5.0
period = 3 / fR

existData = False

fig, axs = plt.subplots(4, len(vRF), sharex=True, sharey="row",
                                                            figsize=(17, 10))

for i in range(len(vRF)):

    nameFile = ("Data/RateEquations_%imA_%imV_%iGHZ.npz"
                %(iBias, vRF[i] *10**(12), fR))

    if os.path.isfile(nameFile) and existData:

        dataPSD = np.load(nameFile)
        print "Opening file " + nameFile

        time = dataPSD['time']
        I = dataPSD['I']
        S = dataPSD['S']
        dPhi = dataPSD['dPhi']
        N = dataPSD['N']

    else:
        laser = Simulation(iBias, vRF[i], fR)
        laser.rateEquations()
        #laser.save()

        time, I, S, dPhi = laser.time, laser.I, laser.S, laser.dPhi
        N = laser.N/nTr

    indexes = np.where((time > 1.2) & (time < period+1.2))

    time = time[indexes]
    I = I[indexes]
    S = S[indexes]
    dPhi = dPhi[indexes]
    N = N[indexes]

    axs[0][i].set_title("$V_{RF} = $%.2f V" %(vRF[i] * 10**9), color=colors[i])

    axs[0][i].plot(time, I, colors[i])
    axs[0][i].axhline(y=14.8, linestyle=":", color='k', linewidth=3)
    axs[0][i].grid(linestyle='-.')
    axs[0][i].annotate(graphLabel[0][i], (0.9, 0.85),
                                            xycoords='axes fraction', size=20)

    axs[1][i].plot(time, S, colors[i])
    axs[1][i].grid(linestyle='-.')
    axs[1][i].annotate(graphLabel[1][i], (0.9, 0.85),
                                            xycoords='axes fraction', size=20)

    axs[2][i].plot(time, dPhi, colors[i], label="N(t)")
    axs[2][i].grid(linestyle='-.')
    axs[2][i].annotate(graphLabel[2][i], (0.9, 0.85),
                                            xycoords='axes fraction', size=20)

    axs[3][i].plot(time, N, colors[i])
    axs[3][i].grid(linestyle='-.')
    axs[3][i].annotate(graphLabel[3][i], (0.9, 0.85),
                                            xycoords='axes fraction', size=20)

axs[0][0].set_ylabel("I(t) [$mA$]")
axs[1][0].set_ylabel("S(t) [$m^{-3}$]")
axs[2][0].set_ylabel("Chirp [GHz]")
axs[3][0].set_ylabel("$N(t) / N_{Tr}$")

for i in range(len(vRF)):
    axs[-1][i].set_xlabel("t [ns]")

plt.tight_layout()
plt.show()
