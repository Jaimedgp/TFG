import matplotlib
import numpy as np
import matplotlib.pyplot as plt
import os.path

import sys
sys.path.insert(0, '../')

from Constants import constP
from simulation import Simulation

font = {'family' : 'serif',
        'weight' : 'normal',
        'size'   : 20}
matplotlib.rc('font', **font)

iBias = [35, 50]  # bias current [mA] / must be in [C ns^-1] by multiplying *10**-12
vRF = 1.0 * 10**(-9) #RMS voltage value of the signal generator [V]
fR = 5.0
period = 3 / fR

existData = True

fig, axs = plt.subplots(1, 2, figsize=(20, 10))

colors = ['#1f77b4', '#ff7f0e']

for i in range(len(iBias)):

    nameFileRateEq = ("Data/RateEquations_%imA_%imV_%iGHZ.npz"
                      %(iBias[i], vRF *10**(12), fR))
    nameFilePSD = "Data/PSD_%imA_%imV_%iGHZ.npz" %(iBias[i], vRF *10**(12), fR)

    if os.path.isfile(nameFilePSD) and existData:

        dataPSD = np.load(nameFilePSD)
        print "Open file " + nameFilePSD

        dataRateEq = np.load(nameFileRateEq)
        print "Open file " + nameFileRateEq

        time, S = dataRateEq['time'], dataRateEq['S']
        fftWL, TFavg = dataPSD['fftWL'], dataPSD['TFavg']

    else:
        laser = Simulation(iBias[i], vRF, fR)
        laser.allSimulation()

        time, S = laser.time, laser.S
        fftWL, TFavg = laser.fftWL, laser.TFavg

    indexes = np.where((time > 1.2) & (time < period+1.2))
    time = time[indexes]
    power = constP * S[indexes] *10**(12)

    axs[0].plot(time, power, colors[i], label="%i mA" %(iBias[i]))
    axs[1].plot(fftWL, TFavg, colors[i], label="%i mA" %(iBias[i]))

axs[0].grid(linestyle='-.')
axs[0].set_xlabel("t [ns]")
axs[0].set_xlim([1.2, 1.8])

axs[1].set_yscale("log")
axs[1].set_xlabel("$\lambda$ [nm]")
axs[1].set_xlim([1546.1, 1547.6])
axs[1].set_ylim([10**(-14), 0.7*10**(-3)])

axs[0].set_ylabel("P(t) [mW]")
axs[1].set_ylabel("PSD")

axs[0].legend()
axs[1].legend()
plt.tight_layout()
plt.show()

