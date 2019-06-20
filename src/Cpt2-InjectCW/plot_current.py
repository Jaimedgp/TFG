import matplotlib
import numpy as np
import matplotlib.pyplot as plt
import os.path

import sys
sys.path.insert(0, '../')

from Constants import constP, f0, c0
from simulation import Simulation

font = {'family' : 'serif',
        'weight' : 'normal',
        'size'   : 15}
matplotlib.rc('font', **font)

iBias = [35]  # bias current [mA] / must be in [C ns^-1] by multiplying *10**-12
vRF = 0.0 * 10**(-9) #RMS voltage value of the signal generator [V]
fR = 5.0
period = 3 / fR

sInjct = float(4 * 10**(19)) # 1 microW -> 4 * 10**(20)
# detuning of the injected laser field with respect to the emission frequency
nuDetng = -2 # [GHz]

existData = True

fig, axs = plt.subplots(1, 2, figsize=(20, 10))

colors = ['b', 'r']

for i in range(len(iBias)):

    oderMg = np.log10(sInjct / 4)
    nameFileRateEq = ("Data/RateEquations_%iS_%iGHz.npz"
                      %(oderMg, nuDetng)
                     )
    nameFilePSD = ("Data/PSD_%iS_%iGHz.npz"
                   %(oderMg, nuDetng)
                  )

    if os.path.isfile(nameFilePSD) and existData:

        dataPSD = np.load(nameFilePSD)
        print "Open file " + nameFilePSD

        dataRateEq = np.load(nameFileRateEq)
        print "Open file " + nameFileRateEq

        time, S = dataRateEq['time'], dataRateEq['S']
        fftWL, TFavg = dataPSD['fftWL'], dataPSD['TFavg']

    else:
        laser = Simulation(iBias[i], vRF, fR, sInjct, nuDetng)
        laser.allSimulation()
        laser.save()

        time, S = laser.time, laser.S
        fftWL, TFavg = laser.fftWL, laser.TFavg

    indexes = np.where((time > 10.2) & (time < period+10.2))
    time = time[indexes]
    power = constP * S[indexes] *10**(12)

    axs[0].plot(time, power, colors[i], label="%i mA" %(iBias[i]))
    axs[1].plot(fftWL, TFavg, colors[i], label="%i mA" %(iBias[i]))

axs[0].grid(linestyle='-.')
axs[0].set_xlabel("time [ns]", fontsize=15)
axs[0].set_xlim([10.2, 10.2+period])

axs[1].set_yscale("log")
axs[1].set_xlabel("WL [nm]", fontsize=15)
axs[1].set_xlim([1546.1, 1547.6])
axs[1].set_ylim([10**(-14), 0.7*10**(-3)])

axs[0].set_ylabel("Power [mW]", fontsize=15)
axs[1].set_ylabel("PSD", fontsize=15)

axs[0].legend()
axs[1].legend()
plt.tight_layout()
plt.show()

