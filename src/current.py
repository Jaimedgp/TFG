import matplotlib
import numpy as np
import matplotlib.pyplot as plt
import os.path

from Constantes import constP
from simulacion import Simulacion

font = {'family' : 'serif',
        'weight' : 'normal',
        'size'   : 15}
matplotlib.rc('font', **font)

iBias = [30, 50]  # bias current [mA] / must be in [C ns^-1] by multiplying *10**-12
vRF = 1.0 * 10**(-9) #RMS voltage value of the signal generator [V]
fR = 5.0
period = 3 / fR

existData = True

fig, axs = plt.subplots(1, 2, figsize=(17, 10))
# Remove horizontal space between axes
fig.subplots_adjust(left=0.05, bottom=0.08, right=0.96, top=0.94, hspace=0.2)

colors = ['b', 'r']

for i in range(len(iBias)):

    nameFileRateEq = "Data/RateEquations_%imA_%imV_%iGHZ.npz" %(iBias[i], vRF *10**(12), fR)
    nameFilePSD = "Data/PSD_%imA_%imV_%iGHZ.npz" %(iBias[i], vRF *10**(12), fR)

    if os.path.isfile(nameFilePSD) and existData:

        dataPSD = np.load(nameFilePSD)
        print "Opening file " + nameFilePSD

        dataRateEq = np.load(nameFileRateEq)
        print "Opening file " + nameFileRateEq

        time, S = dataRateEq['time'], dataRateEq['S']
        fftWL, TFprom = dataPSD['fftWL'], dataPSD['TFprom']

    else:
        laser = Simulacion(iBias[i], vRF, fR)
        laser.allSimulation()

        time, S = laser.time, laser.S
        fftWL, TFprom = laser.fftWL, laser.TFprom

    indexes = np.where((time > 1.2) & (time < period+1.2))
    time = time[indexes]
    power = constP * S[indexes] *10**(12)

    axs[0].plot(time, power, colors[i], label="%i mA" %(iBias[i]))
    axs[1].plot(fftWL, TFprom, colors[i], label="%i mA" %(iBias[i]))

axs[0].grid(linestyle='-.')
axs[0].set_xlabel("time [ns]", fontsize=15)
axs[0].set_xlim([1.2, 1.8])

axs[1].set_yscale("log")
axs[1].set_xlabel("WL [nm]", fontsize=15)
axs[1].set_xlim([1546.1, 1547.6])
axs[1].set_ylim([10**(-14), 0.7*10**(-3)])

axs[0].set_ylabel("Power [mW]", fontsize=15)
axs[1].set_ylabel("PSD", fontsize=15)

axs[0].legend()
axs[1].legend()
plt.show()

