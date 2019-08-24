import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import os.path

import sys
sys.path.insert(0, '../')

font = {'family' : 'serif',
        'weight' : 'normal',
        'size'   : 20}
matplotlib.rc('font', **font)

from Constants import constP
from simulation import Simulation

# bias current [mA] / must be in [C ns^-1] by multiplying *10**-12
iBias = 50
#RMS voltage value of the signal generator [V]
vRF = [0.05 *10**(-9), 0.4 *10**(-9), 1.0 * 10**(-9), 1.2 * 10**(-9)]
fR = 0.5 # [GHz]
period = 3 / fR

existData = True

fig, axs = plt.subplots(2, len(vRF), figsize=(20, 10))

for i in range(len(vRF)):

    nameFileRateEq = ("Data/RateEquations_%imA_%imV_%iGHZ.npz"
                      %(iBias, vRF[i] *10**(12), fR))
    nameFilePSD = ("Data/PSD_%imA_%imV_%iGHZ.npz"
                   %(iBias, vRF[i] *10**(12), fR))

    if os.path.isfile(nameFilePSD) and existData:

        dataPSD = np.load(nameFilePSD)
        print "Open file " + nameFilePSD

        dataRateEq = np.load(nameFileRateEq)
        print "Open file " + nameFileRateEq

        time, S = dataRateEq['time'], dataRateEq['S']
        fftWL, TFavg = dataPSD['fftWL'], dataPSD['TFavg']

    else:
        laser = Simulation(iBias, vRF[i], fR)
        laser.allSimulation()

        time, S = laser.time, laser.S
        fftWL, TFavg = laser.fftWL, laser.TFavg

    indexes = np.where((time > 1.2) & (time < period+1.2))
    time = time[indexes]
    power = constP * S[indexes] *10**(12)

    #########################################
    ##  Representacion de los Datos
    #########################################

    axs[0][i].set_title("%.2f V" %(vRF[i] * 10**9))

    axs[0][i].plot(time, power, '#ff7f0e')
    axs[0][i].grid(linestyle='-.')
    axs[0][i].set_xlabel("t [ns]")
    axs[0][i].set_xlim([1.5, 7.5])

    axs[1][i].plot(fftWL, TFavg, '#1f77b4')
    axs[1][i].set_yscale("log")
    axs[1][i].set_xlabel("$\lambda$ [nm]")
    axs[1][i].set_xlim([1547.15, 1547.37])
    axs[1][i].set_ylim(6.82*10**(-11))

axs[0][0].set_ylabel("P(t) [mW]")
axs[1][0].set_ylabel("PSD")

axs[0][3].get_shared_y_axes().join(axs[0][0], axs[0][1], axs[0][2], axs[0][3])

plt.tight_layout()
plt.show()
