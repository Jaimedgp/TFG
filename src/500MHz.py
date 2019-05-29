import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import os.path

font = {'family' : 'serif',
        'weight' : 'normal',
        'size'   : 15}
matplotlib.rc('font', **font)

from Constantes import constP
from simulacion import Simulacion

iBias = 50  # bias current [mA] / must be in [C ns^-1] by multiplying *10**-12
vRF = [0.05 *10**(-9), 0.4 *10**(-9), 1.0 * 10**(-9), 1.2 * 10**(-9)]#, 1.5 * 10**(-9)] #RMS voltage value of the signal generator [V]
fR = 0.5 # [GHz]
period = 3 / fR

existData = True

fig, axs = plt.subplots(2, len(vRF), figsize=(17, 10))
# Remove horizontal space between axes
fig.subplots_adjust(left=0.05, bottom=0.08, right=0.98, top=0.94, hspace=0.2)

for i in range(len(vRF)):

    nameFileRateEq = "Data/RateEquations_%imA_%imV_%iGHZ.npz" %(iBias, vRF[i] *10**(12), fR)
    nameFilePSD = "Data/PSD_%imA_%imV_%iGHZ.npz" %(iBias, vRF[i] *10**(12), fR)

    if os.path.isfile(nameFilePSD) and existData:

        dataPSD = np.load(nameFilePSD)
        print "Opening file " + nameFilePSD

        dataRateEq = np.load(nameFileRateEq)
        print "Opening file " + nameFileRateEq

        time, S = dataRateEq['time'], dataRateEq['S']
        fftWL, TFprom = dataPSD['fftWL'], dataPSD['TFprom']

    else:
        laser = Simulacion(iBias, vRF[i], fR)
        laser.allSimulation()

        time, S = laser.time, laser.S
        fftWL, TFprom = laser.fftWL, laser.TFprom

    indexes = np.where((time > 1.2) & (time < period+1.2))
    time = time[indexes]
    power = constP * S[indexes] *10**(12)

    #########################################
    ##  Representacion de los Datos
    #########################################

    axs[0][i].set_title("%.2f V" %(vRF[i] * 10**9), fontsize=15)

    axs[0][i].plot(time, power, 'r')
    axs[0][i].grid(linestyle='-.')
    axs[0][i].set_xlabel("time [ns]", fontsize=15)
    axs[0][i].set_xlim([1.5, 7.5])

    axs[1][i].plot(fftWL, TFprom, 'b')
    axs[1][i].set_yscale("log")
    axs[1][i].set_xlabel("WL [nm]", fontsize=15)
    axs[1][i].set_xlim([1547.15, 1547.37])
    axs[1][i].set_ylim(6.82*10**(-11))

axs[0][0].set_ylabel("Power [mW]", fontsize=15)
axs[1][0].set_ylabel("PSD", fontsize=15)
plt.show()
