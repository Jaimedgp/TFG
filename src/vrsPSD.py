import matplotlib
import numpy as np
import matplotlib.pyplot as plt
import os.path

from simulacion import Simulacion
#from Constantes import *

font = {'family' : 'serif',
        'weight' : 'normal',
        'size'   : 15}
matplotlib.rc('font', **font)

iBias = 35  # bias current [mA] / must be in [C ns^-1] by multiplying *10**-12
vRF = [0.05 *10**(-9), 1 *10**(-9), 1.5 * 10**(-9)] #RMS voltage value of the signal generator [V]
fR = 5.0

existData = True

colors = ['g', 'b', 'r']

fig, axs = plt.subplots(1, len(vRF), sharex=True, sharey=True,
                                              figsize=(20, int(20/len(vRF))))
fig.subplots_adjust(left=0.05, bottom=0.08, right=0.96, top=0.94, hspace=0.2)

for i in range(len(vRF)):

    nameFile = "Data/PSD_%imA_%imV_%iGHZ.npz" %(iBias, vRF[i] *10**(12), fR)

    if os.path.isfile(nameFile) and existData:

        dataPSD = np.load(nameFile)
        print "Opening file: " + nameFile

        fftWL, TFprom = dataPSD['fftWL'], dataPSD['TFprom']

    else:
        laser = Simulacion(iBias, vRF[i], fR)
        laser.allSimulation()

        fftWL, TFprom = laser.fftWL, laser.TFprom

    #########################################
    ##  Representacion de los Datos
    #########################################

    #-------------------------------------------------------------------------------
    # Espectro optico frente a la longitud de onda
    #
    #   EL espectro optico se obtiene realizando la transformada de Fourier del
    #   campo optico total (opField). Las frecuencias han de ir en pasos de
    #   1/2Ndelta en el intervalo [-1/2delta, 1/2delta], siendo delta el paso de la
    #   transformada de fourier.  Para la obtencion de la longitud de onda se
    #   obtienen las frecuencias de la transformada de Fourier y se le suma la
    #   frecuencia total (freqTotal) y se pasa a longitud de onda con c0
    #-------------------------------------------------------------------------------

    axs[i].plot(fftWL, TFprom, colors[i])
    axs[i].set_xlabel("$\lambda$ [nm]", fontsize=15)
    axs[i].set_ylabel("PSD", fontsize=15)
    #axs[i].set_xlim([1546.0, 1547.5])
    #axs[i].set_ylim([9*10**(-15), 0.003])
    axs[i].set_yscale("log")
    axs[i].set_title("$V_{RF} = $ %.2f V" %(vRF[i]*10**9), color=colors[i])

plt.show()
