"""
     Optical spectra againe a la longitud de onda

       Optical spectra is obtained by the Fast Fourier transform (FFT) of the
       total optical field (opField). Frequency has to go by steps of 1/2Ndelta in
       interval [-1/2delta, 1/2delta], being delta the FFT's step.In order to
       obtain the wavelength, the FFT's frequencies have been obtained and
       added the total frequency (freqTotal). Frequencies are transform to
       wavelength by c0.
"""
__author__ = 'Jaime Diez G-P'
__version__ = '1.0.0'
__email__ = "jaimediezgp@gmail.com"
__date__ = "Jun 13, 2019"

import matplotlib
import matplotlib.ticker as ticker
from matplotlib.ticker import FormatStrFormatter
import numpy as np
import matplotlib.pyplot as plt
import os.path

import sys
sys.path.insert(0, '../')

from simulation import Simulation

font = {'family' : 'serif',
        'weight' : 'normal',
        'size'   : 20}
matplotlib.rc('font', **font)

# bias current [mA] / must be in [C ns^-1] by multiplying *10**-12
iBias = 30
#RMS voltage value of the signal generator [V]
vRF = [0.05 *10**(-9), 1 *10**(-9), 1.5 * 10**(-9)]
fR = 5.0

existData = True

colors = ['g', '#1f77b4', '#ff7f0e']

fig, axs = plt.subplots(1, len(vRF), sharex=True, sharey=True,
                                              figsize=(20, int(20/len(vRF))))

for i in range(len(vRF)):

    nameFile = "Data/PSD_%imA_%imV_%iGHZ.npz" %(iBias, vRF[i] *10**(12), fR)

    if os.path.isfile(nameFile) and existData:

        dataPSD = np.load(nameFile)
        print "Opening file: " + nameFile

        fftWL, TFavg = dataPSD['fftWL'], dataPSD['TFavg']

    else:
        sys.stdout.write("Simulating a (%imA %imV %iGHz) laser....." %(iBias,
                                                                       vRF[i]
                                                                    *10**(12),
                                                                       fR))
        sys.stdout.flush()
        laser = Simulation(iBias, vRF[i], fR)
        laser.allSimulation()
        print "Done"
        #sys.stdout.write("Done")
        #sys.stdout.flush()

        fftWL, TFavg = laser.fftWL, laser.TFavg


    axs[i].plot(fftWL, TFavg, colors[i])
    axs[i].set_xlabel("$\lambda$ [nm]")
    axs[i].set_ylabel("PSD")
    axs[i].set_xlim([1546.0, 1547.75])
    axs[i].xaxis.set_major_locator(ticker.MultipleLocator(0.5))
    axs[i].xaxis.set_minor_locator(ticker.MultipleLocator(0.1))
    axs[i].set_yscale("log")
    axs[i].set_title("$V_{RF} = $ %.2f V" %(vRF[i]*10**9), color=colors[i])

plt.tight_layout()
plt.show()
