__author__ = 'Jaime Diez G-P'
__version__ = '1.0.0'
__email__ = "jaimediezgp@gmail.com"
__date__ = "May 27, 2019"

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import os.path

from simulacion import Simulacion

font = {'family' : 'serif',
        'weight' : 'normal',
        'size'   : 15}
matplotlib.rc('font', **font)

iBias = 30  # bias current [mA] / must be in [C ns^-1] by multiplying *10**-12
vRF = 1.0 *10**(-9) #RMS voltage value of the signal generator [V]
fR = 5.0

existData = True

nameFile = "Data/PSD_%imA_%imV_%iGHZ.npz" %(iBias, vRF *10**(12), fR)

if os.path.isfile(nameFile) and existData:

    dataPSD = np.load(nameFile)
    print "Opening file: " + nameFile

    fftWL = dataPSD['fftWL']
    TFprom = dataPSD['TFprom']
    TFang = dataPSD['TFang']

else:
    laser = Simulacion(iBias, vRF, fR)
    laser.allSimulation()

    fftWL, TFprom, TFang = laser.fftWL, laser.TFprom, laser.TFang

fig, ax1 = plt.subplots(figsize=(8,6))

ax1.plot(fftWL, TFang, "r")
ax1.set_xlabel("$\lambda$ [nm]", fontsize=15)
ax1.set_ylabel("FFT Phase [rad]", fontsize=15)
ax1.tick_params('y', colors='r')

ax2 = ax1.twinx()
ax2.plot(fftWL, TFprom, "b")
ax2.set_ylabel("PSD", fontsize=15)
ax2.set_yscale("log")
ax2.tick_params('y', colors='b')
fig.tight_layout()
plt.show()
