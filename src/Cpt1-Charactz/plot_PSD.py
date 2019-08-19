__author__ = 'Jaime Diez G-P'
__version__ = '1.0.0'
__email__ = "jaimediezgp@gmail.com"
__date__ = "May 27, 2019"

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import os.path

import sys
sys.path.insert(0, '../')

from simulation import Simulation

font = {'family' : 'serif',
        'weight' : 'normal',
        'size'   : 20}
matplotlib.rc('font', **font)

iBias = 35  # bias current [mA] / must be in [C ns^-1] by multiplying *10**-12
vRF = 0.0 *10**(-9) #RMS voltage value of the signal generator [V]
fR = 5.0

existData = True

nameFile = "Data/PSD_%imA_%imV_%iGHZ.npz" %(iBias, vRF *10**(12), fR)

if os.path.isfile(nameFile) and existData:

    dataPSD = np.load(nameFile)
    print "Open file: " + nameFile

    fftWL, TFavg = dataPSD['fftWL'], dataPSD['TFavg']

else:
    laser = Simulation(iBias, vRF, fR)
    laser.allSimulation()

    fftWL, TFavg = laser.fftWL, laser.TFavg

fig = plt.figure(figsize=(8,6))
fig.subplots_adjust(left=0.05, bottom=0.08, right=0.96, top=0.94, hspace=0.2)

plt.plot(fftWL, TFavg)
plt.xlabel("$\lambda$ [nm]")
plt.ylabel("PSD")
plt.yscale("log")
plt.title("$I_{Bias}$ = %i mA \t $V_{RF} = $ %.1f V" %(iBias, vRF*10**9))
plt.tight_layout()
plt.show()
