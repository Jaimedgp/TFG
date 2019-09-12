import matplotlib
import numpy as np
import matplotlib.pyplot as plt
import os.path

import sys
sys.path.insert(0, '../')

from simulation import Simulation

font = {'family' : 'serif',
        'weight' : 'normal',
        'size'   : 15}
matplotlib.rc('font', **font)

iBias = [i for i in range(15, 40, 5) ] # bias current [C ns^-1]
vRF = 0 #RMS voltage value of the signal generator [V]
fR = 5.0
existData = True
saveIt = False

WL = []

fig = plt.figure(figsize=(11,8))

for i in range(len(iBias)):

    nameFilePSD = "Data/PSD_%imA_%imV_%iGHZ.npz" %(iBias[i], vRF *10**(12), fR)

    if os.path.isfile(nameFilePSD) and existData:

        dataPSD = np.load(nameFilePSD)
        print "Opening file " + nameFilePSD

        fftWL, TFavg = dataPSD['fftWL'], dataPSD['TFavg']

    else:
        existData = False

        if saveIt:
            nWindw = 10
            laser = Simulation(iBias[i], vRF, fR, nWindw)
            laser.allSimulation()
            laser.save()
        else:
            laser = Simulation(iBias[i], vRF, fR)
            laser.allSimulation()

        fftWL, TFavg = laser.fftWL, laser.TFavg

    WLmax = fftWL[np.where(TFavg == max(TFavg))]
    WL.append(WLmax)

    plt.plot(fftWL, TFavg, label="%i mA" %(iBias[i]))

plt.xlabel("$\lambda$ [nm]")
plt.ylabel("PSD [u.a.]")
plt.yscale("log")
plt.xlim(1546.8, 1547.2)
plt.ylim(0.9*10**(-12), 0.003)
plt.legend()
plt.tight_layout()
plt.show()
#fig.savefig("./Graficas/Espectros.png", dpi=300)
if not existData and saveIt:
    with open("../Datos/Landas.txt", 'w') as fw:
        fw.write(
"""######################################################################
##		Emission wavelength data for the laser at different currents
##	I_Bias
##
##	Intensity(mA)	Emission Wavelength(nm)
######################################################################\n"""
)
        for i in range(len(iBias)):
            fw.write("%i \t %.2f \n" %(iBias[i]*10**12, WL[i]))
