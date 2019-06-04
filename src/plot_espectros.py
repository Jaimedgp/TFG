import matplotlib
import numpy as np
import matplotlib.pyplot as plt
import os.path

#from Constantes import *
from simulacion import Simulacion

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

        fftWL, TFprom = dataPSD['fftWL'], dataPSD['TFprom']

    else:
        existData = False

        if saveIt:
            nWindw = 10
            laser = Simulacion(iBias[i], vRF, fR, nWindw)
            laser.allSimulation()
            laser.save()
        else:
            laser = Simulacion(iBias[i], vRF, fR)
            laser.allSimulation()

        fftWL, TFprom = laser.fftWL, laser.TFprom

    WLmax = fftWL[np.where(TFprom == max(TFprom))]
    WL.append(WLmax)

    plt.plot(fftWL, TFprom, label="%i mA" %(iBias[i]))
    #plt.text(WLmax, max(TFprom), "%.2f" %(WLmax))

plt.xlabel("$\lambda$ [nm]", fontsize=15)
plt.ylabel("PSD", fontsize=15)
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
##		Datos de la longitud de onda de emision del laser para
##	diferentes corrientes I_Bias obtenidos por simulacion con un
##	promedio a %s ventanas
##
##	Intensidad(mA)	Longitud de onda de pico(nm)
######################################################################\n"""
%(nWindw))
        for i in range(len(iBias)):
            fw.write("%i \t %.2f \n" %(iBias[i]*10**12, WL[i]))
