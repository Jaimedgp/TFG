import matplotlib
import numpy as np
import matplotlib.pyplot as plt
import os.path

import sys
sys.path.insert(0, '../')

from Constants import constP, f0, c0, nTr, pht2muWatt
from simulation import Simulation

font = {'family' : 'serif',
        'weight' : 'normal',
        'size'   : 15}
matplotlib.rc('font', **font)

iBias = 35  # bias current [mA] / must be in [C ns^-1] by multiplying *10**-12
vRF = 0.0 #RMS voltage value of the signal generator [V]
fR = 5.0
period = 6 / fR

sInjct = 1000 * pht2muWatt # 1 microW -> 3.8979848173986477 * 10**(17)
# detuning of the injected laser field with respect to the emission frequency
nuDetng = 6.2 # [GHz]

existData = True

#fig, axs = plt.subplots(2, 3, figsize=(20, 10))

powr = sInjct / pht2muWatt
nameFileRateEq = ("Data/RateEquations_%smuW_%sGHz"
                   %(powr, nuDetng)
                 )
nameFileRateEq = nameFileRateEq.replace('.',',')
nameFileRateEq += ".npz"

nameFilePSD = ("Data/PSD_%smuW_%iGHz"
                %(powr, nuDetng)
              )
nameFilePSD = nameFilePSD.replace('.',',')
nameFilePSD += ".npz"

if os.path.isfile(nameFilePSD) and existData:

    dataPSD = np.load(nameFilePSD)
    print "Open file " + nameFilePSD

    dataRateEq = np.load(nameFileRateEq)
    print "Open file " + nameFileRateEq

    time, S = dataRateEq['time'], dataRateEq['S']
    N, dPhi = dataRateEq['N'], dataRateEq['dPhi']
    fftWL, TFavg = dataPSD['fftWL'], dataPSD['TFavg']

else:
    laser = Simulation(iBias, vRF, fR, sInjct, nuDetng)
    laser.allSimulation()
    laser.save()

    time, S = laser.time, laser.S
    N, dPhi = laser.N, laser.dPhi
    fftWL, TFavg = laser.fftWL, laser.TFavg

indexes = np.where((time > 10.2) & (time < period+10.2))
time = time[indexes]
power = constP * S[indexes] *10**(12)
N = N[indexes]
dPhi = dPhi[indexes]

gridsize = (3, 6)
fig = plt.figure(figsize=(20, 10))
ax1 = plt.subplot2grid(gridsize, (0, 1), colspan=4, rowspan=2)

ax1.plot(fftWL, TFavg)
ax1.set_yscale("log")
ax1.set_xlabel("$\lambda$ [nm]")
ax1.set_ylabel("PSD")

ax2 = plt.subplot2grid(gridsize, (2, 0), colspan=2, rowspan=1)
ax2.plot(time, power)
ax2.set_xlabel("t [ns]", fontsize=15)
ax2.set_ylabel("Power [mW]", fontsize=15)

ax3 = plt.subplot2grid(gridsize, (2, 2), colspan=2, rowspan=1)
ax3.plot(time, N/nTr)
ax3.set_xlabel("t [ns]", fontsize=15)
ax3.set_ylabel("$N(t) / N_{Tr}$", fontsize=15)

ax4 = plt.subplot2grid(gridsize, (2, 4), colspan=2, rowspan=1)
ax4.plot(time, dPhi/(2*np.pi))
ax4.set_xlabel("t [ns]", fontsize=15)
ax4.set_ylabel("Chirp [GHz]", fontsize=15)

title = ("Inyeccion de Luz con $P_{inj} = $ %s" %(sInjct*constP*10**(15))
         +"$\mu$W y detuning $\delta \\nu = $%s GHz" %(nuDetng))
ax1.set_title(title, fontsize=25)
plt.tight_layout()
plt.show()
