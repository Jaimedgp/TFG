"""
        This script represent Power Spectra Density of the most representative
        case for each non linear zone .
"""

__author__ = 'Jaime Diez G-P'
__version__ = '2.0.0'
__email__ = "jaimediezgp@gmail.com"
__date__ = "Aug 14, 2019"

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import os.path

import sys
sys.path.insert(0, '../')

from Constants import *
from getDictValues import *
from simulation import Simulation

###################################################
##     CHANGES THE FONT APPEARANCE IN GRAPHS     ##
###################################################

font = {'family' : 'serif',
        'weight' : 'normal',
        'size'   : 20
       }
matplotlib.rc('font', **font)

zones = ['CH-IR', "IL", "P1"]
colors = ['g', '#1f77b4', '#ff7f0e']
graphLabel = [
    ['(a)', '(b)', '(c)'],
    ['(d)', '(e)', '(f)'],
    ['(g)', '(h)', '(i)']
]

psdLim = [[[1546.80, 1547.33], [5.920*10**(-11), 0.0409]],
          [[1546.80, 1547.35], [6.574*10**(-12), 0.0389]],
          [[1546.58, 1547.51], [2.427*10**(-12), 0.0718]],
         ]
########################################
##       LASER CHARACTERIZATION       ##
######################################## 

#----------------------------
#   Slave laser
#----------------------------

iBias = 35  # bias current [mA] / must be in [C ns^-1] by multiplying *10**-12
vRF = 0.0 #RMS voltage value of the signal generator [V]
fR = 5.0

#----------------------------
#   Master laser
#----------------------------

pwrInjct = [20, 100, 1000]
# detuning of the injected laser field with respect to the emission frequency

nuDetng = -2 # [GHz]

#----------------------------
#   Study interval
#----------------------------

period = 10 / fR
transient = 1.2

########################################
##          DATA COLLECTION           ##
######################################## 
#
#   Check if the data is already saved in a numpy file (.npz) in folder Data/
#   and if it exists, open the files and obtain the desired data. If it not
#   exists, execute the simulation script in order to calculate the data.
#
#   The data files are identify by the injection parameters in their names
#
#-----------------------------------------------------------------------------
existData = True

fig, axs = plt.subplots(3, len(pwrInjct), figsize=(17, 10))

for i in range(len(pwrInjct)):

    nameFileRateEq = ("Data/RateEquations_%smuW_%sGHz"
                    %(pwrInjct[i], nuDetng)
                    )
    nameFileRateEq = nameFileRateEq.replace('.',',')
    nameFileRateEq += ".npz"

    nameFilePSD = ("Data/PSD_%smuW_%sGHz"
                    %(pwrInjct[i], nuDetng)
                )
    nameFilePSD = nameFilePSD.replace('.',',')
    nameFilePSD += ".npz"

    if os.path.isfile(nameFilePSD) and existData:

        dataPSD = np.load(nameFilePSD)
        print "Opening file " + nameFilePSD

        dataRateEq = np.load(nameFileRateEq)
        print "Opening file " + nameFileRateEq

        time, S = dataRateEq['time'], dataRateEq['S']
        N, Phi = dataRateEq['N'], dataRateEq['Phi']
        fftWL, TFavg = dataPSD['fftWL'], dataPSD['TFavg']

    else:
        laser = Simulation(iBias, vRF, fR, pwrInjct[i], nuDetng)
        laser.rateEquations()

        time, S = laser.time, laser.S
        N, Phi = laser.N, laser.Phi
        fftWL, TFavg = laser.fftWL, laser.TFavg

    #--------------------------------------------
    #  Takes the values of the study inteval 
    #--------------------------------------------

    indexes = np.where((time > transient) & (time < period+transient))

    time = time[indexes]
    power = constP * S[indexes] *10**(12)
    N = N[indexes]
    Phi = Phi[indexes] - 2*np.pi*(nuDetng - f0 + (c0/(1.54705*10**(-6))))*time

########################################
##          PLOT DATA                 ##
######################################## 

    #-----------------------------------------------
    #  Mark with an arrow the injection wavelength
    #-----------------------------------------------

    # Injection wavelength
    wLInject = (c0*10**(9))/(c0/emissnWL[iBias] + nuDetng)

    indexes = np.argmax(np.where((fftWL > wLInject)))
    yValue = max(TFavg[indexes-5:indexes+5])
    length = np.log10(max(TFavg) - min(TFavg)) / 2.0
    yStart = yValue*10**(-length)

    axs[0][i].set_title("$P_{Iny}$ = %i $\mu$W" %(pwrInjct[i]), color=colors[i])
    axs[0][i].plot(fftWL, TFavg, colors[i])
    axs[0][i].set_yscale('log')
    axs[0][i].set_xlim(psdLim[i][0])
    axs[0][i].set_ylim(psdLim[i][1])
    axs[0][i].annotate(zones[i], (0.1, 0.70), xycoords='axes fraction', size=20)
    axs[0][i].annotate(graphLabel[0][i], (0.9, 0.85),
                                           xycoords='axes fraction', size=20)
    axs[0][i].annotate('', xy=(wLInject,yValue),
                            xytext=(wLInject,yStart),
                            arrowprops={'arrowstyle': '-|>', 'color':'r'},
                            va='center'
                        )

    axs[1][i].plot(time, power, colors[i])
    axs[1][i].grid(linestyle='-.')
    axs[1][i].annotate(graphLabel[1][i], (0.9, 0.85),
                                           xycoords='axes fraction', size=20)

    axs[2][i].plot(time, Phi, colors[i], label="N(t)")
    axs[2][i].grid(linestyle='-.')
    axs[2][i].annotate(graphLabel[2][i], (0.9, 0.85),
                                           xycoords='axes fraction', size=20)

axs[1][0].get_shared_x_axes().join(axs[1][0], axs[1][1], axs[1][2])
axs[2][0].get_shared_x_axes().join(axs[2][0], axs[2][1], axs[2][2])

axs[1][2].get_shared_y_axes().join(axs[1][2], axs[1][0], axs[1][1])

axs[0][0].set_ylabel("PSD [u.a.]")
axs[1][0].set_ylabel("P(t) [mW]")
axs[2][0].set_ylabel("$\Phi$")

for i in range(1, len(pwrInjct)):
    for j in range(0, 3):
        axs[i][j].set_xlabel("t [ns]")
        axs[i][j].set_xlim([1.1, 3.5])

for j in range(0, 3):
    axs[0][j].set_xlabel("$\lambda$ [nm]")

plt.tight_layout()
plt.show()
