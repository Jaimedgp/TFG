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

colors = ['g', 'b', 'r']
graphLabel = [
    ['(a)', '(b)', '(c)', '(d)'],
    ['(e)', '(f)', '(g)', '(h)'],
    ['(i)', '(j)', '(k)', '(l)']
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

pwrInjct = [50, 200, 1000]
# detuning of the injected laser field with respect to the emission frequency

nuDetng = 5 # [GHz]

#----------------------------
#   Study interval
#----------------------------

period = 60 / fR
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

fig, axs = plt.subplots(len(pwrInjct), 4, sharex="row", sharey="row",
                                                            figsize=(17, 10))

for i in range(len(pwrInjct)):

    nameFileRateEq = ("Data/RateEquations_%smuW_%sGHz"
                    %(pwrInjct, nuDetng)
                    )
    nameFileRateEq = nameFileRateEq.replace('.',',')
    nameFileRateEq += ".npz"

    nameFilePSD = ("Data/PSD_%smuW_%sGHz"
                    %(pwrInjct, nuDetng)
                )
    nameFilePSD = nameFilePSD.replace('.',',')
    nameFilePSD += ".npz"

    if os.path.isfile(nameFilePSD) and existData:

        dataPSD = np.load(nameFilePSD)
        print "Opening file " + nameFilePSD

        dataRateEq = np.load(nameFileRateEq)
        print "Opening file " + nameFileRateEq

        time, S = dataRateEq['time'], dataRateEq['S']
        N = dataRateEq['N']
        fftWL, TFavg = dataPSD['fftWL'], dataPSD['TFavg']

    else:
        laser = Simulation(iBias, vRF[i], fR, sInyct, nuDetng)
        laser.rateEquations()

        time, S = laser.time, laser.S
        N = laser.N
        fftWL, TFavg = laser.fftWL, laser.TFavg

    #--------------------------------------------
    #  Takes the values of the study inteval 
    #--------------------------------------------

    indexes = np.where((time > 1.2) & (time < period+1.2))

    time = time[indexes]
    powerSp = constP * S[indexes] *10**(12)
    NSp = N[indexes]

########################################
##          PLOT DATA                 ##
######################################## 

    axs[i][0].plot(fftWL, TFavg, colors[i])
    axs[i][0].axhline(y=14.8, linestyle=":", color='k', linewidth=3)
    axs[i][0].grid(linestyle='-.')
    axs[i][0].annotate(graphLabel[0][i], (0.9, 0.85),
                                            xycoords='axes fraction', size=20)

    axs[i][1].plot(time, powerSp, colors[i])
    axs[i][1].grid(linestyle='-.')
    axs[i][1].annotate(graphLabel[1][i], (0.9, 0.85),
                                            xycoords='axes fraction', size=20)

    axs[i][2].plot(time, NSp, colors[i], label="N(t)")
    axs[i][2].grid(linestyle='-.')
    axs[i][2].annotate(graphLabel[2][i], (0.9, 0.85),
                                            xycoords='axes fraction', size=20)

    axs[i][3].plot(N, power colors[i])
    axs[i][3].grid(linestyle='-.')
    axs[i][3].annotate(graphLabel[3][i], (0.9, 0.85),
                                            xycoords='axes fraction', size=20)

#axs[0][0].set_ylabel("I(t) [$mA$]", fontsize=15)
#axs[1][0].set_ylabel("S(t) [$m^{-3}$]", fontsize=15)
#axs[2][0].set_ylabel("Chirp [GHz]", fontsize=15)
#axs[3][0].set_ylabel("$N(t) / N_{Tr}$", fontsize=15)

#for i in range(len(vRF)):
#    axs[-1][i].set_xlabel("t [ns]", fontsize=15)

plt.tight_layout()
plt.show()
