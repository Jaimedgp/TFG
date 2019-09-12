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

colors = ['g', '#1f77b4', '#ff7f0e']
graphLabel = [
    ['(a)', '(b)', '(c)'],
    ['(d)', '(e)', '(f)'],
    ['(g)', '(h)', '(i)']
]

zones = ["P1", "P2", "CH-IR"]

psdLim = [[[1546.55, 1547.44], [2.505*10**(-13), 0.0112]],
          [[1546.31, 1547.76], [1.733*10**(-12), 0.0049]],
          [[1546.42, 1547.58], [9.037*10**(-13), 0.0014]]
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

pwrInjct = [50, 1000, 200]
# detuning of the injected laser field with respect to the emission frequency

nuDetng = 5 # [GHz]

#----------------------------
#   Study interval
#----------------------------

period = 6 / fR
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

fig, axs = plt.subplots(len(pwrInjct), 3, figsize=(17, 10))

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
        N = dataRateEq['N']
        fftWL, TFavg = dataPSD['fftWL'], dataPSD['TFavg']

    else:
        laser = Simulation(iBias, vRF, fR, pwrInjct[i], nuDetng)
        laser.rateEquations()

        time, S = laser.time, laser.S
        N = laser.N
        fftWL, TFavg = laser.fftWL, laser.TFavg

    #--------------------------------------------
    #  Takes the values of the study inteval 
    #--------------------------------------------

    indexes = np.where((time > 1.2) & (time < period+1.2))

    power = constP * S *10**(12)
    time = time[indexes]
    powerSp = power[indexes]
    NSp = N[indexes]

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

    axs[i][0].plot(fftWL, TFavg, colors[i])
    axs[i][0].set_yscale("log")
    axs[i][0].set_ylabel("PSD [u.a.]")
    axs[i][0].set_xlabel("$\lambda$ [nm]")
    axs[i][0].set_xlim(psdLim[i][0])
    axs[i][0].set_ylim(psdLim[i][1])
    axs[i][0].annotate(zones[i], (0.1, 0.70), xycoords='axes fraction', size=20)
    axs[i][0].annotate(graphLabel[i][0], (0.9, 0.85),
                                           xycoords='axes fraction', size=20)
    axs[i][0].annotate('', xy=(wLInject,yValue),
                            xytext=(wLInject,yStart),
                            arrowprops={'arrowstyle': '-|>', 'color':'r'},
                            va='center'
                        )


    axs[i][1].plot(time, powerSp, colors[i])
    axs[i][1].set_ylabel("P(t) [$\mu$W]")
    axs[i][1].set_xlabel("t [ns]")
    axs[i][1].grid(linestyle='-.')
    axs[i][1].annotate(graphLabel[i][1], (0.9, 0.85),
                                           xycoords='axes fraction', size=20)


    axs[i][2].plot(NSp, powerSp, colors[i])
    axs[i][2].set_xlabel("N(t)/$N_{Tr}$")
    axs[i][2].set_ylabel("P(t) [$\mu$W]")
    axs[i][2].annotate(graphLabel[i][2], (0.9, 0.85),
                                            xycoords='axes fraction', size=20)

    axs[0][1].get_shared_x_axes().join(axs[0][1], axs[1][1], axs[2][1])

    center = np.mean(axs[i][2].get_ylim())
    left = axs[i][2].get_xlim()[-1]
    axs[i][2].text(left, center, "$P_{Iny}$ = %i $\mu$W" %(pwrInjct[i]),
                   {'color': colors[i], 'fontsize': 20},
                   horizontalalignment='left',
                   verticalalignment='center',
                   rotation=-90,
                   clip_on=False)
plt.tight_layout()
plt.show()
