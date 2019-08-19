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

from matplotlib.patches import Arrow
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

########################################
##       LASER CHARACTERIZATION       ##
######################################## 

zones = [["P1", "CH-IR", "IL", "IL"],
         ["P1", "P1", "IL", "IL"]
        ]
labels = [["(a)", "(b)", "(c)", "(d)"],
          ["(e)", "(f)", "(g)", "(h)"]
         ]

psdLim = [[
           [[1546.98, 1547.14], [1.729*10**(-11), 0.0034]],#1 microW
           [[1546.88, 1547.22], [5.446*10**(-11), 0.0220]],#20 microW
           [[1546.75, 1547.34], [4.191*10**(-12), 0.0579]],#100 microW
           [[1546.77, 1547.32], [3.607*10**(-12), 0.0591]]
          ], [
           [[1546.70, 1547.40], [3.964*10**(-12), 0.0339]],#300 microW
           [[1546.47, 1547.54], [3.250*10**(-12), 0.0389]],#1000 microW
           [[1546.54, 1547.55], [3.376*10**(-12), 0.0637]],#8000 microW
           [[1546.69, 1547.38], [2.909*10**(-12), 0.1132]]#20000 microW
          ]
         ]

#----------------------------
#   Slave laser
#----------------------------

iBias = 35  # bias current [mA] / must be in [C ns^-1] by multiplying *10**-12
vRF = 0.0 #RMS voltage value of the signal generator [V]
fR = 5.0

#----------------------------
#   Master laser
#----------------------------

pwrInjct = [[1, 20, 100, 200],
            [300, 1000, 8000, 20000]
           ]
# detuning of the injected laser field with respect to the emission frequency
nuDetng = -2 # [GHz]

##----------------------------
##   Study interval
##----------------------------
#
#period = 60 / fR
#transient = 1.2

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

fig, axs = plt.subplots(len(pwrInjct), len(pwrInjct[0]), figsize=(17, 10))

existData = True

for i in range(len(pwrInjct)):
    for j in range(len(pwrInjct[0])):

        nameFilePSD = ("Data/PSD_%smuW_%sGHz"
                        %(pwrInjct[i][j], nuDetng)
                    )
        nameFilePSD = nameFilePSD.replace('.',',')
        nameFilePSD += ".npz"

        if os.path.isfile(nameFilePSD) and existData:

            dataPSD = np.load(nameFilePSD)
            print "Open file " + nameFilePSD

            fftWL, TFavg = dataPSD['fftWL'], dataPSD['TFavg']

        else:
            laser = Simulation(iBias, vRF, fR, pwrInjct[i][j], nuDetng)
            laser.allSimulation()
            #laser.save()

            fftWL, TFavg = laser.fftWL, laser.TFavg

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

        arrow = Arrow(x=wLInject, y=yStart,
                    dx=0, dy=yValue-yStart,
                    width=0.02,
                    color='c'
                    )

        axs[i][j].set_title("$P_{Iny}$ = %i $\mu$W" %(pwrInjct[i][j]))
        axs[i][j].plot(fftWL, TFavg, color='b')
        axs[i][j].set_yscale('log')
        axs[i][j].add_patch(arrow)
        axs[i][j].set_xlabel("$\lambda$ [nm]")
        axs[i][j].set_ylabel("PSD")
        axs[i][j].set_xlim(psdLim[i][j][0])
        axs[i][j].set_ylim(psdLim[i][j][1])
        axs[i][j].annotate(zones[i][j], (0.1, 0.70), xycoords='axes fraction', size=20)
        axs[i][j].annotate(labels[i][j], (0.9, 0.85),
                                            xycoords='axes fraction', size=20)

plt.tight_layout()
plt.show()
