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

zones = ["P1", "CH-IR", "IL", "IL", "P1", "P1", "IL", "IL"]
labels = ["(a)", "(b)", "(c)", "(d)", "(e)", "(f)", "(g)", "(h)"]
pwrInjct = [1, 20, 100, 200, 300, 1000, 8000, 20000]
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

fig = plt.figure()

existData = True

for i in range(len(pwrInjct)):

    nameFilePSD = ("Data/PSD_%smuW_%sGHz"
                    %(pwrInjct[i], nuDetng)
                )
    nameFilePSD = nameFilePSD.replace('.',',')
    nameFilePSD += ".npz"

    if os.path.isfile(nameFilePSD) and existData:

        dataPSD = np.load(nameFilePSD)
        print "Open file " + nameFilePSD

        fftWL, TFavg = dataPSD['fftWL'], dataPSD['TFavg']

    else:
        laser = Simulation(iBias, vRF, fR, pwrInjct[i], nuDetng)
        laser.allSimulation()
        laser.save()

        fftWL, TFavg = laser.fftWL, laser.TFavg

    ########################################
    ##          PLOT DATA                 ##
    ######################################## 

    plt.subplot(2, 1, i+1)
    #axs.set_title("%.2f V" %(vRF[i] * 10**9), fontsize=15)
    plt.plot(fftWL, TFavg, 'b')
    plt.yscale("log")
    plt.xlabel("WL [nm]")
    plt.ylabel("PSD")
    plt.annotate(labels[i], (0.9, 0.85), xycoords='axes fraction', size=20)
    plt.annotate(zones[i], (0.9, 0.85), xycoords='axes fraction', size=20)

plt.tight_layout()
plt.show()
