"""
        This script represent most representative data of the nonlinear dinamic
        of the semiconductor laser in sinusoidal wave (SW) with optical
        injection cahracterized by power injection (pwrInjct) and the detuning
        respect of the emission frequency (nuDetng)
"""

__author__ = 'Jaime Diez G-P'
__version__ = '2.0.0'
__email__ = "jaimediezgp@gmail.com"
__date__ = "Jun 26, 2019"

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
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
        'size'   : 15
       }
matplotlib.rc('font', **font)

########################################
##       LASER CHARACTERIZATION       ##
######################################## 

#----------------------------
#   Slave laser
#----------------------------

iBias = 35  # bias current [mA] / must be in [C ns^-1] by multiplying *10**-12
vRF = 1.0*10**(-9) #RMS voltage value of the signal generator [V]
fR = 5.0 # [GHz]

#----------------------------
#   Master laser
#----------------------------

pwrInjct = 100
# detuning of the injected laser field with respect to the emission frequency
nuDetng = -2 # [GHz]

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
    print "Open file " + nameFilePSD

    dataRateEq = np.load(nameFileRateEq)
    print "Open file " + nameFileRateEq

    time, S = dataRateEq['time'], dataRateEq['S']
    N, Phi = dataRateEq['N'], dataRateEq['Phi']
    fftWL, TFavg = dataPSD['fftWL'], dataPSD['TFavg']

else:
    laser = Simulation(iBias, vRF, fR, pwrInjct, nuDetng)
    laser.allSimulation()
    laser.save()

    time, S = laser.time, laser.S
    N, Phi = laser.N, laser.Phi
    fftWL, TFavg = laser.fftWL, laser.TFavg

#--------------------------------------------
#  Takes the values of the study inteval 
#--------------------------------------------
"""
indexes = np.where((time > transient) & (time < period+transient))
time = time[indexes]
N = N[indexes]
Phi = Phi[indexes] - 2*np.pi*(nuDetng - f0 + (c0/emissnWL[iBias]))*time
"""
power = constP * S *10**(12)
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

ax.plot(N, power, Phi)
plt.show()
