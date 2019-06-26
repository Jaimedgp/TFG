"""
        This script represent most representative data of the nonlinear dinamic
        of the semiconductor laser in continuous wave (CW) with optical
        injection cahracterized by power injection (pwrInjct) and the detuning
        respect of the emission frequency (nuDetng)

        The squeme used to plot the data is the following:

                            +-------------------------+
                            |                         |
                            |                         |
                            |          PSD            |
                            |                         |
                            |                         |
                            +-------------------------+

                +--------------+  +--------------+  +--------------+
                |              |  |              |  |              |
                |  power / N   |  |   N / time   |  | Phase / time |
                |              |  |              |  |              |
                +--------------+  +--------------+  +--------------+
"""

__author__ = 'Jaime Diez G-P'
__version__ = '2.0.0'
__email__ = "jaimediezgp@gmail.com"
__date__ = "Jun 26, 2019"

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import os.path

import sys
sys.path.insert(0, '../')

from Constants import constP, f0, c0, pht2muWatt
from simulation import Simulation

###################################################
##     CHANGES THE FONT APPEARANCE IN GRAPHS     ##
###################################################

font = {'family' : 'serif',
        'weight' : 'normal',
        'size'   : 15}
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

indexes = np.where((time > transient) & (time < period+transient))
time = time[indexes]
power = constP * S[indexes] *10**(12)
N = N[indexes]
Phi = Phi[indexes] - 2*np.pi*(nuDetng - f0 + (c0/(1.54705*10**(-6))))*time

########################################
##          PLOT DATA                 ##
######################################## 

gridsize = (3, 6)
fig = plt.figure(figsize=(20, 10))
ax1 = plt.subplot2grid(gridsize, (0, 1), colspan=4, rowspan=2)

ax1.plot(fftWL, TFavg)
ax1.set_yscale("log")
ax1.set_xlabel("$\lambda$ [nm]")
ax1.set_ylabel("PSD")
ax1.set_xlim(1546.5, 1547.5)

ax2 = plt.subplot2grid(gridsize, (2, 0), colspan=2, rowspan=1)
ax2.plot(N, power)
ax2.set_xlabel("$N(t) / N_{Tr}$", fontsize=15)
ax2.set_ylabel("Power [mW]", fontsize=15)

ax3 = plt.subplot2grid(gridsize, (2, 2), colspan=2, rowspan=1)
ax3.plot(time, N)
ax3.set_xlabel("t [ns]", fontsize=15)
ax3.set_ylabel("$N(t) / N_{Tr}$", fontsize=15)

ax4 = plt.subplot2grid(gridsize, (2, 4), colspan=2, rowspan=1)
ax4.plot(time, Phi/(2*np.pi))
ax4.set_xlabel("t [ns]", fontsize=15)
ax4.set_ylabel("Chirp [GHz]", fontsize=15)

title = ("Inyeccion de Luz con $P_{inj} = $ %s" %(pwrInjct)
         +"$\mu$W y detuning $\delta \\nu = $%s GHz" %(nuDetng))
ax1.set_title(title, fontsize=25)

plt.tight_layout()
plt.show()
