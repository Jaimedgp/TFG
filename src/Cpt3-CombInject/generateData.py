"""
        This script represent Power Spectra Density of the most representative
        case for each non linear zone .
"""

__author__ = 'Jaime Diez G-P'
__version__ = '2.0.0'
__email__ = "jaimediezgp@gmail.com"
__date__ = "Aug 14, 2019"

import sys
sys.path.insert(0, '../')

from simulation import Simulation

########################################
##       LASER CHARACTERIZATION       ##
######################################## 

#----------------------------
#   Slave laser
#----------------------------

iBias = 35  # bias current [mA] / must be in [C ns^-1] by multiplying *10**-12
vRF = 1.0*10**(-9) #RMS voltage value of the signal generator [V]
fR = 5.0

#----------------------------
#   Master laser
#----------------------------

pwrInjct = [1, 100, 200, 300, 1000, 20000]
# detuning of the injected laser field with respect to the emission frequency
nuDetng = -2 # [GHz]

########################################
##       CREATE THREADS               ##
######################################## 

jobs = []
lasers = []

for i in range(len(pwrInjct)):
    lasers.append(Simulation(iBias, vRF, fR, pwrInjct[i], nuDetng))

    lasers[i].allSimulation()

for i in range(len(lasers)):
    lasers[i].save()
