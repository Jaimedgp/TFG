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
vRF = 0.0 #RMS voltage value of the signal generator [V]
fR = 5.0

#----------------------------
#   Master laser
#----------------------------

pwrInjct = [1, 20, 100, 200, 300, 1000, 8000, 20000]
# detuning of the injected laser field with respect to the emission frequency
nuDetng = -2 # [GHz]

########################################
##       CREATE THREADS               ##
######################################## 

jobs = []
lasers = []

for i in range(len(pwrInjct)):
    lasers.append(Simulation(iBias, vRF, fR, pwrInjct[i], nuDetng))

    thread = threading.Thread(target=lasers[i].allSimulation())
    jobs.append(thread)

#----------------------------
#   Master laser
#----------------------------

pwrInjct = [50, 200, 1000]
# detuning of the injected laser field with respect to the emission frequency
nuDetng = 5 # [GHz]

for i in range(len(pwrInjct), len(pwrInjct)+2):
    lasers.append(Simulation(iBias, vRF, fR, pwrInjct[i], nuDetng))

    thread = threading.Thread(target=lasers[i].allSimulation())
    jobs.append(thread)

#start the threads (i.e. calculate the random number lists)
for j in jobs:
    j.start()

# Ensure all of the threads have finished
for j in jobs:
    j.join()


for i in range(len(lasers)):
    lasers[i].save()
