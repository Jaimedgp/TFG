"""
    Simulation of the rate equations used to discribe the dynamics of the
    carrier density (N(t)), the photon density (S(t)) and optical phase (phi(t))

    The simulation is done without taking into account the inyection terms

    author: JaimeDGP
    latest version: 19 April 2019
"""

import numpy as np
import matplotlib.pyplot as plt
import cmath

from Constantes import *
from getTempValues import getDeltaT

iBias = 35 *10**(-12) # bias current [C ns^-1]
vRF = 1.2 *10**(-9) #RMS voltage value of the signal generator [V]

deltaT = getDeltaT(int(iBias * 10**(12)))

faseTerm = faseConstant - pi2t * deltaT

nWindw = 1 # numero de ventanas (para promediar) N natural
delta = 0.0025 # tiempo de muestreo para la FFT [ns]
ndelta = int(delta / tIntev) # ndelta*tIntev=delta

nTrans = int(tTrans / delta)

periodos = 3 / fR
nPeriodo = int(periodos / delta)

tTotal = tTrans + periodos
nTotal = int(tTotal / tIntev)


#                  INTENSIDAD
#
#                 2 sqrt(2) vRF 
# I_bias + cLoss --------------- sin(2 pi fR t)
#                   z0 + zL
current = lambda t: (iBias + (cLoss * 2.0 * np.sqrt(2) * vRF * np.sin(2
                                                * np.pi * fR * t)) / (z0 + zL))

################################################################################
##  Inicializar los vectores de tiempo (time), de la densidad de portadores (N)
##      densidad de fotones (S) y de la fase optica (Phi)
################################################################################

time = np.linspace(0, tTotal, nTotal)
N = np.zeros(nPeriodo)
S = np.zeros(nPeriodo)
Phi = np.zeros(nPeriodo)

############################
##  Iniciar Simulacion
############################

TFprom = 0

currentTerm = eVinv * current(time)

for win in range(0, nWindw):

    #---------------------------------------------------------
    # Vectores Gaussianos N(0,1) para el  Ruido
    #---------------------------------------------------------

    X = np.random.normal(0, 1, nTotal)
    Y = np.random.normal(0, 1, nTotal)

    # Se definen las condiciones iniciales para resolver la EDO
    tempN = nTr
    tempS = float(10**20)
    tempPhi = 0

    for q in range(0, nTrans):
        for k in range(0, ndelta):

            index = (q*ndelta) + k

            bTN = bTIntv * tempN * tempN

            invS = 1 / ((1/tempS) + epsilon)

            tempPhi = (tempPhi + aphvgTGmm*tempN - faseTerm +
                                    ruidoPhi*tempN*Y[index]/np.sqrt(abs(tempS)))

            tempS = (tempS + vgTGmm*tempN*invS - vgTGmmN*invS - intTtau*tempS +
                        btGmm*bTN + ruidoS*tempN*np.sqrt(abs(tempS))*X[index])

            tempN = (tempN + currentTerm[index] - aTIntv*tempN - bTN -
                                (cTIntv*tempN**3) - vgT*tempN*invS + vgtN*invS)

    N[0] = tempN
    S[0] = tempS
    Phi[0] = tempPhi
    usedIndex = nTrans*ndelta

    for q in range(1, nPeriodo):
        for k in range(0, ndelta):

            index = (q-1)*ndelta + k + usedIndex

            bTN = bTIntv * tempN * tempN

            invS = 1 / ((1/tempS) + epsilon)

            tempPhi = (tempPhi + aphvgTGmm*tempN - faseTerm +
                                        ruidoPhi*tempN*Y[index]/np.sqrt(tempS))

            tempS = (tempS + vgTGmm*tempN*invS - vgTGmmN*invS - intTtau*tempS +
                              btGmm*bTN + ruidoS*tempN*np.sqrt(tempS)*X[index])

            tempN = (tempN + currentTerm[index] - aTIntv*tempN - bTN -
                                (cTIntv*tempN**3) - vgT*tempN*invS + vgtN*invS)

        N[q] += tempN/float(nWindw)
        S[q] += tempS/float(nWindw)
        Phi[q] += tempPhi/float(nWindw)

#########################################
##  Representacion de los Datos
#########################################

timePeriod = np.linspace(tTrans, tTotal, nPeriodo)

fig, ax1 = plt.subplots()
ax1.plot(timePeriod, N, 'r', label="N(t)")
ax1.set_xlabel("tiempo [ns]", fontsize=20)
ax1.set_xlim(0, 1)
ax1.set_ylabel("N(t) [$m^3$]", color='r', fontsize=20)
ax1.tick_params('y', colors='r')

ax2 = ax1.twinx()
ax2.plot(timePeriod, S, 'b', label="S(t)")
ax2.set_ylabel("S(t) [$m^3$]", color='b', fontsize=20)
ax2.set_yscale('log')
ax2.tick_params('y', colors='b')
fig.tight_layout()
plt.show()
