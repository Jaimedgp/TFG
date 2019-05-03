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

iBias = 30  # bias current [mA] / must be in [C ns^-1] by multiplying *10**-12
vRF = 0 #RMS voltage value of the signal generator [V]

deltaT = getDeltaT(iBias)

faseTerm = faseConstant - pi2t * deltaT

tIntev = 1 *10**(-5) # tiempo de integracion [ns]
delta = 0.0025 # tiempo de muestreo para la FFT [ns]
ndelta = int(delta / tIntev)

tWindw = 1.2
nTime = int(tWindw / delta) # numero de pasos de integracion

tTrans = 0.2
nTrans = int(tTrans / delta)

tTotal = tTrans + tWindw
nTotal = int(tTotal / tIntev)
nTotalD = int(tTotal / delta)

#                  INTENSIDAD
#
#                 2 sqrt(2) vRF 
# I_bias + cLoss --------------- sin(2 pi fR t)
#                   z0 + zL
current = lambda t, vRFi: (iBias*10**(-12) + (cLoss * 2.0 * np.sqrt(2) * vRFi *
                                            np.sin(2 * np.pi * fR * t)) / rInt)

################################################################################
##  Inicializar los vectores de tiempo (time), de la densidad de portadores (N)
##      densidad de fotones (S) y de la fase optica (Phi)
################################################################################

time = np.linspace(0, tTotal, nTotal)
I = np.zeros(nTotalD)
N = np.zeros(nTotalD)
S = np.zeros(nTotalD)
dPhi = np.zeros(nTotalD)

############################
##  Iniciar Simulacion
############################

derivAphvgTGmm = aphvgTGmm / tIntev
derivFaseTerm = faseTerm / tIntev
derivRuidoPhi = ruidoPhi / tIntev

fig, axs = plt.subplots(4, 1, sharex=True, figsize=(10, 20))
fig.subplots_adjust(left=0.08, bottom=0.06, right=0.96, top=0.94, hspace=0.1)

inten = current(time, vRF)
currentTerm = eVinv * inten

#---------------------------------------------------------
# Vectores Gaussianos N(0,1) para el  Ruido
#---------------------------------------------------------

X = np.random.normal(0, 1, nTotal)
Y = np.random.normal(0, 1, nTotal)

# Se definen las condiciones iniciales para resolver la EDO
tempN = nTr
tempS = float(10**15)
tempPhi = 0

for q in range(0, nTrans):
    for j in range(0, ndelta):
        totalIndex = q*ndelta + j

        bTN = bTIntv * tempN * tempN
        invS = 1 / ((1/tempS) + epsilon)
        sqrtS = np.sqrt(abs(tempS))

        tempPhi = (tempPhi + aphvgTGmm*tempN - faseTerm +
                                            ruidoPhi*tempN*Y[totalIndex]/sqrtS)

        tempS = (tempS + vgTGmm*tempN*invS - vgTGmmN*invS - intTtau*tempS +
                                btGmm*bTN + ruidoS*tempN*sqrtS*X[totalIndex])

        tempN = (tempN + currentTerm[totalIndex] - aTIntv*tempN - bTN -
                                (cTIntv*tempN**3) - vgT*tempN*invS + vgtN*invS)
    I[q] = iBias
    N[q] = tempN
    S[q] = tempS
    dPhi[q] = (1/(2*np.pi))*(derivAphvgTGmm*tempN - derivFaseTerm +
                                        derivRuidoPhi*tempN*Y[totalIndex]/sqrtS)

for q in range(nTrans, nTotalD):
    for k in range(0, ndelta):

        index = (q-nTrans)*ndelta + k

        bTN = bTIntv * tempN * tempN
        invS = 1 / ((1/tempS) + epsilon)
        sqrtS = np.sqrt(tempS)

        tempPhi = (tempPhi + aphvgTGmm*tempN - faseTerm +
                                                ruidoPhi*tempN*Y[index]/sqrtS)

        tempS = (tempS + vgTGmm*tempN*invS - vgTGmmN*invS - intTtau*tempS +
                                        btGmm*bTN + ruidoS*tempN*sqrtS*X[index])

        tempN = (tempN + currentTerm[index] - aTIntv*tempN - bTN -
                                (cTIntv*tempN**3) - vgT*tempN*invS + vgtN*invS)

    I[q] = iBias
    N[q] = tempN
    S[q] = tempS
    dPhi[q] = (1/(2*np.pi))*(derivAphvgTGmm*N[q] - derivFaseTerm +
                                    derivRuidoPhi*N[q]*Y[index]/np.sqrt(S[q]))
I[0] = 0

#########################################
##  Representacion de los Datos
#########################################

timePeriod = np.linspace(0, tTotal, nTotalD)

axs[0].plot(timePeriod, I, 'r')
axs[0].axhline(y=13.2, linestyle=":", color='k', linewidth=3)
axs[0].grid(linestyle='-.')

axs[1].plot(timePeriod, S, 'b')
axs[1].set_yscale('log')
axs[1].grid(linestyle='-.')

axs[2].plot(timePeriod, N/nTr, 'r')
axs[2].grid(linestyle='-.')

axs[3].plot(timePeriod, dPhi, 'r', label="N(t)")
axs[3].grid(linestyle='-.')
axs[3].set_ylim([-40, 20])

axs[0].set_ylabel("I(t) [$mA$]", fontsize=15)
axs[1].set_ylabel("S(t) [$m^3$]", fontsize=15)
axs[2].set_ylabel("$N(t) / N_{Tr}$", fontsize=15)
axs[3].set_ylabel("Chirp [GHz]", fontsize=15)

axs[3].set_xlabel("t [ns]", fontsize=15)

plt.show()
