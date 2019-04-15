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
vRF = 1.0 *10**(-9) #RMS voltage value of the signal generator [V]

deltaT = getDeltaT(int(iBias * 10**(12)))

faseTerm = faseConstant - pi2t * deltaT

nWindw = 1 # numero de ventanas (para promediar) N natural
delta = 0.0025 # tiempo de muestreo para la FFT [ns]
ndelta = int(delta / tIntev) # ndelta*tIntev=delta

#nTrans = int(tTrans / delta)

periodos = 4 / fR
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
I = np.zeros(nPeriodo)
N = np.zeros(nPeriodo)
S = np.zeros(nPeriodo)
dPhi = np.zeros(nPeriodo)

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
        bTN = bTIntv * tempN * tempN

        invS = 1 / ((1/tempS) + epsilon)

        tempPhi = (tempPhi + aphvgTGmm*tempN - faseTerm +
                                ruidoPhi*tempN*Y[q]/np.sqrt(abs(tempS)))

        tempS = (tempS + vgTGmm*tempN*invS - vgTGmmN*invS - intTtau*tempS +
                    btGmm*bTN + ruidoS*tempN*np.sqrt(abs(tempS))*X[q])

        tempN = (tempN + currentTerm[q] - aTIntv*tempN - bTN -
                            (cTIntv*tempN**3) - vgT*tempN*invS + vgtN*invS)
    I[0] = currentTerm[q]
    N[0] = tempN
    S[0] = tempS
    dPhi[0] = tempPhi

    for q in range(1, nPeriodo):
        for k in range(0, ndelta):

            index = (q-1)*ndelta + k + nTrans

            bTN = bTIntv * tempN * tempN

            invS = 1 / ((1/tempS) + epsilon)

            tempPhi = (tempPhi + aphvgTGmm*tempN - faseTerm +
                                        ruidoPhi*tempN*Y[index]/np.sqrt(tempS))

            tempS = (tempS + vgTGmm*tempN*invS - vgTGmmN*invS - intTtau*tempS +
                              btGmm*bTN + ruidoS*tempN*np.sqrt(tempS)*X[index])

            tempN = (tempN + currentTerm[index] - aTIntv*tempN - bTN -
                                (cTIntv*tempN**3) - vgT*tempN*invS + vgtN*invS)

        I[q] = currentTerm[index]
        N[q] += tempN/float(nWindw)
        S[q] += tempS/float(nWindw)
        dPhi[q] += tempPhi/float(nWindw)

#########################################
##  Representacion de los Datos
#########################################

timePeriod = np.linspace(tTrans, tTotal, nPeriodo)

fig, axs = plt.subplots(5, 1, sharex=True)
# Remove horizontal space between axes
fig.subplots_adjust(hspace=0)

axs[0].plot(timePeriod, I, 'r')
axs[0].set_ylabel("I(t) [$C ns^-1$]", color='r', fontsize=15)

axs[1].plot(timePeriod, N, 'r')
axs[1].set_ylabel("N(t) [$m^3$]", color='r', fontsize=15)

axs[2].plot(timePeriod, S, 'b')
axs[2].set_ylabel("S(t) [$m^3$]", color='b', fontsize=15)
axs[2].set_yscale('log')

axs[3].plot(timePeriod, N, 'r', label="N(t)")
axs[3].set_ylabel("N(t) [$m^3$]", color='r', fontsize=15)
axs[3].tick_params('y', colors='r')

ax4 = axs[3].twinx()
ax4.plot(timePeriod, S, 'b', label="S(t)")
ax4.set_ylabel("S(t) [$m^3$]", color='b', fontsize=15)
ax4.set_yscale('log')
ax4.tick_params('y', colors='b')

axs[4].plot(timePeriod, dPhi, 'r', label="N(t)")
axs[4].set_ylabel("$\Phi$(t) [$m^3$]", color='r', fontsize=15)
axs[4].tick_params('y', colors='r')
plt.show()
