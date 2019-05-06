"""
    Simulation of the rate equations used to discribe the dynamics of the
    carrier density (N(t)), the photon density (S(t)) and optical phase (phi(t))

    The simulation is done without taking into account the inyection terms

    author: JaimeDGP
    latest version: 29 April 2019
"""

import numpy as np
import matplotlib.pyplot as plt
import cmath

from Constantes import *
from getTempValues import getDeltaT

iBias = 30  # bias current [mA] / must be in [C ns^-1] by multiplying *10**-12
vRF = [0.05 *10**(-9), 1 *10**(-9), 1.5 * 10**(-9)] #RMS voltage value of the signal generator [V]
colors = ['g', 'b', 'r']

deltaT = getDeltaT(iBias)

faseTerm = faseConstant - pi2t * deltaT

delta = 0.0025 # tiempo de muestreo para la FFT [ns]
ndelta = int(delta / tIntev) # ndelta*tIntev=delta

periodos = 3 / fR
nPeriodo = int(periodos / delta)

tTotal = tTrans + periodos
nTotal = int(tTotal / tIntev)

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
I = np.zeros(nPeriodo)
N = np.zeros(nPeriodo)
S = np.zeros(nPeriodo)
dPhi = np.zeros(nPeriodo)

############################
##  Iniciar Simulacion
############################

#----------------------------------------
#   TERMS USED TO CALCULATE THE CHIRP
#----------------------------------------
derivAphvgTGmm = aphvgTGmm / tIntev
derivFaseTerm = faseTerm / tIntev
derivRuidoPhi = ruidoPhi / tIntev

fig, axs = plt.subplots(4, len(vRF), sharex=True, sharey="row",
                                                            figsize=(17, 10))
# Remove horizontal space between axes
fig.subplots_adjust(left=0.05, bottom=0.08, right=0.96, top=0.94, hspace=0.2)

timePeriod = np.linspace(tTrans, tTotal, nPeriodo)

for i in range(len(vRF)):
    inten = current(time, vRF[i])
    currentTerm = eVinv * inten

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
        sqrtS = np.sqrt(abs(tempS))

        tempPhi = (tempPhi + aphvgTGmm*tempN - faseTerm +
                                                    ruidoPhi*tempN*Y[q]/sqrtS)

        tempS = (tempS + vgTGmm*tempN*invS - vgTGmmN*invS - intTtau*tempS +
                                            btGmm*bTN + ruidoS*tempN*sqrtS*X[q])

        tempN = (tempN + currentTerm[q] - aTIntv*tempN - bTN -
                                (cTIntv*tempN**3) - vgT*tempN*invS + vgtN*invS)

    I[0] = inten[q]*10**(12)
    N[0] = tempN
    S[0] = tempS
    dPhi[0] = (1/(2*np.pi))*(derivAphvgTGmm*tempN - derivFaseTerm +
                                                derivRuidoPhi*tempN*Y[q]/sqrtS)

    for q in range(1, nPeriodo):
        for k in range(0, ndelta):

            index = (q-1)*ndelta + k + nTrans

            bTN = bTIntv * tempN * tempN
            invS = 1 / ((1/tempS) + epsilon)
            sqrtS = np.sqrt(tempS)

            tempPhi = (tempPhi + aphvgTGmm*tempN - faseTerm +
                                                ruidoPhi*tempN*Y[index]/sqrtS)

            tempS = (tempS + vgTGmm*tempN*invS - vgTGmmN*invS -
                        intTtau*tempS + btGmm*bTN + ruidoS*tempN*sqrtS*X[index])

            tempN = (tempN + currentTerm[index] - aTIntv*tempN - bTN -
                                (cTIntv*tempN**3) - vgT*tempN*invS + vgtN*invS)

        I[q] = inten[index]*10**(12)
        N[q] = tempN
        S[q] = tempS
        dPhi[q] = (1/(2*np.pi))*(derivAphvgTGmm*N[q] - derivFaseTerm +
                                    derivRuidoPhi*N[q]*Y[index]/np.sqrt(S[q]))

    #########################################
    ##  Representacion de los Datos
    #########################################

    axs[0][i].set_title("$V_{RF} = $%.2f V" %(vRF[i] * 10**9), color=colors[i], fontsize=15)

    axs[0][i].plot(timePeriod, I, colors[i])
    axs[0][i].axhline(y=14.8, linestyle=":", color='k', linewidth=3)
    axs[0][i].grid(linestyle='-.')

    axs[1][i].plot(timePeriod, S, colors[i])
    axs[1][i].grid(linestyle='-.')

    axs[2][i].plot(timePeriod, dPhi, colors[i], label="N(t)")
    axs[2][i].grid(linestyle='-.')

    axs[3][i].plot(timePeriod, N/nTr, colors[i])
    axs[3][i].grid(linestyle='-.')

axs[0][0].set_ylabel("I(t) [$mA$]", fontsize=15)
axs[1][0].set_ylabel("S(t) [$m^3$]", fontsize=15)
axs[2][0].set_ylabel("Chirp [GHz]", fontsize=15)
axs[3][0].set_ylabel("$N(t) / N_{Tr}$", fontsize=15)

for i in range(len(vRF)):
    axs[-1][i].set_xlabel("t [ns]", fontsize=15)

plt.show()
