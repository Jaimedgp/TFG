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
vRF = [0.05 *10**(-9), 1 *10**(-9)]#, 1.5 * 10**(-9)] #RMS voltage value of the signal generator [V]

deltaT = getDeltaT(int(iBias * 10**(12)))

faseTerm = faseConstant - pi2t * deltaT

nWindw = 1 # numero de ventanas (para promediar) N natural
delta = 0.0025 # tiempo de muestreo para la FFT [ns]
ndelta = int(delta / tIntev) # ndelta*tIntev=delta

#nTrans = int(tTrans / delta)

periodos = 3 / fR
nPeriodo = int(periodos / delta)

tTotal = tTrans + periodos
nTotal = int(tTotal / tIntev)


#                  INTENSIDAD
#
#                 2 sqrt(2) vRF 
# I_bias + cLoss --------------- sin(2 pi fR t)
#                   z0 + zL
current = lambda t, vRFi: (iBias + (cLoss * 2.0 * np.sqrt(2) * vRFi * np.sin(2
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

derivAphvgTGmm = aphvgTGmm / tIntev
derivFaseTerm = faseTerm / tIntev
derivRuidoPhi = ruidoPhi / tIntev

fig, axs = plt.subplots(len(vRF), 4, sharex=True, figsize=(17, 10))
# Remove horizontal space between axes
#fig.subplots_adjust(vspace=1)
fig.subplots_adjust(hspace=0.2)
timePeriod = np.linspace(tTrans, tTotal, nPeriodo)

for i in range(len(vRF)):
    inten = current(time, vRF[i])
    currentTerm = eVinv * inten

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
        I[0] = inten[q]*10**(12)
        N[0] = tempN
        S[0] = tempS
        dPhi[0] = (1/(2*np.pi))*(derivAphvgTGmm*tempN - derivFaseTerm +
                                    derivRuidoPhi*tempN*Y[q]/np.sqrt(abs(tempS)))

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

            I[q] = inten[index]*10**(12)
            N[q] = tempN
            S[q] = tempS
            dPhi[q] = (1/(2*np.pi))*(derivAphvgTGmm*N[q] - derivFaseTerm +
                                            derivRuidoPhi*N[q]*Y[index]/np.sqrt(S[q]))

    #########################################
    ##  Representacion de los Datos
    #########################################

    axs[i][0].set_ylabel("%.2f V" %(vRF[i] * 10**9), fontsize=15)

    axs[i][0].plot(timePeriod, I, 'r')
    axs[i][0].axhline(y=13.2, linestyle=":", color='k', linewidth=3)
    axs[i][0].grid(linestyle='-.')

    axs[i][1].plot(timePeriod, S, 'b')
    axs[i][1].set_yscale('log')
    axs[i][1].grid(linestyle='-.')

    axs[i][2].plot(timePeriod, dPhi, 'r', label="N(t)")
    axs[i][2].grid(linestyle='-.')
    #axs[2][0].tick_params('y', colors='r')

    axs[i][3].plot(timePeriod, N/nTr, 'r')
    axs[i][3].grid(linestyle='-.')

axs[0][0].set_title("I(t) [$mA$]", fontsize=15)
axs[0][1].set_title("S(t) [$m^3$]", fontsize=15)
axs[0][2].set_title("Chirp [GHz]", fontsize=15)
axs[0][3].set_title("$N(t) / N_{Tr}$", fontsize=15)

axs[-1][0].set_xlabel("t [ns]", fontsize=15)
axs[-1][1].set_xlabel("t [ns]", fontsize=15)
axs[-1][2].set_xlabel("t [ns]", fontsize=15)
axs[-1][3].set_xlabel("t [ns]", fontsize=15)

#fig.savefig("./Graficas/rateEquations.png", dpi=500)
plt.show()
