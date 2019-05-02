import numpy as np
import matplotlib.pyplot as plt
import cmath

from Constantes import *
from getTempValues import *

iBias = 50  # bias current [mA] / must be in [C ns^-1] by multiplying *10**-12

deltaT = getDeltaT(iBias)
deltaF = getConstante(iBias)

faseTerm = faseConstant - pi2t * deltaT

vRF = [0.05 *10**(-9), 0.4 *10**(-9), 1.0 * 10**(-9), 1.2 * 10**(-9)]#, 1.5 * 10**(-9)] #RMS voltage value of the signal generator [V]

fR = 0.5 # [GHz]
rInt = 103.36 # [Ohm]

nWindw = 1 # numero de ventanas (para promediar) N natural

delta = 0.0025 # tiempo de muestreo para la FFT [ns]
nFFT = int(tWindw / delta) # numero de puntos de la FFT (potencia de 2)
ndelta = int(delta / tIntev) # ndelta*tIntev=delta

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
opField = np.zeros(nFFT, dtype=complex)
P = np.zeros(nFFT)
tmP = np.zeros(nFFT)

#-------------------------------------------------------------------------------
# Espectro optico frente a la longitud de onda
#
#   EL espectro optico se obtiene realizando la transformada de Fourier del
#   campo optico total (opField). Las frecuencias han de ir en pasos de
#   1/2Ndelta en el intervalo [-1/2delta, 1/2delta], siendo delta el paso de la
#   transformada de fourier.  Para la obtencion de la longitud de onda se
#   obtienen las frecuencias de la transformada de Fourier y se le suma la
#   frecuencia total (freqTotal) y se pasa a longitud de onda con c0
#-------------------------------------------------------------------------------
frecuencyLimits = 1 / (2*delta)
fftTime = np.linspace(-frecuencyLimits, frecuencyLimits, nFFT)#, endpoint=True)

fftTime += f0 - (deltaF/(2.0*np.pi))

fftWL = (c0/fftTime) *10**(9) # longitud de onda [nm]

############################
##  Iniciar Simulacion
############################

fig, axs = plt.subplots(2, len(vRF), figsize=(17, 10))
# Remove horizontal space between axes
fig.subplots_adjust(hspace=0.2)

for i in range(len(vRF)):
    inten = current(time, vRF[i])
    currentTerm = eVinv * inten
    TFprom = 0

    for win in range(0, nWindw):

        #---------------------------------------------------------
        # Vectores Gaussianos N(0,1) para el  Ruido
        #---------------------------------------------------------
        X = np.random.normal(0, 1, nTotal)
        Y = np.random.normal(0, 1, nTotal)

        # Se definen las condiciones iniciales para resolver la EDO
        tempN = nTr
        tempS = float(10**(20))
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
                                    cTIntv*tempN**3 - vgT*tempN*invS + vgtN*invS)

        opField[0] = np.sqrt(constP * tempS) * np.exp(1j*tempPhi)
        P[0] += (constP * tempS)/float(nWindw)
        tmP[0] = time[q]

        for q in range(1, nFFT):
            for k in range(0, ndelta):

                index = (q-1)*ndelta + k + nTrans

                bTN = bTIntv * tempN * tempN
                invS = 1 / ((1/tempS) + epsilon)
                sqrtS = np.sqrt(tempS)

                tempPhi = (tempPhi + aphvgTGmm*tempN - faseTerm +
                                            ruidoPhi*tempN*Y[index]/sqrtS)

                tempS = (tempS + vgTGmm*tempN*invS - vgTGmmN*invS - intTtau*tempS
                                + btGmm*bTN + ruidoS*tempN*sqrtS*X[index])

                tempN = (tempN + currentTerm[index] - aTIntv*tempN - bTN -
                                    (cTIntv*tempN**3) - vgT*tempN*invS + vgtN*invS)

            opField[q] = np.sqrt(constP * tempS) * np.exp(1j*tempPhi)
            P[q] += (constP * tempS)/float(nWindw)
            tmP[q] = time[index]

        transFourier = np.fft.fft(opField)
        TFprom += (abs(np.fft.fftshift(transFourier)) *
                abs(np.fft.fftshift(transFourier))/float(nWindw))

    #########################################
    ##  Representacion de los Datos
    #########################################

    axs[0][i].set_title("%.2f V" %(vRF[i] * 10**9), fontsize=15)

    axs[0][i].plot(tmP, P, 'r')
    axs[0][i].grid(linestyle='-.')
    axs[0][i].set_xlabel("time [ns]", fontsize=15)

    axs[1][i].plot(fftWL, TFprom, 'b')
    axs[1][i].set_yscale("log")
    axs[1][i].set_xlabel("WL [nm]", fontsize=15)

axs[0][0].set_ylabel("Power [$J ns^{-1}$]", fontsize=15)
axs[1][0].set_ylabel("PSD", fontsize=15)
plt.show()

