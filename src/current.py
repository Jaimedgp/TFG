import numpy as np
import matplotlib.pyplot as plt
import cmath

from Constantes import *
from getTempValues import *

iBias = [30, 50]  # bias current [mA] / must be in [C ns^-1] by multiplying *10**-12
colors = ['b', 'r']

vRF = 1.0 * 10**(-9) #RMS voltage value of the signal generator [V]

nWindw = 1
delta = 0.0025 # tiempo de muestreo para la FFT [ns]
nFFT = int(tWindw / delta) # numero de puntos de la FFT (potencia de 2)
ndelta = int(delta / tIntev) # ndelta*tIntev=delta

#                  INTENSIDAD
#
#                 2 sqrt(2) vRF 
# I_bias + cLoss --------------- sin(2 pi fR t)
#                   z0 + zL
current = lambda t, currnt: (currnt*10**(-12) + (cLoss * 2.0 * np.sqrt(2) * vRF *
                                            np.sin(2 * np.pi * fR * t)) / rInt)

################################################################################
##  Inicializar los vectores de tiempo (time), de la densidad de portadores (N)
##      densidad de fotones (S) y de la fase optica (Phi)
################################################################################

time = np.linspace(0, tTotal, nTotal)

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

opField = np.zeros(nFFT, dtype=complex)
tmP = np.zeros(nFFT)

############################
##  Iniciar Simulacion
############################

fig, axs = plt.subplots(1, 2, figsize=(17, 10))
# Remove horizontal space between axes
fig.subplots_adjust(hspace=0.2)

for i in range(len(iBias)):
    deltaT = getDeltaT(iBias[i])
    deltaF = getConstante(iBias[i])

    faseTerm = faseConstant - pi2t * deltaT

    inten = current(time, iBias[i])
    currentTerm = eVinv * inten
    P = np.zeros(nFFT)
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

            tempPhi = (tempPhi + aphvgTGmm*tempN - faseTerm +
                                            ruidoPhi*tempN*Y[q]/np.sqrt(abs(tempS)))

            tempS = (tempS + vgTGmm*tempN*invS - vgTGmmN*invS - intTtau*tempS +
                                btGmm*bTN + ruidoS*tempN*np.sqrt(abs(tempS))*X[q])

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

                tempPhi = (tempPhi + aphvgTGmm*tempN - faseTerm +
                                            ruidoPhi*tempN*Y[index]/np.sqrt(tempS))

                tempS = (tempS + vgTGmm*tempN*invS - vgTGmmN*invS - intTtau*tempS
                                + btGmm*bTN + ruidoS*tempN*np.sqrt(tempS)*X[index])

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
    fftTimeI = fftTime + f0 - (deltaF/(2.0*np.pi))

    fftWL = (c0/fftTimeI) *10**(9) # longitud de onda [nm]

    axs[0].plot(tmP, P, colors[i], label="%i mA" %(iBias[i]))
    axs[1].plot(fftWL, TFprom, colors[i], label="%i mA" %(iBias[i]))

axs[0].grid(linestyle='-.')
axs[0].set_xlabel("time [ns]", fontsize=15)

axs[1].set_yscale("log")
axs[1].set_xlabel("WL [nm]", fontsize=15)

axs[0].set_ylabel("Power [$J ns^{-1}$]", fontsize=15)
axs[1].set_ylabel("PSD", fontsize=15)

axs[0].legend()
axs[1].legend()
plt.show()

