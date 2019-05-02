import numpy as np
import matplotlib.pyplot as plt
import cmath

from Constantes import *
from getTempValues import *

iBias = 30  # bias current [mA] / must be in [C ns^-1] by multiplying *10**-12

deltaT = getDeltaT(iBias)
deltaF = getConstante(iBias)

faseTerm = faseConstant - pi2t * deltaT

vRF = 1.0 *10**(-9) #RMS voltage value of the signal generator [V]

nWindw = 10 # numero de ventanas (para promediar) N natural

delta = 0.0025 # tiempo de muestreo para la FFT [ns]
nFFT = int(tWindw / delta) # numero de puntos de la FFT (potencia de 2)
ndelta = int(delta / tIntev) # ndelta*tIntev=delta

#                  INTENSIDAD
#
#                 2 sqrt(2) vRF 
# I_bias + cLoss --------------- sin(2 pi fR t)
#                   z0 + zL
current = lambda t: (iBias*10**(-12) + (cLoss * 2.0 * np.sqrt(2) * vRF *
                                            np.sin(2 * np.pi * fR * t)) / rInt)

################################################################################
##  Inicializar los vectores de tiempo (time), de la densidad de portadores (N)
##      densidad de fotones (S) y de la fase optica (Phi)
################################################################################

time = np.linspace(0, tTotal, nTotal)
opField = np.zeros(nFFT, dtype=complex)

############################
##  Iniciar Simulacion
############################

TFprom = 0
angProm = 0

currentTerm = eVinv * current(time)

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

    transFourier = np.fft.fft(opField)
    TFprom += (abs(np.fft.fftshift(transFourier)) *
               abs(np.fft.fftshift(transFourier))/float(nWindw))
    angProm += np.angle(np.fft.fftshift(transFourier))/float(nWindw)

#########################################
##  Representacion de los Datos
#########################################

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

fig, ax1 = plt.subplots()
ax1.plot(fftWL, angProm, "r")
ax1.set_xlabel("$\lambda$ [nm]", fontsize=15)
ax1.set_ylabel("FFT Phase [rad]", fontsize=15)
ax1.tick_params('y', colors='r')

ax2 = ax1.twinx()
ax2.plot(fftWL, TFprom, "b")
ax2.set_ylabel("PSD", fontsize=15)
ax2.set_yscale("log")
ax2.tick_params('y', colors='b')
fig.suptitle("$I_{Bias}$ = %i mA \t $V_{RF} = $ %.1f V" %(iBias, vRF*10**9),
                                                                fontsize = 20)
fig.tight_layout()
plt.show()
