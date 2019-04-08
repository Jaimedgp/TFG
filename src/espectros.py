import numpy as np
import matplotlib.pyplot as plt
import cmath

from Constantes import *
from getTempValues import *

frqMv = getConstante()
tmpMv = getDeltaT()


iBias = [i *10**(-12) for i in range(15, 40, 5) ] # bias current [C ns^-1]

nWindw = 1 # numero de ventanas (para promediar) N natural

delta = 0.0025 # tiempo de muestreo para la FFT [ns]
nFFT = int(tWindw / delta) # numero de puntos de la FFT (potencia de 2)
ndelta = int(delta / tIntev) # ndelta*tIntev=delta

nTrans = int(tTrans / delta)

############################
##  Iniciar Simulacion
############################
fig = plt.figure(figsize=(8,6))

WL = []

frecuencyLimits = 1 / (2*delta)
fftTime = np.linspace(-frecuencyLimits, frecuencyLimits, nFFT, endpoint=True)
fftTime += f0

for inten in iBias:
    TFprom = 0
    intensidad = int(inten *10**(12))
    opField = np.zeros(nFFT, dtype=complex)

    tmp = pi2t * tmpMv[intensidad]

    # Fase constant
    faseTerm = faseConstant - tmp

    currentTerm = eVinv * inten

    for win in range(nWindw):

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
            for k in range(0, ndelta):

                index = (q*ndelta) + k

                bTN = bTIntv * tempN * tempN

                invS = 1 / ((1/tempS) + epsilon)

                tempPhi = (tempPhi + aphvgTGmm*tempN - faseTerm +
                                    ruidoPhi*tempN*Y[index]/np.sqrt(abs(tempS)))

                tempS = (tempS + vgTGmm*tempN*invS - vgTGmmN*invS -
                                                    intTtau*tempS + btGmm*bTN +
                                    ruidoS*tempN*np.sqrt(abs(tempS))*X[index])

                tempN = (tempN + currentTerm - aTIntv*tempN - bTN -
                                (cTIntv*tempN**3) - vgT*tempN*invS + vgtN*invS)

        opField[0] = np.sqrt(constP * tempS) * np.exp(1j*tempPhi)

        usedIndex = nTrans*ndelta

        for q in range(1, nFFT):
            for k in range(0, ndelta):

                index = (q-1)*ndelta + k + usedIndex

                bTN = bTIntv * tempN * tempN

                invS = 1 / ((1/tempS) + epsilon)

                tempPhi = (tempPhi + aphvgTGmm*tempN - faseTerm +
                                        ruidoPhi*tempN*Y[index]/np.sqrt(tempS))

                tempS = (tempS + vgTGmm*tempN*invS - vgTGmmN*invS - intTtau*tempS
                            + btGmm*bTN + ruidoS*tempN*np.sqrt(tempS)*X[index])

                tempN = (tempN + currentTerm - aTIntv*tempN - bTN -
                                (cTIntv*tempN**3) - vgT*tempN*invS + vgtN*invS)

            opField[q] = np.sqrt(constP * tempS) * np.exp(1j*tempPhi)

        transFourier = np.fft.fft(opField)
        TFprom += abs(np.fft.fftshift(transFourier)/float(nWindw))

    fftWL = c0 *10**(9) / (fftTime - (frqMv[intensidad]/(2.0*np.pi))) # longitud de onda [nm]

    WL.append(fftWL[np.where(TFprom == max(TFprom))][0])

    plt.plot(fftWL, TFprom, label=str(intensidad))

plt.xlabel("$\lambda$ [nm]", fontsize=15)
plt.ylabel("PSD", fontsize=15)
plt.yscale("log")
plt.legend()
plt.show()
#fig.savefig("./Graficas/"+str(int(vRF*10**10))+"dV/EfftWL.png")

with open("./Datos/Landas.txt", 'w') as fw:
    fw.write("""######################################################################
##		Datos de la longitud de onda de emision del laser para
##	diferentes corrientes I_Bias obtenidos por simulacion con un
##	promedio a %s ventanas
##
##	Intensidad(mA)	Longitud de onda de pico(nm)
######################################################################\n""" %(nWindw))
    for i in range(len(iBias)):
        fw.write("%i \t %.2f \n" %(iBias[i]*10**12, WL[i]))
