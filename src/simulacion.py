"""
    Programa principal para 

"""
__author__ = 'Jaime Diez G-P'
__version__ = '1.0.0'
__email__ = "jaimediezgp@gmail.com"
__date__ = "May 27, 2019"

import numpy as np

from Constantes import *
from getTempValues import *

class Simulacion():

    def __init__(self, iBias, vRF, fR, numWindw=1):
        self.iBias = iBias
        self.vRF = vRF
        self.fR = fR
        self.numWindw = numWindw

        rInt = rIntLists[fR]

        deltaT = getDeltaT(self.iBias)
        self.faseTerm = faseConstant - pi2t * deltaT

        #                 2 sqrt(2) vRF 
        # I_bias + cLoss --------------- sin(2 pi fR t)
        #                   z0 + zL
        self.current = lambda t: (self.iBias*10**(-12)
                                  + (cLoss*2.0*np.sqrt(2)*self.vRF
                                  * np.sin(2*np.pi*fR*t)) / rInt)

        self.tWindw = 40.96
        self.tTrans = 2.2

    def transitorio(self):

        self.tWindw = 2.2
        self.tTrans = 0.2

        self.allSimulation()

    def rateEquations(self):

        self.tWindw = 3 / self.fR
        self.tTrans = 2.2

        self.allSimulation()

    def setArrays(self):
        """  Inicializar los vectores de tiempo (time), de la densidad """
        """     de portadores (N) densidad de fotones (S) y de la fase """
        """      optica (Phi)                                          """

        global mTrans, mTotal, nTotal, t

        mWindw = int(self.tWindw / delta)
        mTrans = int(self.tTrans / delta)
        tTotal = self.tWindw + self.tTrans
        nTotal = int(tTotal / tIntev)
        mTotal = int(tTotal / delta)

        t = np.linspace(0, tTotal, nTotal)
        self.time = np.linspace(0, tTotal, mTotal)
        self.I = np.zeros(mTotal)
        self.N = np.zeros(mTotal)
        self.S = np.zeros(mTotal)
        self.dPhi = np.zeros(mTotal)
        self.opField = np.zeros(mWindw, dtype=complex)
        self.currentTerm = eVinv * self.current(t)

    def allSimulation(self):
        self.setArrays()

        derivFaseTerm = self.faseTerm / tIntev
        self.TFprom = 0
        self.TFang = 0

        frecuencyLimits = 1 / (2*delta)
        self.fftFreq = np.linspace(-frecuencyLimits, frecuencyLimits,
                                   len(self.opField))#, endpoint=True)

        deltaF = getConstante(self.iBias)
        self.fftFreq += f0 - (deltaF/(2.0*np.pi))
        self.fftWL = (c0/self.fftFreq) *10**(9) # longitud de onda [nm]

        for win in range(0, self.numWindw):

            #---------------------------------------------------------
            # Vectores Gaussianos N(0,1) para el  Ruido
            #---------------------------------------------------------
            X = np.random.normal(0, 1, nTotal)
            Y = np.random.normal(0, 1, nTotal)

            # Se definen las condiciones iniciales para resolver la EDO
            tempN = nTr
            tempS = float(10**(15))
            tempPhi = 0

            for q in range(0, mTrans):
                for k in range(0, ndelta):

                    index = q*ndelta + k

                    bTN = bTIntv * tempN * tempN
                    invS = 1 / ((1/tempS) + epsilon)
                    sqrtS = np.sqrt(abs(tempS))

                    tempPhi = (tempPhi + aphvgTGmm*tempN - self.faseTerm
                               + ruidoPhi*tempN*Y[index]/sqrtS)

                    tempS = (tempS + vgTGmm*tempN*invS - vgTGmmN*invS
                             - intTtau*tempS + btGmm*bTN
                             + ruidoS*tempN*sqrtS*X[index])

                    tempN = (tempN + self.currentTerm[index] - aTIntv*tempN
                             - bTN - cTIntv*tempN**3 - vgT*tempN*invS + vgtN*invS)

                self.I[q] = self.currentTerm[index] *10**(12) / eVinv
                self.N[q] = tempN
                self.S[q] = tempS
                self.dPhi[q] = (1/(2*np.pi))*(derivAphvgTGmm*tempN
                                - derivFaseTerm + derivRuidoPhi*tempN*Y[q]/sqrtS)

            self.opField[0] = np.sqrt(constP * tempS) * np.exp(1j*tempPhi)

            for q in range(mTrans, mTotal):
                for k in range(0, ndelta):

                    index = q*ndelta + k + mTrans

                    bTN = bTIntv * tempN * tempN
                    invS = 1 / ((1/tempS) + epsilon)
                    sqrtS = np.sqrt(tempS)

                    tempPhi = (tempPhi + aphvgTGmm*tempN - self.faseTerm
                               + ruidoPhi*tempN*Y[index]/sqrtS)

                    tempS = (tempS + vgTGmm*tempN*invS - vgTGmmN*invS
                             - intTtau*tempS + btGmm*bTN
                             + ruidoS*tempN*sqrtS*X[index])

                    tempN = (tempN + self.currentTerm[index]
                             - aTIntv*tempN - bTN - (cTIntv*tempN**3)
                             - vgT*tempN*invS + vgtN*invS)

                self.I[q] = self.currentTerm[index]*10**(12) / eVinv
                self.N[q] = tempN
                self.S[q] = tempS
                self.dPhi[q] = ((derivAphvgTGmm*tempN - derivFaseTerm
                                 + derivRuidoPhi*tempN*Y[index]/np.sqrt(tempS))
                                 / (2*np.pi))
                self.opField[q-mTrans] = (np.sqrt(constP*tempS)
                                         * np.exp(1j*tempPhi))

            self.I[0] = 0
            self.N[0] = nTr
            self.S[0] = float(10**(15))
            transFourier = np.fft.fft(self.opField)
            self.TFprom += (abs(np.fft.fftshift(transFourier))
                            * abs(np.fft.fftshift(transFourier))
                            / float(self.numWindw))

            self.TFang += (np.angle(np.fft.fftshift(transFourier))
                           / float(self.numWindw))

    def save(self):
        nameRtEqtins = ("Data/RateEquations_%imA_%imV_%iGHZ"
                        %(self.iBias, self.vRF*10**(12), self.fR))
        np.savez(
                    nameRtEqtins, time=self.time,
                                  I=self.I,
                                  N=self.N,
                                  S=self.S,
                                  dPhi=self.dPhi
                )

        namePSD = ("Data/PSD_%imA_%imV_%iGHZ"
                   %(self.iBias, self.vRF*10**(12), self.fR))
        np.savez(
                    namePSD, fftWL=self.fftWL,
                             fftFreq=self.fftFreq,
                             TFprom=self.TFprom,
                             TFang=self.TFang
                )

if __name__ == '__main__':
    iBias = 35
    vRF = 1.0 *10**(-9)
    fR = 5.0
    numWindw = 1

    laser = Simulacion(iBias, vRF, fR, numWindw)
    laser.allSimulation()
    laser.save()
