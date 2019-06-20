"""
        Main program for the simulation of semiconductor laser by the
    resolution of the sochastic diferential equations. The dynamics of a
    semiconductor laser will be simulated in a discrete way subjected to a
    modulation of high frequency current and to the injection of light from an
    external laser. The optical frequency combs generated in the previous
    system are light sources with a spectrum formed by lines regularly spaced
    in frequency. These combs have applications in optical communications and
    spectroscopy. The work will consist of writing the laser simulation code
    including current modulation, intrinsic system noise and light injection
    from an external laser.

"""
__author__ = 'Jaime Diez G-P'
__version__ = '1.0.0'
__email__ = "jaimediezgp@gmail.com"
__date__ = "Jun 13, 2019"

import numpy as np

from Constants import *
from getTempValues import *

class Simulation():

    def __init__(self, iBias, vRF, fR, sInjct=0, nuDetng=0, numWindw=1):
        self.iBias = iBias
        self.vRF = vRF
        self.fR = fR
        self.sInjct = sInjct
        self.nuDetngEmssFrq = nuDetng
        # detuning of the injected laser field with respect to the frequency of the
        # slaver laser (SL) at threshold
        # = detuning respect to emission freq - emission freq + freq at threshold
        self.nuDetng = nuDetng - f0 + (c0/(1.54705*10**(-6))) # nu - nuTH + nuI
        self.numWindw = numWindw

        rInt = rIntLists[fR]

        deltaT = getDeltaT(self.iBias)
        self.phaseTerm = phaseConstant - pi2t * deltaT

        #                 2 sqrt(2) vRF 
        # I_bias + cLoss --------------- sin(2 pi fR t)
        #                   z0 + zL
        self.current = lambda t: (self.iBias*10**(-12)
                                  + (cLoss*2.0*np.sqrt(2)*self.vRF
                                  * np.sin(2*np.pi*fR*t)) / rInt
                                 )

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
        """ Initialize values and array used in simulation to store    """
        """     data as time (time), carrier density (N), photon       """
        """       density (S), derivative optical phase (dPhi),        """
        """            and total optical field (opField)               """

        global mTrans, mTotal, nTotal, t, angInject, senInject, cosInject

        mWindw = int(self.tWindw / delta)
        mTrans = int(self.tTrans / delta)
        tTotal = self.tWindw + self.tTrans
        mTotal = int(tTotal / delta)
        nTotal = mTotal*ndelta #int(tTotal / tIntev)

        t = np.linspace(0, tTotal, nTotal)
        angInject = 2*np.pi*self.nuDetng*t
        senInject = np.sin(angInject)
        cosInject = np.cos(angInject)
        self.time = np.linspace(0, tTotal, mTotal)
        self.I = np.zeros(mTotal)
        self.N = np.zeros(mTotal)
        self.S = np.zeros(mTotal)
        self.dPhi = np.zeros(mTotal)
        self.opField = np.zeros(mWindw, dtype=complex)
        self.currentTerm = eVinv * self.current(t)

    def allSimulation(self):
        self.setArrays()

        ampInject = kc*np.sqrt(self.sInjct)*tIntev
        opFldInject = np.sqrt(rSL*constP*self.sInjct)*np.exp(1j*np.pi)

        derivPhaseTerm = self.phaseTerm / tIntev
        self.TFavg = 0
        self.TFang = 0

        frecuencyLimits = 1 / (2*delta)
        self.fftFreq = np.linspace(-frecuencyLimits, frecuencyLimits,
                                   len(self.opField)
                                  )

        deltaF = getConstants(self.iBias)
        self.fftFreq += f0 - (deltaF/(2.0*np.pi))
        self.fftWL = (c0/self.fftFreq) *10**(9) # Wavelength [nm]

        for win in range(0, self.numWindw):

            #---------------------------------------------------------
            # Gaussian arrays N(0,1) for the Noise
            #---------------------------------------------------------
            X = np.random.normal(0, 1, nTotal)
            Y = np.random.normal(0, 1, nTotal)

            # Initial conditions are defined in order to resolved the SDE
            tempN = nTr
            tempS = float(10**(15))
            tempPhi = 0

            for q in range(0, mTrans):
                for k in range(0, ndelta):

                    index = q*ndelta + k

                    bTN = bTIntv * tempN * tempN
                    invS = 1 / ((1/tempS) + epsilon)
                    sqrtS = np.sqrt(abs(tempS))
                    cosPhi = np.cos(tempPhi)
                    senPhi = np.sin(tempPhi)

                    tempPhi = (tempPhi + aphvgTGmm*tempN - self.phaseTerm
                               + noisePhi*tempN*Y[index]/sqrtS
                               - (ampInject/sqrtS)*senPhi*cosInject[index]
                               + (ampInject/sqrtS)*cosPhi*senInject[index]
                              )

                    tempS = (tempS + vgTGmm*tempN*invS - vgTGmmN*invS
                             - intTtau*tempS + btGmm*bTN
                             + noiseS*tempN*sqrtS*X[index]
                             + 2*ampInject*sqrtS*cosPhi*cosInject[index]
                             + 2*ampInject*sqrtS*senPhi*senInject[index]
                            )

                    tempN = (tempN + self.currentTerm[index] - aTIntv*tempN
                             - bTN - cTIntv*tempN**3 - vgT*tempN*invS + vgtN*invS
                            )

                self.I[q] = self.currentTerm[index] *10**(12) / eVinv
                self.N[q] = tempN
                self.S[q] = tempS
                self.dPhi[q] = (1/(2*np.pi))*(derivAphvgTGmm*tempN
                                - derivPhaseTerm + derivNoisePhi*tempN*Y[q]/sqrtS
                               )

            for q in range(mTrans, mTotal):
                for k in range(0, ndelta):

                    index = q*ndelta + k

                    bTN = bTIntv * tempN * tempN
                    invS = 1 / ((1/tempS) + epsilon)
                    sqrtS = np.sqrt(tempS)
                    cosPhi = np.cos(tempPhi)
                    senPhi = np.sin(tempPhi)

                    tempPhi = (tempPhi + aphvgTGmm*tempN - self.phaseTerm
                               + noisePhi*tempN*Y[index]/sqrtS
                               - (ampInject/sqrtS)*senPhi*cosInject[index]
                               + (ampInject/sqrtS)*cosPhi*senInject[index]
                              )

                    tempS = (tempS + vgTGmm*tempN*invS - vgTGmmN*invS
                             - intTtau*tempS + btGmm*bTN
                             + noiseS*tempN*sqrtS*X[index]
                             + 2*ampInject*sqrtS*cosPhi*cosInject[index]
                             + 2*ampInject*sqrtS*senPhi*senInject[index]
                            )

                    tempN = (tempN + self.currentTerm[index]
                             - aTIntv*tempN - bTN - (cTIntv*tempN**3)
                             - vgT*tempN*invS + vgtN*invS
                            )

                self.I[q] = self.currentTerm[index]*10**(12) / eVinv
                self.N[q] = tempN
                self.S[q] = tempS
                self.dPhi[q] = ((derivAphvgTGmm*tempN - derivPhaseTerm
                                 + derivNoisePhi*tempN*Y[index]/np.sqrt(tempS))
                                 / (2*np.pi)
                               )
                self.opField[q-mTrans] = (np.sqrt(constP*tempS)
                                          * np.exp(1j*tempPhi)
                                          + opFldInject
                                          * np.exp(1j*angInject[index])
                                         )

            self.I[0] = 0
            self.N[0] = nTr
            self.S[0] = float(10**(15))
            transFourier = np.fft.fft(self.opField)
            self.TFavg += (abs(np.fft.fftshift(transFourier))
                            * abs(np.fft.fftshift(transFourier))
                            / float(self.numWindw)
                           )

            self.TFang += (np.angle(np.fft.fftshift(transFourier))
                           / float(self.numWindw)
                          )

    def save(self):
        if self.sInjct == 0.0:
            laserCharactz = ("_%imA_%imV_%iGHZ"
                             %(self.iBias, self.vRF*10**(12), self.fR)
                            )
        else:
            oderMg = np.log10(self.sInjct / 4)
            laserCharactz = ("_%iS_%iGHz"
                             %(oderMg, self.nuDetngEmssFrq)
                            )

        nameRtEqtins = "Data/RateEquations"+laserCharactz
        np.savez(
                    nameRtEqtins, time=self.time,
                                  I=self.I,
                                  N=self.N,
                                  S=self.S,
                                  dPhi=self.dPhi
                )

        namePSD = "Data/PSD"+laserCharactz
        np.savez(
                    namePSD, fftWL=self.fftWL,
                             fftFreq=self.fftFreq,
                             TFavg=self.TFavg,
                             TFang=self.TFang
                )

if __name__ == '__main__':
    iBias = 35
    vRF = 1.0 *10**(-9)
    fR = 5.0
    numWindw = 1

    laser = Simulation(iBias, vRF, fR, numWindw)
    laser.allSimulation()
    #laser.save()
