"""
    Simulation of the rate equations used to discribe the dynamics of the
    carrier density (N(t)), the photon density (S(t)) and optical phase (phi(t))

    The simulation is done without taking into account the noise and the
    inyection terms

    author: JaimeDGP
    latest version: 11 March 2019

                    ---------------
                    |    INDICE   |
                    ---------------

    CONSTANTES.................................................41

        Fisicas...........................................41
        Articulo..........................................52

            Tabla.....................................56
            Recopilado................................75
            Angel Valle...............................83
            Diego Chaves..............................94

        Muestreo..........................................112
        Simulacion........................................140

    Vectores..................................................212
    Simulacion................................................231
    Representacion............................................261

        Espectro Optico...................................272
"""

import numpy as np
import matplotlib.pyplot as plt
import cmath

#####################################
##       CONSTANTES FISICAS (PDB)
#####################################

c0 = 0.299792458 # speed of light in vacuum [m ns^-1]
e = 1.6021766208 *10**(-19) # electron charge [C]
h = 6.626070040 *10**(-25) # Plank's constant [J ns]

################################################################################
##    Datos tomados del articulo cientifico
##
##       NUMERICAL AND EXPERIMENTAL ANALYSIS OF OPTICAL FREQUECY COMB GENERATION
##                    IN GAIN-SWITCHED SEMICONDUCTOR LASERS
################################################################################

#---------------------------------------------------
# Table 1: Summary of main simulation parameters
#---------------------------------------------------

vAct = 1.53 *10**(-17) # active volume [m^3]
gamma = 0.06 # optical confinement factor
nTr = 1.3 *10**(24) # transparency carrier density [m^-3]
B = 1.5 *10**(-25) # spontaneous coefficient [m^3 ns^-1]
dGdN = 4.38 *10**(-20) # differential gain [m^2]
tauP = 2.17*10**(-3) # photon lifetime [ns]
A = 0.28 # non-radiative coefficient [ns^-1]
C = 9.0 *10**(-50) # Auger recombination coefficient [m^6 ns^-1]
beta = 5.3 *10**(-6) #fraction of spontaneous emission coupled in2 lasing mode
epsilon = 1.97 *10**(-23) # non-linear gain coefficient [m^3]
alpha = 3.0 # linewidth engancement factor

etaF = 0.17 # in-fiber external quantum efficiency
f0 = c0 / (1.546823 * 10**(-6))# emission frequency at threshold [GHz]

#---------------------------------------------------
# Recopilado por el articulo
#---------------------------------------------------

iBias = 34 *10**(-12) # bias current [C ns^-1]
fR = 5.0 #  [GHz]
vRF = 1 *10**(-9) #RMS voltage value of the signal generator [V]

#---------------------------------------------------
# Facilitados por Angel Valle
#---------------------------------------------------

ng = 3.5 # index of the 
vg = c0/ng #*10**(-9)# group velocity [m ns^-1]

cLoss = 1 # loss coeficient accounting for the frequency
zL = 50.0 # impedance of the laser module [ohms]
z0 = 50.0 # generator output impedance [ohms]

#---------------------------------------------------
# Facilitados por las medidas de Diego Chaves
#---------------------------------------------------

dfdT = -12.7 # temperature coefficient of the emission frequency [GHz/K]

"""
        the difference between the active region temp at the
    operating conditiona and at the threshold for Ibias = 34 mA

        the value has been obtain by a polynomial regression of a table of data
    (deltaT.txt) using a python script (getTemperature.py)
"""
tempIntev = 1.95791805783#2.09268916658557#

################################################################################
##  Valores del muestreo para la simulacion
##
##      datos para el calculo de las ventanas de estudio, tiempos de
##      muestre para la FFT y tiempos de integracion
################################################################################
"""
    Datos tomados del programa en Fortran de Angel Valle

nv = 20 # N de ventanas (para promediar) N natural
ttran  = 80 # Tiempo necesario para superar el transitorio(ns)
tventana = 40.96 # Tiempo de la ventana (ns)
tfinal = ttran + nv*tventana# Tiempo total que simulamos
delta = 0.25 #???d-2???? # Tiempo de muestreo para la FFT (ns)
tr = tfinal-ttran # Tiempo real que se utiliza para la FFT
no = int(tventana/delta) # N de valores de DFT (potencia de 2)
"""

nWindw = 2 # numero de ventanas (para promediar) N natural
tWindw = 40.96 / 2.0 # tiempo de la ventana [ns]

tIntev = 1 *10**(-5) # tiempo de integracion [ns]
nTime = int(tWindw / tIntev) # numero de pasos de integracion

delta = 0.0025 # tiempo de muestreo para la FFT [ns]
nFFT = int(tWindw / delta) # numero de puntos de la FFT (potencia de 2)
ndelta = int(delta / tIntev) # ndelta*tIntev=delta

tTrans = 1.2 # tiempo del transitorio [ns]
nTrans = int(tTrans / delta)

tTotal = tWindw + tTrans
nTotal = int(tTotal / tIntev)

################################################################################
##  Constantes a user durante la simulacion
##
##      se computan antes para disminuir el numero de calculos
################################################################################

#                  INTENSIDAD
#
#                 2 sqrt(2) vRF 
# I_bias + cLoss --------------- sin(2 pi fR t)
#                   z0 + zL
current = lambda t: (iBias + (cLoss * 2.0 * np.sqrt(2) * vRF * np.sin(2
                                                * np.pi * fR * t)) / (z0 + zL))
#   tIntev
# ----------
#   e Vact
eVinv = tIntev / (e*vAct)

aTIntv = A * tIntev
bTIntv = B * tIntev
cTIntv = C * tIntev

#       dg
#   vg ---- tIntev
#       dN
vgT = vg * dGdN * tIntev

#       dg
#   vg ---- tIntev Ntr
#       dN
vgtN = vgT * nTr

#       dg
#   vg ---- tIntev Gamma
#       dN
vgTGmm = vgT * gamma

#       dg
#   vg ---- tIntev Gamma Ntr
#       dN
vgTGmmN = vgTGmm * nTr

intTtau = tIntev / tauP

btGmm = beta * gamma

#    alpha      dg
#   ------- vg ---- tIntev Gamma
#      2        dN
aphvgTGmm = (alpha / 2.0) * vgTGmm

#    alpha      dg
#   ------- vg ---- tIntev Gamma Ntr
#      2        dN
aphvgTGmmN = aphvgTGmm * nTr

aphintTtau = (alpha / 2.0) * intTtau

#       TEMPERATURE TERM
#
#        df
#  2 pi ---- TempIntev tIntev
#        dT
tmp = 2 * np.pi *dfdT * tempIntev * tIntev

# Fase constant
faseConstant = aphvgTGmmN + aphintTtau - tmp

#        h f0 Vact
# etaF -------------
#       Gamma tauP
constP = (etaF * h * f0 * vAct) / (gamma * tauP)

#---------------------------------------------------------
# Terminos de Ruido
#---------------------------------------------------------

#parte constante del termino de ruido del S(t)
ruidoS = np.sqrt(2 * beta * gamma * bTIntv)

# Parte constante del termino de ruido del Phi
ruidoPhi = np.sqrt(beta * gamma * bTIntv / 2.0)

################################################################################
##  Inicializar los vectores de tiempo (time), de la densidad de portadores (N)
##      densidad de fotones (S) y de la fase optica (Phi)
################################################################################

time = np.linspace(0, tTotal, nTotal)
N = np.zeros(nFFT)
S = np.zeros(nFFT)
Phi = np.zeros(nFFT)

opField = np.zeros(nFFT, dtype=complex)

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
    tempS = float(10**(20))
    tempPhi = 0

    for q in range(0, nTrans):
        for k in range(0, ndelta):

            index = (q*ndelta) + k

            bTN = bTIntv * tempN * tempN

            invS = 1 / ((1/tempS) + epsilon)

            tempPhi = (tempPhi + aphvgTGmm*tempN - faseConstant +
                                            ruidoPhi*tempN*Y[index]/np.sqrt(abs(tempS)))

            tempS = (tempS + vgTGmm*tempN*invS - vgTGmmN*invS - intTtau*tempS
                                + btGmm*bTN +
                     ruidoS*tempN*np.sqrt(abs(tempS))*X[index])

            tempN = (tempN + currentTerm[index] - aTIntv*tempN - bTN -
                                                                (cTIntv*tempN**3) -
                                                        vgT*tempN*invS + vgtN*invS)

    N[0] = tempN
    S[0] = tempS
    Phi[0] = tempPhi
    opField[0] = np.sqrt(constP * tempS)

    usedIndex = nTrans*ndelta

    for q in range(1, nFFT):
        for k in range(0, ndelta):

            index = (q-1)*ndelta + k + usedIndex

            bTN = bTIntv * tempN * tempN

            invS = 1 / ((1/tempS) + epsilon)

            tempPhi = (tempPhi + aphvgTGmm*tempN - faseConstant +
                                                ruidoPhi*tempN*Y[index]/np.sqrt(tempS))

            tempS = (tempS + vgTGmm*tempN*invS - vgTGmmN*invS - intTtau*tempS
                                + btGmm*bTN +
                     ruidoS*tempN*np.sqrt(tempS)*X[index])

            tempN = (tempN + currentTerm[index] - aTIntv*tempN - bTN -
                                                                (cTIntv*tempN**3) -
                                                        vgT*tempN*invS + vgtN*invS)

        N[q] += tempN/float(nWindw)
        S[q] += tempS/float(nWindw)
        Phi[q] += tempPhi/float(nWindw)

        opField[q] = np.sqrt(constP * tempS) * np.exp(1j*tempPhi)

    transFourier = np.fft.fft(opField)
    TFprom += np.fft.fftshift(transFourier)/float(nWindw)

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

plt.plot(S)
plt.show()

frecuencyLimits = 1 / (2*delta)
fftTime = np.linspace(-frecuencyLimits, frecuencyLimits, nFFT)

# freqTotal = frecuancia de emision del laser +
#                          + frecuencia de la fase (freq = 1/2pi dPhi/dt)
#freqTotal = f0# + dfdT*tempIntev
fftTime += f0

fftWL = c0/fftTime *10**(9) # longitud de onda [nm]

fig = plt.figure(figsize=(8,6))
plt.plot(fftWL, abs(TFprom))
plt.xlabel("$\lambda$ [nm]", fontsize=15)
plt.ylabel("PSD", fontsize=15)
plt.yscale("log")
plt.title("$I_{Bias}$ = "+str(iBias*10**12)+" mA \t $V_{RF} = $"+str(vRF*10**9)+" V")
plt.show()
#fig.savefig("./Graficas/"+str(int(vRF*10**10))+"dV/EfftWL.png")

