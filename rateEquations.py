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

    CONSTANTES.................................................40

        Fisicas...........................................40
        Articulo..........................................51

            Tabla.....................................55
            Recopilado................................74
            Angel Valle...............................82
            Diego Chaves..............................93

        Muestreo..........................................111
        Simulacion........................................140

    Vectores..................................................212
    Simulacion................................................231
    Representacion............................................260

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
alpha = 3 # linewidth engancement factor

etaF = 0.17 # in-fiber external quantum efficiency
f0 = c0 / (1.546823 * 10**(-6))# emission frequency [GHz]

#---------------------------------------------------
# Recopilado por el articulo
#---------------------------------------------------

iBias = 34 *10**(-12) # bias current [C ns^-1]
fR = 5 #  [GHz]
vRF = 1.0 *10**(-9) #RMS voltage value of the signal generator [V]

#---------------------------------------------------
# Facilitados por Angel Valle
#---------------------------------------------------

ng = 3.5 # index of the 
vg = c0/ng #*10**(-9)# group velocity [m ns^-1]

cLoss = 1 # loss coeficient accounting for the frequency
zL = 50 # impedance of the laser module [ohms]
z0 = 50 # generator output impedance [ohms]

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
tempIntev = 1.9579185783

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

nWindw = 1 # numero de ventanas (para promediar) N natural
tWindw = 40.96 # tiempo de la ventana [ns]
tFinal = nWindw * tWindw # tiempo total simulado

tIntev = 1 *10**(-5) # tiempo de integracion [ns]
nTime = int(tWindw / tIntev) # numero de pasos de integracion

delta = 0.0025 # tiempo de muestreo para la FFT [ns]
nFFT = int(tWindw / delta) # numero de puntos de la FFT (potencia de 2)
ndelta = 250.0 # ndelta*tIntev=delta

################################################################################
##  Constantes a user durante la simulacion
##
##      se computan antes para disminuir el numero de calculos
################################################################################

#   Delta_t I_bias
# -----------------
#      e Vact
tIeV = (tIntev * iBias) / (e*vAct)

#   tIntev cLoss 2 sqrt(2) vRF
# ------------------------------
#       e vAct (z0 + zL)
amplit = (tIntev * cLoss * 2 * np.sqrt(2) * vRF) / (e * vAct * (z0 + zL))

angFreq = 2 * np.pi * fR

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
aphvgTGmm = (alpha / 2) * vgTGmm

#    alpha      dg
#   ------- vg ---- tIntev Gamma Ntr
#      2        dN
aphvgTGmmN = aphvgTGmm * nTr

aphintTtau = (alpha / 2) * intTtau

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

################################################################################
##  Inicializar los vectores de tiempo (time), de la densidad de portadores (N)
##      densidad de fotones (S) y de la fase optica (Phi)
################################################################################

time = np.linspace(0, tFinal, nTime)
N = np.zeros(nTime)
S = np.zeros(nTime)
Phi = np.zeros(nTime)

opField = np.zeros(nFFT, dtype=complex)
topField = np.linspace(0, tFinal, nFFT)

# Se definen las condiciones iniciales para resolver la EDO
N[0] = nTr
S[0] = 10**(20)
Phi[0] = 0

opField[0] = np.sqrt(constP * S[0])

############################
##  Iniciar Simulacion
############################

sint = np.sin(angFreq*time)
TFprom = 0

for win in range(0, nWindw):
    for i in range(0, nTime-1):

        bTN = bTIntv * N[i] * N[i]

        invS = 1 / ((1/S[i]) + epsilon)

        N[i+1] = (N[i] + tIeV + amplit*sint[i] - aTIntv*N[i] - bTN -
                                    (cTIntv*N[i]**3) - vgT*N[i]*invS + vgtN*invS)

        S[i+1] = S[i] + vgTGmm*N[i]*invS - vgTGmmN*invS - intTtau*S[i] + btGmm*bTN

        Phi[i+1] = Phi[i] + aphvgTGmm*N[i] - faseConstant

        if ((i+1) % ndelta) == 0:

            index = int((i+1)/ndelta)
            opField[index] = np.sqrt(constP * S[i+1]) * np.exp(1j*Phi[i+1])

    transFourier = np.fft.fft(opField)
    TFprom += np.fft.fftshift(transFourier)/nWindw

#########################################
##  Representacion de los Datos
#########################################

fftTime = np.fft.fftfreq(nFFT, d=1/ndelta)
fftTime = np.fft.fftshift(fftTime)

fig = plt.figure(figsize=(8,6))
plt.plot(fftTime, abs(TFprom))
plt.xlabel("$\\nu$ [GHz]", fontsize=15)
plt.ylabel("EOP", fontsize=15)
plt.title("Transformada de Fourier de $E(t)$")
plt.show()
fig.savefig("./Graficas/Efft.png")
