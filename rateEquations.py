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

    CONSTANTES.................................................44

        Fisicas...........................................44
        Articulo..........................................55

            Tabla.....................................59
            Recopilado................................78
            Angel Valle...............................86
            Diego Chaves..............................97

        Muestreo..........................................115
        Simulacion........................................143

    Vectores..................................................215
    Simulacion................................................234
    Representacion............................................264

        Espectro Optico...................................273
        Frecuencia........................................299
        Both N y S........................................314

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

iBias = 35 *10**(-12) # bias current [C ns^-1]
fR = 5.0 #  [GHz]
vRF = 0 *10**(-9) #RMS voltage value of the signal generator [V]

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
tempIntev = 2.09268916658557#1.95791805783

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

TFprom = 0

currentTerm = eVinv * current(time)

for win in range(0, nWindw):
    for i in range(0, nTime-1):

        bTN = bTIntv * N[i] * N[i]

        invS = 1 / ((1/S[i]) + epsilon)

        N[i+1] = (N[i] + currentTerm[i] - aTIntv*N[i] - bTN -
                                    (cTIntv*N[i]**3) - vgT*N[i]*invS + vgtN*invS)

        S[i+1] = S[i] + vgTGmm*N[i]*invS - vgTGmmN*invS - intTtau*S[i] + btGmm*bTN

        Phi[i+1] = Phi[i] + aphvgTGmm*N[i] - faseConstant

        if ((i+1) % ndelta) == 0:

            index = int((i+1)/ndelta)
            opField[index] = np.sqrt(constP * S[i+1]) * np.exp(1j*Phi[i+1])

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

frecuencyLimits = 1 / (2*delta)
fftTime = np.linspace(-frecuencyLimits, frecuencyLimits, nFFT)

# freqTotal = frecuancia de emision del laser +
#                          + frecuencia de la fase (freq = 1/2pi dPhi/dt)
freqTotal = f0# + dfdT*tempIntev
fftTime += freqTotal

fftWL = c0/fftTime *10**(9) # longitud de onda [nm]

fig = plt.figure(figsize=(8,6))
plt.plot(fftWL, abs(TFprom))
plt.xlabel("$\lambda$ [nm]", fontsize=15)
plt.ylabel("PSD", fontsize=15)
plt.yscale("log")
plt.title("$V_{RF} = $"+str(vRF*10**9)+" V")
plt.show()
#fig.savefig("./Graficas/"+str(int(vRF*10**10))+"dV/EfftWL.png")

#-----------------------------------------------------------------------
# Variacion de la frecuencia en funcion del tiempo a partir de la fase
#
#                   1     d Phi
#       freq(t) = ------ -------
#                  2 pi    dt
#-----------------------------------------------------------------------
"""
freqTime = 1/(2*np.pi) *( (alpha / 2.0) * ((gamma * vg * dGdN * (N - nTr)) -
                                        (1/tauP)) + 2 * np.pi * dfdT * tempIntev )

fig = plt.figure(figsize=(8,6))
plt.plot(time, freqTime)
plt.xlabel("tiempo [ns]", fontsize=15)
plt.ylabel("$\\frac{d \Phi}{d t} [GHz]$", fontsize=15)
plt.title("Frecuancia a partir de la fase optica $\\nu \\propto \\frac{d\Phi}{dt}$")
plt.show()
fig.savefig("./Graficas/FrecPhi.png")
"""
#--------------------------------------------------------
#   Both carrier and photon density versus time
#--------------------------------------------------------
"""
fig, ax1 = plt.subplots()
ax1.plot(time, N, 'r', label="N(t)")
ax1.set_xlabel("tiempo [ns]", fontsize=20)
ax1.set_xlim(0, 1)
ax1.set_ylabel("N(t) [$m^3$]", color='r', fontsize=20)
ax1.tick_params('y', colors='r')

ax2 = ax1.twinx()
ax2.plot(time, S, 'b', label="S(t)")
ax2.set_ylabel("S(t) [$m^3$]", color='b', fontsize=20)
ax2.tick_params('y', colors='b')
fig.tight_layout()
plt.show()
fig.savefig("./Graficas/BothNetS.png")
"""
