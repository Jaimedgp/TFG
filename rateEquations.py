"""
    Simulation of the rate equations used to discribe the dynamics of the
    carrier density (N(t)), the photon density (S(t)) and optical phase (phi(t))

    The simulation is done without taking into account the noise and the
    inyection terms

    author: JaimeDGP
    latest version: 26 February 2019
"""

import numpy as np
import matplotlib.pyplot as plt

#####################################
##       CONSTANTES FISICAS
#####################################

c0 = 299792458 # speed of light in vacuum [m s^-1]
e = 1.6021766208 *10**(-19) # electron charge[C]

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
B = 1.5 *10**(-16) # spontaneous coefficient [m^3 s^-1]
dGdN = 4.38 *10**(-20) # differential gain [m^2]
tauP = 2.17 # photon lifetime [ps]
A = 2.8 *10**(8) # non-radiative coefficient [s^-1]
C = 9.0 *10**(-41) # Auger recombination coefficient [m^6 s^-1]
beta = 5.3 *10**(-6) #fraction of spontaneous emission coupled in2 lasing mode
epsilon = 1.97 *10**(-23) # non-linear gain coefficient [m^3]
alpha = 3 # linewidth engancement factor

#---------------------------------------------------
# Recopilado por el articulo
#---------------------------------------------------

iBias = 34 # bias current [mA]
fR = 5 #  [GHz]
vRF = 1.8 #RMS voltage value of the signal generator [V]

#---------------------------------------------------
# Facilitados por Angel Valle
#---------------------------------------------------

ng = 3.5 # index of the 
vg = c0/ng # group velocity [m s^-1

cLoss = 1 # loss coeficient accounting for the frequency
zL = 50 # impedance of the laser module [ohms]
z0 = 50 # generator output impedance [ohms]

#---------------------------------------------------
# Facilitados por las medidas de Diego Chaves
#---------------------------------------------------

dfdT = -13.5 # temperature coefficient of the emission frequency [GHz/K]
tempIntev = 1.3 # the difference between the active region temp at the
                # operating conditiona and at the threshold

################################################################################
##  Valores del muestre para la simulacion
##
##      datos tanto para el calculo de las ventanas de estudio, tiempos de
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

tIntev = 1 *10**(-9)

################################################################################
##  Constantes a user durante la simulacion
##
##      se computan antes para disminuir el numero de claculos
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

intTtau = tIntev * tauP

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
