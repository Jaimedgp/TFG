import numpy as np

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
f0 = c0 / (1.546843 * 10**(-6)) # emission frequency at threshold [GHz]

# TERMINOS DE INYECCION
rSL = 0.1 # front facet reflectivity of the SL
kc = 42.3 # master-slave coupling coefficient [ns^-1]

#---------------------------------------------------
# Recopilado por el articulo
#---------------------------------------------------

fR = 5.0 # [GHz]
rInt = 142.25 # [Ohm]
#fR = 0.5 # [GHz]
#rInt = 103.36 # [Ohm]

#---------------------------------------------------
# Facilitados por Angel Valle
#---------------------------------------------------

ng = 3.5 # index of the 
vg = c0/ng #*10**(-9)# group velocity [m ns^-1]

cLoss = 1 # loss coeficient accounting for the frequency

# TERMINOS DE INYECCION
phiInyct = 0 # Fase de la inyeccion
sInyct = 3 *10**(21)# densidad de fotones de inyeccion
deltaNu = 5 # detunning of the injected laser field with respect to the frequency of the SL

################################################################################
##  Valores del muestreo para la simulacion
##
##      datos para el calculo de las ventanas de estudio, tiempos de
##      muestre para la FFT y tiempos de integracion
################################################################################

tWindw = 40.96 # tiempo de la ventana [ns]

tIntev = 1 *10**(-5) # tiempo de integracion [ns]
nTime = int(tWindw / tIntev) # numero de pasos de integracion

tTrans = 1.2 # tiempo del transitorio [ns]
nTrans = int(tTrans / tIntev)

tTotal = tWindw + tTrans
nTotal = int(tTotal / tIntev)

################################################################################
##  Constantes a user durante la simulacion
##
##      se computan antes para disminuir el numero de calculos
################################################################################

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

pi2t = np.pi * 2 * tIntev
# Fase constant
faseConstant = aphvgTGmmN + aphintTtau

#        h f0 Vact
# etaF -------------
#       Gamma tauP
constP = (etaF * h * f0 * vAct) / (gamma * tauP)

#---------------------------------------------------------
# Terminos de Ruido
#---------------------------------------------------------

#parte constante del termino de ruido del S(t)
ruidoS = np.sqrt(2 * beta * gamma * bTIntv)

# Parte constante del termino de ruido del Phi(t)
ruidoPhi = np.sqrt(beta * gamma * bTIntv / 2.0)

#---------------------------------------------------------
# Terminos de Ruido
#---------------------------------------------------------

ampltdS = 2 * kc * np.sqrt(sInyct) * tIntev

ampltdPhi = -kc * np.sqrt(sInyct) * tIntev

angConstnt = 2 * np.pi * deltaNu
