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

iBias = [i *10**(-12) for i in range(15, 70, 5) ] # bias current [C ns^-1]

#---------------------------------------------------
# Facilitados por Angel Valle
#---------------------------------------------------

ng = 3.5 # index of the 
vg = c0/ng #*10**(-9)# group velocity [m ns^-1]

################################################################################
##  Valores del muestreo para la simulacion
##
##      datos para el calculo de las ventanas de estudio, tiempos de
##      muestre para la FFT y tiempos de integracion
################################################################################

nWindw = 1 # numero de ventanas (para promediar) N natural
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
eVinvTerm = tIntev / (e*vAct)

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

# Fase constant
faseConstant = aphvgTGmmN + aphintTtau

#---------------------------------------------------------
# Terminos de Ruido
#---------------------------------------------------------

#parte constante del termino de ruido del S(t)
ruidoS = np.sqrt(2 * beta * gamma * bTIntv)

# Parte constante del termino de ruido del Phi
#ruidoPhi = np.sqrt(beta * gamma * bTIntv / 2.0)

################################################################################
##  Inicializar los vectores de tiempo (time), de la densidad de portadores (N)
##      densidad de fotones (S) y de la fase optica (Phi)
################################################################################

TotalC = []
############################
##  Iniciar Simulacion
############################

for inten in iBias:

    eVinv = eVinvTerm*inten
    N = np.zeros(nTotal)
    S = np.zeros(nTotal)
    C = np.zeros(nTotal)

    #---------------------------------------------------------
    # Vectores Gaussianos N(0,1) para el  Ruido
    #---------------------------------------------------------

    X = np.random.normal(0, 1, nTotal)
    Y = np.random.normal(0, 1, nTotal)

    # Se definen las condiciones iniciales para resolver la EDO
    tempN = nTr
    tempS = float(10**(20))
    tempC = 0

    for q in range(0, nTrans):

        bTN = bTIntv * tempN * tempN

        invS = 1 / ((1/tempS) + epsilon)

        tempC = aphvgTGmm*tempN - faseConstant

        tempS = (tempS + vgTGmm*tempN*invS - vgTGmmN*invS - intTtau*tempS
                            + btGmm*bTN)# +
                   # ruidoS*tempN*np.sqrt(abs(tempS))*X[q])

        tempN = (tempN + eVinv - aTIntv*tempN - bTN - (cTIntv*tempN**3) -
                                                vgT*tempN*invS + vgtN*invS)

        """
        N[q] = tempN
        S[q] = tempS
        C[q] = tempC
        """

    for q in range(nTrans, nTotal):

        bTN = bTIntv * tempN * tempN

        invS = 1 / ((1/tempS) + epsilon)

        tempC = aphvgTGmm*tempN - faseConstant

        tempS = (tempS + vgTGmm*tempN*invS - vgTGmmN*invS - intTtau*tempS
                            + btGmm*bTN )#+
                   # ruidoS*tempN*np.sqrt(tempS)*X[q])

        tempN = (tempN + eVinv - aTIntv*tempN - bTN -
                                                            (cTIntv*tempN**3) -
                                                    vgT*tempN*invS + vgtN*invS)

        """
        N[q] = tempN
        S[q] = tempS
        C[q] = tempC
        """

    TotalC.append(tempC)#np.mean(C[nTrans:]))

curr=[]
exp=[]
with open("./tabla.dat", "r") as fr:
    line = fr.readline()
    while line != "":
        current, experimental = line.split("\t")
        curr.append(current)
        exp.append(experimental)
        line = fr.readline()

with open("./Table.txt", "w") as fw:
    for i in range(len(TotalC)):
        fw.write(str(curr[i]) + "\t" + str(TotalC[i]) + "\t" + str(exp[i]))
