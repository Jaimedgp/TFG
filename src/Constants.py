import numpy as np

#####################################
##      PHYSICAL CONSTANTS  (PDB)
#####################################

c0 = 0.299792458 # speed of light in vacuum [m ns^-1]
e = 1.6021766208 *10**(-19) # electron charge [C]
h = 6.626070040 *10**(-25) # Plank's constant [J ns]

################################################################################
##    Data taken from the scientific report
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

kc = 42.3 # master slave coupling coefficient [s^-1]
rSL = 0.1 # front facet reflectivity of the SL

etaF = 0.17 # in-fiber external quantum efficiency
f0 = c0 / (1.546843 * 10**(-6))# emission frequency at threshold [GHz]

#---------------------------------------------------
# Provided by Dr. Angel Valle
#---------------------------------------------------

# rInt = {fR [GHz]: rInt [Ohms]}
rIntLists = {5.0: 142.25, 0.5: 103.36} # [Ohm]

ng = 3.5 # index of the group
vg = c0/ng # group velocity [m ns^-1]

cLoss = 1 # loss coeficient accounting for the frequency

################################################################################
##  Sampling values for simulation
##
##      data for the calculation of the study windows, sample times
##      for the FFT and integration times
################################################################################

tIntev = 1 *10**(-5) # integration time [ns]

delta = 0.0025 # sample time for FFT [ns]
ndelta = int(delta / tIntev) # ndelta*tIntev=delta

################################################################################
##  Constants used during the simulation
##
##      some constants are computed before simultaion in order to 
##      decrease computation time
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

# Constant Phase
phaseConstant = aphvgTGmmN + aphintTtau

#        h f0 Vact
# etaF -------------
#       Gamma tauP
constP = (etaF * h * f0 * vAct) / (gamma * tauP)

#---------------------------------------------------------
# Noise Terms
#---------------------------------------------------------

# contant part of S(t) noise term
noiseS = np.sqrt(2 * beta * gamma * bTIntv)

# contant part of Phi(t) noise term
noisePhi = np.sqrt(beta * gamma * bTIntv / 2.0)

#---------------------------------------------------------
# CHIRP Term
#---------------------------------------------------------

derivAphvgTGmm = aphvgTGmm / tIntev
derivNoisePhi = noisePhi / tIntev
