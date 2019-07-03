
"""
    Programa para guardar los datos en un fichero txt para Angel Valle

"""

from simulation import Simulation
from Constants import *
from getDictValues import *

laser = Simulation(
                    iBias = 35, # bias current [mA]
                    vRF = 1.0*10**(-9), #RMS voltage value of the signal generator [V]
                    fR = 5.0, # [GHz]
                    pwrInjct = 100,
                    nuDetng = 5, # [GHz]
                    numWindw = 1
                  )

laser.allSimulation()

#--------------------------------------
#          primer fichero
#   primera tiempo
#   segunda corriente
#   tercera potencia
#   cuarta N(t)/N_tr.
#--------------------------------------
with open("./PrimerArchivo.txt", "w") as fw1:
    fw1.write("Tiempo"+"\t"+"Corriente"+"\t"+"Potencia"+"\t"+"N(t)/N_tr"+"\n")

    for i in range(len(laser.time)):
        fw1.write(str(laser.time[i])+"\t"
                 +str(laser.I[i])+"\t"
                 +str(laser.S[i]*constP)+"\t"
                 +str(laser.N[i]/nTr)+"\n")
    fw1.close()

#--------------------------------------
#          segundo fichero
#   primera tiempo
#   segunda chirping 
#   tercera fase
#--------------------------------------
with open("./SegundoArchivo.txt", "w") as fw2:
    fw2.write("Tiempo"+"\t"+"Chirping"+"\t"+"Fase"+"\n")

    for i in range(len(laser.time)):
        fw2.write(str(laser.time[i])+"\t"
                 +str((1/(2*np.pi))*laser.dPhi[i])+"\t"
                 +str(laser.Phi[i] - 2*np.pi*(laser.nuDetngEmssFrq - f0 +
                                           (c0/emissnWL[laser.iBias]))*laser.time[i])+"\n")
    fw2.close()

#--------------------------------------
#          Tercer fichero
#   primera lambda
#   segunda PSD
#--------------------------------------
with open("./TercerArchivo.txt", "w") as fw3:
    fw3.write("Lambda"+"\t"+"PSD"+"\n")

    for i in range(len(laser.TFavg)):
        fw3.write(str(laser.fftWL[i])+"\t"
                 +str(laser.TFavg[i])+"\n")
    fw3.close()
