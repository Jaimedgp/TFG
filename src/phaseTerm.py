import numpy as np
import matplotlib.pyplot as plt
import cmath

from Constantes import *


iBias = [i *10**(-12) for i in range(15, 70, 5) ] # bias current [C ns^-1]

############################
##  Iniciar Simulacion
############################

TotalC = []

for inten in iBias:

    current = eVinv*inten

    # Se definen las condiciones iniciales para resolver la EDO
    tempN = nTr
    tempS = float(10**(20))
    tempC = 0

    for q in range(0, nTotal):

        bTN = bTIntv * tempN * tempN

        invS = 1 / ((1/tempS) + epsilon)

        tempC = aphvgTGmm*tempN - faseConstant

        tempS = (tempS + vgTGmm*tempN*invS - vgTGmmN*invS - intTtau*tempS
                                                                 + btGmm*bTN)

        tempN = (tempN + current - aTIntv*tempN - bTN - (cTIntv*tempN**3) -
                                                    vgT*tempN*invS + vgtN*invS)

    TotalC.append(tempC)

curr=[]
exp=[]
with open("./Datos/tabla.dat", "r") as fr:
    line = fr.readline()
    while line[0] == "#":
        line = fr.readline()
    while line != "":
        current, experimental = line.split("\t")
        curr.append(current)
        exp.append(experimental)
        line = fr.readline()

with open("./Datos/Table.txt", "w") as fw:
    fw.write("""################################################################
##      Datos obtenidos de la simulacion (phaseTerm.py) en la
##  que se obtiene el valor del termino de la fase que no
##  depende ni de la diferencia de frecuencias ni de ruido y
##  que tiende a un valor constante.
##
##                  alpha  T                   1    T
##  Constante(C) = ------- | Gamma vg g(N) - ------ |
##                    2    L                  tauP  J
##
##      Los valores de la segunda columna son los obtenidos
##  a partir de la simulacion (phaseTerm.py) mientras que
##  los datos de las columnas 1 y 3 son de la tabla (tabla.dat)
##
##  Intensidad(mA)  Constante(ns^-1)    Delta T
################################################################
""")
    for i in range(len(TotalC)):
        fw.write(curr[i] + "\t" + str(TotalC[i]) + "\t" + exp[i])

with open("./src/getTempValues.py", 'w') as fw:
    fw.write("""def getDeltaT(intensidad='all'):
    diccionario = {""")

    for i in range(len(TotalC)-1):
        fw.write(curr[i] + " : " + exp[i] + ",\t\t\t\t\t")
    fw.write(curr[i] + " : " + exp[-1] + "}\n")
    fw.write("""    if intensidad == "all":
        return diccionario
    else:
        return diccionario[intensidad]\n""")

    fw.write("""def getConstante(intensidad='all'):
    diccionario = {""")

    for i in range(len(TotalC)-1):
        fw.write(curr[i] + " : " + str(TotalC[i]) + ",\n\t\t\t\t\t")
    fw.write(curr[i] + " : " + str(TotalC[-1]) + "}\n\n")
    fw.write("""    if intensidad == "all":
        return diccionario
    else:
        return diccionario[intensidad]\n""") 
