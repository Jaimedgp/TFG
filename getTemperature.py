import numpy as np
import matplotlib.pyplot as plt

def getTemperature(nameFile, intensity):
    with open(nameFile, "r") as fr:
        # read titles
        fr.readline()
        fr.readline()

        intensidad = []
        deltaT = []
        line = fr.readline()

        while line != '':
            values = line.split("\t")
            intensidad.append(float(values[0]))
            deltaT.append(float(values[2]))

            line = fr.readline()

    parameters = np.polyfit(intensidad, deltaT, len(deltaT)-1)

    result = 0
    xTh = np.linspace(15, 65, 100)
    for i in range(len(parameters)):
        result += parameters[len(parameters)-(i+1)]*xTh**(i)

    plt.plot(xTh, result)
    plt.plot(intensidad, deltaT, "ro")
    plt.show()
    result = 0
    for i in range(len(parameters)):
        result += parameters[len(parameters)-(i+1)]*intensity**(i)
    return result

s = getTemperature("DeltaT.txt", 34)
print s
