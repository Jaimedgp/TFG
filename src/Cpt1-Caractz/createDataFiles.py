"""
    Programa principal para 

"""
__author__ = 'Jaime Diez G-P'
__version__ = '1.0.0'
__email__ = "jaimediezgp@gmail.com"
__date__ = "May 27, 2019"

from simulacion import Simulacion

def espectrosData():
    iBias = [i for i in range(15, 40)]
    vRf = 0
    fR = 5.0
    numWindw = 10

    for intBias in iBias:
        laser = Simulacion(intBias, vRF, fR, numWindw)
        laser.allSimulation()
        laser.save()

def psdData():
    iBias = 30
    vRf = [0.05*10**(-9), 1.0*10**(-9), 1.5*10**(-9)]
    fR = 5.0
    numWindw = 10

    for voltRF in vRF:
        laser = Simulacion(intBias, voltRF, fR, numWindw)
        laser.allSimulation()
        laser.save()

def frequencyData():
    iBias = 50
    vRf = [0.05*10**(-9), 0.4*10**(-9), 1.0*10**(-9), 1.2*10**(-9)]
    fR = 0.5
    numWindw = 10

    for voltRF in vRF:
        laser = Simulacion(intBias, voltRF, fR, numWindw)
        laser.allSimulation()
        laser.save()

def currentData():
    iBias = [30, 50]
    vRf = 1.0 *10**(-9)
    fR = 5.0
    numWindw = 10

    for intBias in iBias:
        laser = Simulacion(intBias, vRF, fR, numWindw)
        laser.allSimulation()
        laser.save()

if __name__ == '__main__':

    espectrosData()
    psdData()
    frequencyData()
    currentData()
