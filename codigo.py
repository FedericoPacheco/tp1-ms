import csv

from matplotlib import pyplot
import numpy

#---------------------------------------------------------------------------------------------------------------------------------------------------------
# CSV
def leerSenial():
    senial = []
    # El archivo cardio.csv debe estar en el mismo
    # directorio en el que se encuentra este codigo
    with open(__file__ + '\..\cardio.csv', newline='') as archivoCSV:                                      
        lectorCSV = csv.reader(archivoCSV)
        for fila in lectorCSV:
            valor = float(fila[1])
            senial.append(valor)

    return senial

#---------------------------------------------------------------------------------------------------------------------------------------------------------
# Salida
def imprimir(senial):
    for i in range(len(senial)):
        print(senial[i])

def graficar(senial, tieneParteImaginaria, titulo, nombreAbscisa, nombreOrdenada):
    pyplot.figure(figsize=(12,8))
    pyplot.title(titulo)
    pyplot.xlabel(nombreAbscisa)
    pyplot.ylabel(nombreOrdenada)
    pyplot.axhline(0, color="black")
    pyplot.axvline(0, color="black")
    pyplot.grid()
    
    if (tieneParteImaginaria == True):
        pyplot.plot(numpy.real(senial), label = "Parte real", color = "blue", linewidth = 1) #, marker='o', markerfacecolor='blue', markersize=5)
        pyplot.plot(numpy.imag(senial), label = "Parte imaginaria", color = "red", linewidth = 1) #, marker='o', markerfacecolor='red', markersize=5)
        pyplot.legend()
    else:
        pyplot.plot(senial, color = 'blue', linewidth = 1) #, marker='o', markerfacecolor='blue', markersize=5)
        
    pyplot.show()

#---------------------------------------------------------------------------------------------------------------------------------------------------------
# Funciones auxiliares

coeficientesSavitzkyGolayDerivacion = [2/10, 1/10,0, -1/10, -2/10] 
coeficientesSavitzkyGolaySuavizado = [-21/231, 14/231, 39/231, 54/231, 59/231, 54/231, 39/231, 14/231, -21/231]         
                                    #[-3/35, 12/35, 17/35, 12/35, -3/35] 
                                    #[-2/21, 3/21, 6/21, 7/21, 6/21, 3/21, -2/21]

def lomo (inicioInclusivo, finInclusivo, largoSenial):
    lomo = []
    for i in range(0, largoSenial):
        lomo.append(0.0)

    for i in range(inicioInclusivo, finInclusivo + 1):
        lomo[i] = 1.0

    return lomo

def lomoInverso(inicioInclusivo, finInclusivo, largoSenial):
    lomoInverso = []
    for i in range(0, largoSenial):
        lomoInverso.append(1.0)

    lomoNormal = lomo(inicioInclusivo + 1, finInclusivo - 1, largoSenial)

    return restar(lomoInverso, lomoNormal)

#---------------------------------------------------------------------------------------------------------------------------------------------------------
# Operaciones basicas
def restar(senial1, senial2):
    senialResultado = []
    for i in range(len(senial1)): # Suponiendo len(senial1) = len(senial2)
        senialResultado.append(senial1[i] - senial2[i])

    return senialResultado

def multiplicar(senial1, senial2):
    senialResultado = []
    for i in range(len(senial1)): # Suponiendo len(senial1) = len(senial2)
        senialResultado.append(senial1[i] * senial2[i])

    return senialResultado

def sumar(senial1, senial2):
    senialResultado = []
    for i in range(len(senial1)): # Suponiendo len(senial1) = len(senial2)
        senialResultado.append(senial1[i] + senial2[i])

    return senialResultado

def filtrarPasabandaInverso(fourierSenial, piso, techo):
    return multiplicar(fourierSenial, lomoInverso(piso, techo, len(fourierSenial)))

#---------------------------------------------------------------------------------------------------------------------------------------------------------
# Inciso1
def acondicionar(senial): 
    senialAcondicionada = hacerHorizontal(senial)
    senialAcondicionada = eliminarOndasIndeseadas(senialAcondicionada)
    senialAcondicionada = suavizar(senialAcondicionada)

    return senialAcondicionada

def hacerHorizontal(senial):
    senialAcondicionada = []
    for i in range(len(senial)):
        # Restar la recta y = 54 + (7 / 400)*x
        senialAcondicionada.append(senial[i] - (senial[0] + (7 / 400) * i))
  
    return senialAcondicionada

def eliminarOndasIndeseadas(senial):
    # Eliminadas a mano / "a ojo"
    senialAcondicionada = multiplicar(senial, lomoInverso(100, 155, len(senial)))
    senialAcondicionada = multiplicar(senialAcondicionada, lomoInverso(250, 280, len(senial)))
    senialAcondicionada = multiplicar(senialAcondicionada, lomoInverso(375, len(senialAcondicionada), len(senial)))
    
    return senialAcondicionada

def suavizar(senial):
    return numpy.convolve(coeficientesSavitzkyGolaySuavizado, senial)

#---------------------------------------------------------------------------------------------------------------------------------------------------------
#Inciso2,3
def derivadas(senial): 
    graficar(derivar(senial, 1), False, "Aproximación de primera derivada", "n", "~x'[n]")
    graficar(derivar(senial, 2), False, "Aproximación de segunda derivada" , "n", "~x''[n]")

def derivar(senial, ordenDerivada):  #ordenDerivada > 0
    auxSenial = senial
    for i in range(ordenDerivada):
        auxSenial = numpy.convolve(coeficientesSavitzkyGolayDerivacion, auxSenial)

    return auxSenial

#---------------------------------------------------------------------------------------------------------------------------------------------------------
#Inciso4

#Filtros
pisoP = 25
techoP = 390

pisoQRS = 45
techoQRS = 375

pisoT = 25
techoT = 385

def reconstruccionPonSeparacionDeOndas(senial): 
    POriginal = obtenerOndaP(senial)
    QRSOriginal = obtenerOndaQRS(senial)
    TOriginal = obtenerOndaT(senial)

    fourierP = numpy.fft.fft(POriginal)
    fourierQRS = numpy.fft.fft(QRSOriginal)
    fourierT = numpy.fft.fft(TOriginal)

    fourierPFiltrada = filtrarPasabandaInverso(fourierP, pisoP, techoP)
    fourierQRSFiltrada = filtrarPasabandaInverso(fourierQRS, pisoQRS, techoQRS)
    fourierTFiltrada = filtrarPasabandaInverso(fourierT, pisoT, techoT)

    PAntitransformado = numpy.fft.ifft(fourierPFiltrada)
    QRSAntitransformado = numpy.fft.ifft(fourierQRSFiltrada)
    TAntitransformado = numpy.fft.ifft(fourierTFiltrada)

    graficar(POriginal, False, "Onda P original", "n", "p[n]")
    graficar(fourierP, True, "DFT onda P", "k", "P[k]")
    graficar(fourierPFiltrada, True, "DFT onda P filtrada", "k", "P[k]")
    graficar(PAntitransformado, False, "Onda P reconstruida con recorte de frecuencia", "n", "p[n]")

    graficar(QRSOriginal, False, "Onda QRS original", "n", "qrs[n]")
    graficar(fourierQRS, True, "DFT onda QRS", "k", "QRS[k]")
    graficar(fourierQRSFiltrada, True, "DFT onda QRS filtrada", "k", "QRS[k]")
    graficar(QRSAntitransformado, False, "Onda QRS reconstruida con recorte de frecuencia", "n", "qrs[n]")
    
    graficar(TOriginal, False, "Onda T original", "n", "t[n]") 
    graficar(fourierT, True, "DFT onda T", "k", "T[k]")
    graficar(fourierTFiltrada, True, "DFT onda T filtrada", "k", "T[k]")
    graficar(TAntitransformado, False, "Onda T reconstruida con recorte de frecuencia", "n", "t[n]")
    
    graficar(sumar(sumar(PAntitransformado, QRSAntitransformado), TAntitransformado), False, "Reconstrucción por separación de ondas", "n", "x[n]")

def obtenerOndaP(senial):
    lomo1 = lomo(0, 32, len(senial))
    lomo2 = lomo(150, 178, len(senial))
    lomo3 = lomo(280, 309, len(senial))
    trenDeLomos = sumar(sumar(lomo1, lomo2), lomo3)

    return multiplicar(senial, trenDeLomos)

def obtenerOndaQRS(senial):
    lomo1 = lomo(33, 50, len(senial))
    lomo2 = lomo(179, 200, len(senial))
    lomo3 = lomo(310, 330, len(senial))
    trenDeLomos = sumar(sumar(lomo1, lomo2), lomo3)

    return multiplicar(senial, trenDeLomos)

def obtenerOndaT(senial):
    lomo1 = lomo(51, 149, len(senial))
    lomo2 = lomo(201, 279, len(senial))
    lomo3 = lomo(331, len(senial) - 1, len(senial))
    trenDeLomos = sumar(sumar(lomo1, lomo2), lomo3)

    return multiplicar(senial, trenDeLomos)

#---------------------------------------------------------------------------------------------------------------------------------------------------------
#Inciso5
def reconstruccionConPulsoCompleto(senial): 
    fourierSenialCompleta = numpy.fft.fft(senial)
    fourierSenialCompletaFiltrada = filtrarPasabandaInverso(fourierSenialCompleta, min([pisoP, pisoQRS, pisoT]), max([techoP, techoQRS, techoT]))
    
    graficar(fourierSenialCompleta, True, "DFT pulso completo", "k", "X[k]")
    graficar(fourierSenialCompletaFiltrada, True, "DFT pulso completo filtrado", "k", "X[k]")
    graficar(numpy.fft.ifft(fourierSenialCompletaFiltrada), False, "Reconstrucción con pulso completo", "n", "x[n]")

#---------------------------------------------------------------------------------------------------------------------------------------------------------
if __name__ == "__main__":
    graficar(coeficientesSavitzkyGolaySuavizado, False, "Coeficientes Savitzky–Golay para suavizado", "n", "")
    graficar(coeficientesSavitzkyGolayDerivacion, False, "Coeficientes Savitzky–Golay para derivación", "n", "")

    senial = leerSenial()
    graficar(senial, False, "Señal original", "n", "x[n]")

    senial = acondicionar(senial)
    graficar(senial, False, "Señal acondicionada", "n", "x[n]")
    
    derivadas(senial)
    reconstruccionPonSeparacionDeOndas(senial)
    reconstruccionConPulsoCompleto(senial)
    