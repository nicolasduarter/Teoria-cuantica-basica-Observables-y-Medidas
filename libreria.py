import numpy as np
def amplitudTransicion(v1,v2):
    normaV1 = np.linalg.norm(v1)
    v1 = v1/normaV1
    normaV2 = np.linalg.norm(v2)
    v2 = v2/normaV2
    if len(v1) != len(v2) :
        return "las longitudes no coinsiden"
    else:
        rta = np.dot(np.conjugate(v2.T),v1)
        return rta

def media(matriz,vKet):
    Hermitania = np.array_equal(matriz,np.conjugate(np.transpose(matriz)))
    normaMatriz = np.linalg.norm(matriz)
    matriz = matriz/normaMatriz
    vKet = np.array([[1],[0]])
    normaV = np.linalg.norm(vKet)
    vKet = vKet/normaV
    if not Hermitania:
        return "La matriz no es hermitania"
    else:
        media = np.dot(np.conjugate(matriz.T),vKet)
        media = np.dot(np.conjugate(media.T),vKet)
        return media

def varianza(matriz,vKet):
    Hermitania = np.array_equal(matriz,np.conjugate(np.transpose(matriz)))
    normaMatriz = np.linalg.norm(matriz)
    matriz = matriz/normaMatriz
    vKet = np.array([[1],[0]])
    normaV = np.linalg.norm(vKet)
    vKet = vKet/normaV
    if not Hermitania:
        return "La matriz no es hermitania"
    else:
        longitud = len(vKet)
        matrizI = np.eye(int(longitud))
        media = media(matriz,vKet)
        op1 = matriz
        op2 = media* matrizI
        opDelta = op1 - op2
        op = np.dot(opDelta, opDelta)
        varianza = np.dot(np.conjugate(op.T),vKet)
        varianza = np.dot(np.conjugate(varianza.T),vKet)
        return varianza

def TransitoVectoresPropios(matriz,vKet):
    valPropios, vPropios = np.linalg.eig(matriz)
    probs=[]
    for rep in range (len(vPropios)):
        probabilidad = np.dot(np.conjugate(vPropios[rep].T),vKet)
        probabilidad = probabilidad**2
        probs.append(probabilidad)
        
    return probs,valPropios,vPropios

def eFinal(eInicial,matrices):
    numMatrices = len(matrices)
    matrices = map(np.array,matrices) 
    for rep in range(numMatrices):
        eInicial = np.dot(eInicial,matrices[rep])
        
    return eInicial

