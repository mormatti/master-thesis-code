from matplotlib import pyplot as plt
import scipy.sparse as ss
import scipy.sparse.linalg as ssl
import scipy.linalg as sclin
from bandStates import *
from minimizationTools import optimizedThetas

print("L = ", L)

# print(matrixElement(bandState(0), localHamiltonian(L,4), bandState(0)))

H = totalHamiltonian(L)

T = translationOperator(L)

R = xReflectionOperator(L)

F = yReflectionOperator(L)

reciprocalSpaceSites = [int(x - (L - 1)/2) for x in range(0,L)]

# print(reciprocalSpaceSites)

# print("Sqrt of -i = ", np.sqrt(complex(0,-1)))

def symmetricBandState(k):
    bk = bandState(k)
    bmk = bandState(-k)
    return bk / np.sqrt(matrixElement(bmk, R, bk))

def rSquaredHamiltonian():
    jCentered = - (L - 1)/2
    result = localHamiltonian(L,0) * jCentered**2
    for j in range(1,L):
        jCentered = j - (L - 1)/2
        result += localHamiltonian(L,j) * jCentered**2
    return result

varH = rSquaredHamiltonian()

MatrixM = [[matrixElement(symmetricBandState(x), varH, symmetricBandState(y)) for x in reciprocalSpaceSites] for y in reciprocalSpaceSites]

OT = optimizedThetas(MatrixM)
result = OT.x
functionResult = OT.fun
resultUnivoque = (result - result[0]) % (2 * np.pi)
phases = np.exp(complex(0,1) * resultUnivoque)

# print(resultUnivoque)
# print(phases)

def localizedWannier():
    W = symmetricBandState(0) * 0
    for i in reciprocalSpaceSites:
        W += symmetricBandState(i) * phases[int(i + (L - 1)/2)]
    return W / np.sqrt(L)

locW = localizedWannier()

print("")

printMatrix(H)

locW = np.matmul(T,locW)
locW = np.matmul(T,locW)
locW = np.matmul(T,locW)

Sites = [j for j in range(0,L)]
Energies = [expectObsValue(localHamiltonian(L,j), locW) for j in range(0,L)]

# plt.plot(Sites, Energies, 'bo-')
# plt.show()

print("")