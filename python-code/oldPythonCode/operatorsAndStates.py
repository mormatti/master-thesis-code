from configurationsNoFermions import *

# OPERATORS
# An operator is a complex matrix, which is a numpy array.

def printMatrix(L: list[list[float]]):
    print('{', end='')
    for i in range(len(L)):
        print('{', end='')
        for j in range(len(L[i])):
            print(str(L[i][j]), end='')
            if (j+1 != len(L[i])):
                print(',', end='')
        print('}', end='')
        if (i+1 != len(L)):
            print(',', end='')
    print('}')

def electricEnergy(configuration: list[list[int]]):
    totalEnergy = 0
    L = len(configuration)
    for i in range(L):
        totalEnergy += (plaquetteDoubleSpin("L", configuration, i) / 2.0)**2 \
                        + (plaquetteDoubleSpin("B", configuration, i) / 2.0)**2 \
                        + (plaquetteDoubleSpin("T", configuration, i) / 2.0)**2
    return np.array(totalEnergy)

def localElectricEnergy(j: int, configuration: list[list[int]]):
    localEnergy = 0
    localEnergy += ((plaquetteDoubleSpin("L", configuration, j) / 2.0)**2) / 2.0 \
                    + (plaquetteDoubleSpin("B", configuration, j) / 2.0)**2 \
                    + (plaquetteDoubleSpin("T", configuration, j) / 2.0)**2 \
                    + ((plaquetteDoubleSpin("R", configuration, j) / 2.0)**2) / 2.0
    return np.array(localEnergy)

def translationOperator(N: int, dsh: int = dshGlobal, dsv: int = dsvGlobal):
    CC = chainConfigurations(N,dsh,dsv)
    N = len(CC)
    T = [[0 for x in range(N)] for y in range(N)]
    for i in range(N):
        T[i][CC.index(translated(CC[i]))] = 1
    return np.array(T)

def fieldReflectionOperator(N: int, dsh: int = dshGlobal, dsv: int = dsvGlobal):
    CC = chainConfigurations(N,dsh,dsv)
    N = len(CC)
    S = [[0 for x in range(N)] for y in range(N)]
    for i in range(N):
        S[i][CC.index(invertedField(CC[i]))] = 1
    return np.array(S)

def xReflectionOperator(N: int, dsh: int = dshGlobal, dsv: int = dsvGlobal):
    CC = chainConfigurations(N,dsh,dsv)
    N = len(CC)
    S = [[0 for x in range(N)] for y in range(N)]
    for i in range(N):
        S[i][CC.index(xReflected(CC[i]))] = 1
    return np.array(S)

def yReflectionOperator(N: int, dsh: int = dshGlobal, dsv: int = dsvGlobal):
    CC = chainConfigurations(N,dsh,dsv)
    N = len(CC)
    F = [[0 for x in range(N)] for y in range(N)]
    for i in range(N):
        F[i][CC.index(yReflected(CC[i]))] = 1
    return np.array(F)

def totalHamiltonian(N: int, dsh: int = dshGlobal, dsv: int = dsvGlobal):
    CC = chainConfigurations(N,dsh,dsv)
    L = len(CC)
    M = [[0 for x in range(L)] for y in range(L)]
    for c in range(L):
        conf = CC[c]
        for p in range(N):
            newConf1 = applyMagnetic(conf,p,1)
            newConf2 = applyMagnetic(conf,p,-1)
            if (coherentInSpin(newConf1,p,dsh,dsv)):
                M[c][CC.index(newConf1)] += magneticPrefactor(conf,p,1,dsh,dsv) * gB
            if (coherentInSpin(newConf2,p,dsh,dsv)):
                M[c][CC.index(newConf2)] += magneticPrefactor(conf,p,-1,dsh,dsv) * gB
        M[c][c] += electricEnergy(conf) * gE
    return np.array(M)


def localHamiltonian(L: int, j: int, dsh: int = dshGlobal, dsv: int = dsvGlobal):
    CC = chainConfigurations(L,dsh,dsv)
    confNum = len(CC)
    M = [[0 for x in range(confNum)] for y in range(confNum)]
    for c in range(confNum):
        conf = CC[c]
        newConf1 = applyMagnetic(conf,j,1)
        newConf2 = applyMagnetic(conf,j,-1)
        if (coherentInSpin(newConf1,j,dsh,dsv)):
            M[c][CC.index(newConf1)] += magneticPrefactor(conf,j,1,dsh,dsv) * gB
        if (coherentInSpin(newConf2,j,dsh,dsv)):
            M[c][CC.index(newConf2)] += magneticPrefactor(conf,j,-1,dsh,dsv) * gB
        M[c][c] += localElectricEnergy(j, conf) * gE
    return np.array(M)


def magneticHamiltonian(N: int, dsh: int = dshGlobal, dsv: int = dsvGlobal):
    CC = chainConfigurations(N,dsh,dsv)
    L = len(CC)
    M = [[0 for x in range(L)] for y in range(L)]
    for c in range(L):
        conf = CC[c]
        for p in range(N):
            newConf1 = applyMagnetic(conf,p,1)
            newConf2 = applyMagnetic(conf,p,-1)
            if (coherentInSpin(newConf1,p,dsh,dsv)):
                M[c][CC.index(newConf1)] += magneticPrefactor(conf,p,1,dsh,dsv)
            if (coherentInSpin(newConf2,p,dsh,dsv)):
                M[c][CC.index(newConf2)] += magneticPrefactor(conf,p,-1,dsh,dsv)
    return np.array(M)

def electricHamiltonian(N: int, dsh: int = dshGlobal, dsv: int = dsvGlobal):
    CC = chainConfigurations(N,dsh,dsv)
    ListOfConfigs = len(CC)
    M = [[0 for x in range(ListOfConfigs)] for y in range(ListOfConfigs)]
    for c in range(ListOfConfigs):
        conf = CC[c]
        M[c][c] += electricEnergy(conf)
    return np.array(M)

# KET STATES
# An ket operator is a unitary-norm complex (column) vector, which is a numpy array.

def stateNormSquared(state) -> float:
    return np.real(np.vdot(state, state))

def stateNorm(state) -> float:
    return math.sqrt(np.real(np.vdot(state, state)))

def braket(state1, state2) -> complex:
    return np.vdot(state1, state2)

def matrixElement(state1, operator, state2):
    return np.vdot(state1,np.matmul(operator,state2))

# Expectation value of an operator
def expectValue(operator, state):
    return matrixElement(state, operator, state)

# Expectation value of an operator
def expectObsValue(operator, state):
    return np.real(expectValue(operator, state))