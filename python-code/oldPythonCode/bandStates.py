from operatorsAndStates import *
import time

def coordinate(z : complex):
    angle = np.angle(z)
    radius = math.sqrt(z.real**2 + z.imag**2)
    n = round(L * angle / np.pi)
    if (n % 2 == 0):
        return [int(n/2), radius]
    else:
        return [int((n - np.sign(n) * L)/2), -radius]

start = time.time()
T = translationOperator(L)
end = time.time()
print("Done constructing T in " + str(end - start) + " seconds")
print("")

# Computing the hamiltonian
start = time.time()
H = totalHamiltonian(L)
end = time.time()
print("Done costructing H in " + str(end - start) + " seconds")
print("")

# Computing the product of H and T
start = time.time()
TH = np.matmul(T,H)
end = time.time()
print("Done the product TH in " + str(end - start) + " seconds")
print("")

# Computing eigenvalues and eigenvectors of HT
start = time.time()
eigensystemTH = np.linalg.eig(TH)
end = time.time()
print("Done Eigensystem of TH in " + str(end - start) + " seconds")
print("")

# Transposing the matrix to obtain column-eigenvectors
blochStates = np.ndarray.transpose(eigensystemTH[1])

eigenvaluesTH = eigensystemTH[0]

# Computing the list of bloch numbers
blochNumber = []
for i in range(len(eigensystemTH[0])):
    blochNumber.append(coordinate(eigensystemTH[0][i])[0])

# Computing the list of eigenvalues of H
eigenvaluesH = []
for i in range(len(eigensystemTH[0])):
    eigenvaluesH.append(coordinate(eigensystemTH[0][i])[1])

def groundStateGenerator():
    for i in range(len(eigenvaluesH)):
        energy = eigenvaluesH[i]
        if (energy < -2.0 / g**2 * (L - 1)):
            return blochStates[i]

def groundStateEnergyGenerator():
    for i in range(len(eigenvaluesH)):
        energy = eigenvaluesH[i]
        if (energy < -2.0 / g**2 * (L - 1)):
            return energy

groundState = groundStateGenerator()

groundStateEnergy = groundStateEnergyGenerator()

# Selecting the first band
def blochEigAndStatesGenerator():
    blochList = []
    for i in range(len(eigenvaluesH)):
        energy = eigenvaluesH[i]
        if (1 < g**2 * energy / 2 + L < 3):
            blochList.append([blochNumber[i], eigenvaluesH[i], blochStates[i]])
    return blochList

blochEigAndStates = blochEigAndStatesGenerator()

def bandstateEnergy(k : int):
    for i in range(len(blochEigAndStates)):
        if (blochEigAndStates[i][0] == k):
            return blochEigAndStates[i][1]

def bandState(k : int):
    for i in range(len(blochEigAndStates)):
        if (blochEigAndStates[i][0] == k):
            return blochEigAndStates[i][2]