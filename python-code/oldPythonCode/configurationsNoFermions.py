import math
import numpy as np
import copy
from plaquettesNoFermions import *

# A configuration is a list of plaquettes. It represents a given chain configuration in
# the spin representation, i.e. a vector of the spin representation basis of the Hilbert
# space of the system.

# The list of the configurations composed of a single plaquette
def startingConfigurations(dsh: int = dshGlobal) -> list[list[list[int]]]:
    P = plaquetteList(dsh)
    C = []
    for p in P:
        C.append([p])
    return C

# Adds a chain element to the list of configurations. If 'periodic' it 'closes' the
# chain with periodic boundary conditions. If periodic = true, use it only for the last
# plaquette!
def addChainElement(configuration: list[list[list[int]]], periodic: bool, dsh: int = dshGlobal, 
                                        dsv: int = dsvGlobal):
    C = []
    P = plaquetteList(dsh)
    for chain in configuration:
        for p in P:
            if follows(chain[-1], p, dsv):
                if (periodic):
                    if follows(p, chain[0], dsv):
                        newChain = chain + [p]
                        C.append(newChain)
                else:
                    newChain = chain + [p]
                    C.append(newChain)
    return C


def chainConfigurations(chainLength: int, dsh: int = dshGlobal, dsv: int = dsvGlobal):
    if chainLength == 0:
        return []
    if chainLength == 1:
        return startingConfigurations()
    else:
        C = startingConfigurations(dsh)
        for i in range(chainLength-2):
            C = addChainElement(C, False, dsh, dsv)
        C = addChainElement(C, True, dsh, dsv)
        return C

# Returns the translated configuration of 'index' plaquette forward
def translated(configuration: list[list[int]], index = 1) -> list[list[int]]:
    return np.roll(configuration, index, axis = 0).tolist()

def invertedField(configuration: list[list[int]]) -> list[list[int]]:
    configurationCopy = copy.deepcopy(configuration)
    for c in configurationCopy:
        for i in range(len(c)):
            c[i] = -1 * c[i]
    return configurationCopy

# Returns the swapped configuration (x reflection)
def xReflected(configuration: list[list[int]]):
    return np.flip(invertedField(configuration), axis = 0).tolist()

def yReflected(configuration: list[list[int]]):
    return np.flip(configuration, axis = 1).tolist()

def plaquetteDoubleSpin(link: str, configuration: list[list[int]], plaquetteIndex: int):
    a = (plaquetteIndex - 1) % len(configuration)
    b = plaquetteIndex % len(configuration)
    c = (plaquetteIndex + 1) % len(configuration)
    if (link == "B"):
        return configuration[b][1]
    elif (link == "T"):
        return configuration[b][0]
    elif (link == "L"):
        return configuration[b][0] - configuration[a][0]
    elif (link == "R"):
        return configuration[c][0] - configuration[b][0]
    else:
        print("Error: link", link, "does not exist: R returned.")
        return configuration[c][0] - configuration[b][0]

def coherentInSpin(configuration: list[list[int]], plaquetteIndex: int, 
                                        dsh: int = dshGlobal, dsv: int = dsvGlobal):
    dsB = plaquetteDoubleSpin("B", configuration, plaquetteIndex)
    dsT = plaquetteDoubleSpin("T", configuration, plaquetteIndex)
    dsL = plaquetteDoubleSpin("L", configuration, plaquetteIndex)
    dsR = plaquetteDoubleSpin("R", configuration, plaquetteIndex)
    return abs(dsB) <= dsh and abs(dsT) <= dsh and abs(dsL) <= dsv and abs(dsR) <= dsv

def applyMagnetic(configuration: list[list[int]], plaquetteIndex: int, 
                                                                ladderIndex: int = 1):
    configurationCopy = copy.deepcopy(configuration)
    plaquetteIndex = plaquetteIndex % len(configurationCopy)
    configurationCopy[plaquetteIndex][1] += 2 * ladderIndex
    configurationCopy[plaquetteIndex][0] -= 2 * ladderIndex
    return configurationCopy

def magneticPrefactor(configuration: list[list[int]], plaquetteIndex: int, 
                    ladderIndex: int = 1, dsh: int = dshGlobal, dsv: int = dsvGlobal):
    sh = dsh / 2.0
    sv = dsv / 2.0
    sB = plaquetteDoubleSpin("B",configuration,plaquetteIndex) / 2.0
    sT = plaquetteDoubleSpin("T",configuration,plaquetteIndex) / 2.0
    sL = plaquetteDoubleSpin("L",configuration,plaquetteIndex) / 2.0
    sR = plaquetteDoubleSpin("R",configuration,plaquetteIndex) / 2.0
    prefactorB = math.sqrt(sh*(sh + 1) - sB*(sB + ladderIndex))
    prefactorT = math.sqrt(sh*(sh + 1) - sT*(sT - ladderIndex))
    prefactorL = math.sqrt(sv*(sv + 1) - sL*(sL - ladderIndex))
    prefactorR = math.sqrt(sv*(sv + 1) - sR*(sR + ladderIndex))
    return prefactorB * prefactorT * prefactorL * prefactorR