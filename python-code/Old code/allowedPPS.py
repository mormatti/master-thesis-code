from ast import List
import plaquette as plq
import utilities as utils

def allowedPlaquettePeriodicSequence(S: List(int)):
    n = len(S)
    for i in range(n):
        if 0 <= S[i] < plq.plaquetteAllowedNumber:
            j = (i+1) % n
            if 0 <= S[j] < plq.plaquetteAllowedNumber:
                if plq.plaquetteAdjacencyMatrix[S[i]][S[j]] != 0:
                    pass
                else:
                    return False
            else:
                return False
        else:
            return False
    return True

allowedPlaquettePeriodicSequence([2,2,2,2,2])

def allowedSequencesList(n : int):
    ASL = []
    z = [0] * n
    for i in range(plq.plaquetteAllowedNumber**n):
        if (allowedPlaquettePeriodicSequence(z)):
            ASL.append(z)
        z = utils.addOne(z,plq.plaquetteAllowedNumber)
    return ASL