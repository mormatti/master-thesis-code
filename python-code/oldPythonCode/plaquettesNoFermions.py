from globalVariables import *

def follows(p: list[int], q: list[int], dsv: int = dsvGlobal):
    if (len(p) == 2 and len(q) == 2):
        if (p[0] + p[1] == q[0] + q[1]):
            if (abs(q[0] - p[0]) <= dsv):
                return True
    return False

def plaquetteList(dsh: int = dshGlobal):
    L = []
    for ds1 in range(dsh + 1):
        for ds2 in range(dsh + 1):
            L.append([2*dsh - 2*ds1 - dsh, 2*dsh - 2*ds2 - dsh])
    return L

def adjacencyMatrix(dsh: int = dshGlobal, dsv: int = dsvGlobal):
    P = plaquetteList(dsh)
    M = [[0 for x in range((dsh+1)**2)] for y in range((dsh+1)**2)]
    for i in range((dsh+1)**2):
        for j in range((dsh+1)**2):
            if follows(P[i],P[j],dsv):
                M[i][j] = 1
    return M