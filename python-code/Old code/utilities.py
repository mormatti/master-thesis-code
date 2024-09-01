from ast import List
from plaquette import *

def addOne(L: list[int], base: int):
    L.reverse()
    L[0] += 1
    for i in range(len(L)-1):
        if (L[i] == base):
            L[i] = 0
            L[i+1] += 1
    L.reverse()
    return L

def lastConcatenation(V : List(int), W : List(int)):
    L = []
    for w in W:
        L.append(V + [w])
    return L

def generateNext(L : List(List(int)), book : List(List(int))):
    GN = []
    for l in L:
        GN += lastConcatenation(l,book[l[-1]])
    return GN

def generateNextPeriodic(L : List(List(int)), book : List(List(int))):
    GNP = List(List(int))
    GNP = []
    for l in L:
        Q = List(int)
        Q = []
        b = book[l[-1]]
        for k in b:
            if l[0] in book[k]:
                Q.append(l + [k])
        GNP += Q
    return GNP


def horizArrowRep(integer: int):
    return ("--" if integer == 0 else (str(integer) + ">" if integer >= 0 else "<" + str(abs(integer))))

def vertArrowRepTop(integer: int):
    return ("|" if integer == 0 else ("Λ" if integer >= 0 else str(abs(integer))))

def vertArrowRepBottom(integer: int):
    return ("|" if integer == 0 else (str(integer) if integer >= 0 else "V"))