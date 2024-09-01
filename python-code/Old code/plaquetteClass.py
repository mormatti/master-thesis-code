import spin
import numpy as np
import scipy.sparse.linalg as sla
import parameters as params
import utilities as utils

class Plaquette:

    # Plaquette initializer
    def __init__(self, dszB: int, dszR: int, dszT: int, dszL: int):
        self._dszL = dszL
        self._dszB = dszB
        self._dszT = dszT
        self._dszR = dszR
        self._allowed = self.checkIfAllowed()
        self._allowedIDNumber = None

    # Left link property
    def get_dszL(self):
        return self._dszL
    def set_dszL(self, dszL):
        self._dszL = dszL
        self.updateProperties()
    property(get_dszL, set_dszL)

    # Bottom link property
    def get_dszB(self):
        return self._dszB
    def set_dszB(self, dszB):
        self._dszB = dszB
        self.updateProperties()
    property(get_dszB, set_dszB)

    # Top link property
    def get_dszT(self):
        return self._dszT
    def set_dszT(self, dszT):
        self._dszT = dszT
        self.updateProperties()
    property(get_dszT, set_dszT)

    # Right link property
    def get_dszR(self):
        return self._dszR
    def set_dszR(self, dszR):
        self._dszR = dszR
        self.updateProperties()
    property(get_dszR, set_dszR)

    def updateAllowed(self):
        self._allowed = self.checkIfAllowed()

    def updateAllowedIDNumber(self):
        self._allowedIDNumber = self.computeAllowedIDNumber()
    
    def updateProperties(self):
        self.updateAllowed()
        self.updateAllowedIDNumber()

    # External top-left link property
    @property  
    def dszExtTL(self) -> int:   
        return self._dszT - self._dszL

    # External top-right link property
    @property  
    def dszExtTR(self) -> int:   
        return self._dszT + self._dszR
    
    # External bottom-left link property
    @property  
    def dszExtBL(self) -> int:   
        return self._dszB + self._dszL
    
    # External bottom-right link property
    @property  
    def dszExtBR(self) -> int:   
        return self._dszB - self._dszR

    @property
    def allowedIDNumber(self) -> int:
        return self._allowedIDNumber

    def checkIfAllowed(self) -> bool:
        for dsz in [self._dszB, self._dszT, self.dszExtTL, self.dszExtTR, self.dszExtBL, self.dszExtBR]:
            if not spin.spinCoherent(params.dsh, dsz):
                return False
        for dsz in [self._dszL, self._dszR]:
            if not spin.spinCoherent(params.dsv, dsz):
                return False
        return True
    
    def equalsTo(self, plaq: 'Plaquette') -> bool:
        return self._dszL == plaq._dszL and self._dszB == plaq._dszB \
             and self._dszT == plaq._dszT and self._dszR == plaq._dszR

    def precedes(self, plaq: 'Plaquette') -> bool:
        return self._dszT == plaq.dszExtTL and self._dszB == plaq.dszExtBL \
             and self._dszR == plaq._dszL and self.dszExtTR == plaq._dszT and self.dszExtBR == plaq._dszB

    def follows(self, plaq: 'Plaquette') -> bool:
        return plaq.precedes(self)

    def clone(self) -> 'Plaquette':
        P = Plaquette(self._dszB, self._dszR, self._dszT, self._dszL)
        P.updateProperties()
        return P

    @staticmethod    
    def computeAllowedList() -> list['Plaquette']:
        L = list[Plaquette]()
        for dszL in range(-params.dsv,params.dsv+1):
            for dszB in range(-params.dsh,params.dsh+1):
                for dszT in range(-params.dsh,params.dsh+1):
                    for dszR in range(-params.dsv,params.dsv+1):
                        P = Plaquette(dszB,dszR,dszT,dszL)
                        if (P._allowed):
                            L.append(P)
        return L

    # The Adjacency matrix A has the property that A[i][j] = 1 if the plaquette i precedes 
    # the plaquette j, 0 otherwise.
    @staticmethod
    def computeAdjacencyMatrix() -> list[list[int]]:
        L = Plaquette.computeAllowedList()
        n = len(L)
        M = [[0 for x in range(0,n)] for y in range(0,n)]
        for i in range(n):
            for j in range(n):
                if L[i].precedes(L[j]):
                    M[i][j] = 1
        return M

    @staticmethod
    def computeAdjacencyBook():
        L = Plaquette.computeAllowedList()
        B = []
        for i in range(len(L)):
            subL = []
            for j in range(len(L)):
                if L[i].precedes(L[j]):
                    subL.append(j)
            B.append(subL)
        return B
    
    def computeAllowedIDNumber(self) -> int:
        for i in range(len(allowedList)):
            if allowedList[i].equalsTo(self):
                return i
        return len(allowedList)
    
    def printAsSuperPlaquette(self):
        print(utils.horizArrowRep(self.dszExtTL), "*", utils.horizArrowRep(self._dszT), "*", utils.horizArrowRep(self.dszExtTR))
        print("  ", utils.vertArrowRepTop(self._dszL), "  ", utils.vertArrowRepTop(self._dszR), " ")
        print("  ", utils.vertArrowRepBottom(self._dszL), "  ", utils.vertArrowRepBottom(self._dszR), " ")
        print(utils.horizArrowRep(self.dszExtBL), "*", utils.horizArrowRep(self._dszB), "*", utils.horizArrowRep(self.dszExtBR))

    # Action of the operators

    # The sum of all electric field squared. 
    # The left and right links are counted half (for the chain).
    def electricOperatorEigenvalue(self) -> float:
        return (self._dszL / 2.0)**2 + (self._dszB / 2.0)**2 \
                    + (self._dszR / 2.0)**2 + (self._dszT / 2.0)**2
    
    # The action of the magnetic operator.
    def magneticOperatorAction(self, position: str = "central", dagger: bool = False):
        c = (-2 if dagger else 2)
        if position == "left":
            self._dszL = self._dszL + c
        elif position == "right":
            self._dszR = self._dszR - c
        else:
            self._dszL = self._dszL - c
            self._dszB = self._dszB + c
            self._dszT = self._dszT - c
            self._dszR = self._dszR + c
        self.updateProperties()

    def magneticOperatorMatrixElement(self, dagger: bool = False) -> float:
        return spin.ladderOperatorEigenvalue("+" if dagger else "-", params.dsv, self._dszL) \
                * spin.ladderOperatorEigenvalue("-" if dagger else "+", params.dsh, self._dszB) \
                * spin.ladderOperatorEigenvalue("+" if dagger else "-", params.dsv, self._dszT) \
                * spin.ladderOperatorEigenvalue("-" if dagger else "+", params.dsh, self._dszR)
    

allowedList : list['Plaquette'] = Plaquette.computeAllowedList()
adjacencyMatrix : list[list[int]] = Plaquette.computeAdjacencyMatrix()
adjacencyBook : list[list[int]] = Plaquette.computeAdjacencyBook()
allowedNumber : int = len(allowedList)