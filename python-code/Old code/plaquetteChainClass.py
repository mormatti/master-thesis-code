import plaquette as plq
import spin
import numpy as np
import scipy.sparse.linalg as sla

class PlaquetteChain:
    def __init__(self, plaquettes: list[plq.Plaquette]):
        self._plaquettes = plaquettes

    @property
    def length(self):
        return len(self._plaquettes)

    def compatibleOrder(self) -> bool:
        for i in range(len(self._plaquettes)):
            if not self._plaquettes[i-1].precedes[self._plaquettes[i]]:
                return False
        return True
    
    def allowedIDChain(self) -> list[int]:
        return [p.allowedIDNumber for p in self._plaquettes]
    
    def plaquetteChainFromAllowedIdChain(self, allowedIDChain: list[int]) -> 'PlaquetteChain':
        return [plq.allowedList[i] for i in allowedIDChain]
