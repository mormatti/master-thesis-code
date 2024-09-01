import numpy as np
import scipy.optimize as opt
import random

I = complex(1,0)

def functional(inputMatrix : list[list[complex]], thetas : list[float]) -> float:
    N = len(inputMatrix)
    totalSum = complex(0,0)
    for i in range(N):
        for j in range(N):
            phaseTerm1 = complex(np.cos(thetas[i]), np.sin(thetas[i]))
            phaseTerm2 = complex(np.cos(thetas[j]), -np.sin(thetas[j]))
            matrixElement = inputMatrix[i][j]
            totalSum += phaseTerm1 * phaseTerm2 * matrixElement
    return np.real(totalSum)

def optimizedThetas(inputMatrix : list[list[complex]]):
    N = len(inputMatrix)
    theta0 = [random.randint(-4,4) for i in range(N)]

    def func(thetas : list[float]) -> float:
        return functional(inputMatrix, thetas)
    return opt.minimize(func, theta0, method="Nelder-Mead")

M = [[complex(1,0), complex(0,1), complex(0,0), complex(0,0), complex(0,0), complex(1,1)],
     [complex(0,-1),complex(2,0), complex(0,0), complex(0,0), complex(2,8), complex(0,0)],
     [complex(0,0), complex(0,0), complex(1,0), complex(2,2), complex(0,0), complex(0,0)],
     [complex(0,0), complex(0,0),complex(2,-2), complex(5,0), complex(1,1), complex(1,5)],
     [complex(0,0),complex(2,-8), complex(0,0),complex(1,-1), complex(4,0), complex(0,0)],
    [complex(1,-1), complex(0,0), complex(0,0),complex(1,-5), complex(0,0), complex(3,0)]]

OT = optimizedThetas(M)
result = OT.x
functionResult = OT.fun
resultUnivoque = (result - result[0]) % (2 * np.pi)
print(resultUnivoque, functionResult)
# print(OT)

