from scipy.sparse import identity
from scipy.sparse.linalg import eigs
import plaquette as plq
import numpy

#A = identity(10, format='csc')
#A.setdiag(range(1, 11))
# print(eigs(A, 4, sigma=0)) # find three eigenvalues near zero using shift-invert mode
#P = A.power(2)
#A[2,5] = 7
#print(A[3,5])

#import scipy.sparse as sp
#A = sp.identity(10, format='csc') - sp.identity(10, format='csc')

A = plq.Plaquette(-1,-2,1,2)

A.printAsSuperPlaquette()

B = A.clone()

B.printAsSuperPlaquette()

B.magneticOperatorAction("left", False)

B.printAsSuperPlaquette()