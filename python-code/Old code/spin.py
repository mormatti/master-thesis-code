from cmath import sqrt

def spinCoherent(ds: int, dsz: int):
    if (ds < 0):
        #print('Spin negative error')
        return False
    else:
        if (abs(dsz) > ds):
            #print('Sz > S spin error')
            return False
        else:
            if (ds%2 == 0):
                return True if dsz%2 == 0 else False
            else:
                return True if dsz%2 == 1 else False

def ladderOperatorEigenvalue(sign: str, ds: int, dsz: int) -> float:
    c = (-1 if sign == "-" else 1)
    s = ds/2.0
    sz = dsz/2.0
    return sqrt(s * (s + 1) - sz * (sz + c))