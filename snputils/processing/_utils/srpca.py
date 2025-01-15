import numpy as np 
from scipy.sparse.linalg import svds


def prox_fro(X, c):
    n = np.linalg.norm(X, 'fro')
    if n <= c:
        return np.zeros_like(X)
    else:
        return (1 - c / n) * X


def prox_l1(Y, c):
    return np.sign(Y) * np.maximum(np.abs(Y) - c, 0)


def proj_rank_r(Y, k):

    U, S, Vt = svds(Y, k=k)
    
    S[k:len(S)] = 0

    #print(S)

    return U @ np.diag(S) @ Vt


