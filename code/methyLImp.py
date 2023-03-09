import numpy as np
from scipy.special import expit

def methyLImp(M):
    M_prime = M.copy()
    R = set(range(M.shape[0]))
    C = set(range(M.shape[1]))
    NA = set(c for c in range(M.shape[1]) if np.any(np.isnan(M[:,c])))
    L = NA.copy()
    while len(L) > 0:
        col = L.pop()
        R_NA = set(r for r in range(M.shape[0]) if np.isnan(M[r,col]))
        C_NA = set(c for c in range(M.shape[1]) if set(r for r in range(M.shape[0]) if np.isnan(M[r, c]))==R_NA)
        if len(R - R_NA) > 0 and len(C - NA) > 0:
            A = M[list(R - R_NA)][:, list(C - NA)]
            B = M[list(R - R_NA)][:, list(C_NA)]
            X = M[list(R_NA)][:, list(C - NA)]
            A_inv = np.linalg.pinv(A)
            for row, values in zip(list(R_NA), expit(X @ A_inv @ np.log((B + 10e-06) / (1 - B))).reshape(len(R_NA), len(C_NA))):
                for col, val in zip(list(C_NA), values):
                    M_prime[row, col] = val
            # M_prime[list(R_NA)][:, list(C_NA)] = expit(X @ A_inv @ np.log(B / (1 - B)))
        L -= C_NA
    return M_prime