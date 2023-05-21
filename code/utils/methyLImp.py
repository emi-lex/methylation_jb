import numpy as np
from scipy.special import expit
from tqdm import tqdm
from .preprocess import get_chromosome_coords


def methyLImpSingleChromosome(M, cols=None):
    """
    @param M: methylation data, features are in columns
    @param cols: in which columns missing values whould be imputed
    """
    if cols is None:
        cols = np.ones(M.shape[1], dtype=bool)
    M_imputed = M.copy()
    M_NA = np.isnan(M)
    NA = np.any(M_NA, axis=0)
    imputed_number = 0
    for R_NA in np.unique(M_NA[:, cols & NA], axis=1).T:
        C_NA = np.invert(np.any(R_NA[:, None] ^ M_NA, axis=0)) & cols
        R_wo_NA = np.invert(R_NA)
        C_wo_NA = np.invert(NA)
        A = M[R_wo_NA][:, C_wo_NA]
        B = M[R_wo_NA][:, C_NA]
        X = M[R_NA][:, C_wo_NA]
        A_pinv = np.linalg.pinv(A)
        row_idx, col_idx = np.ix_(R_NA, C_NA)
        M_imputed[row_idx, col_idx] = expit(X @ A_pinv @ np.log((B + 1e-10) / (1 - B + 1e-10))).reshape(R_NA.sum(), C_NA.sum())
        imputed_number += C_NA.sum() * R_NA.sum()
    return M_imputed, imputed_number


def methyLImp(data, rows=None, eps=None, verbose=1):
    """
    @param data: columns are samples, rows are features.
    data contains chromosome and position in first 2 columns
    @param rows: rows that have to be imputed
    @param eps: for each missing values
    the nearest eps features up and down are used for the imputation
    @param verbose: controls the amount of the printed output of this function 
    """
    methylations = data.values[:, 2:].astype(float)
    methylations_imputed = methylations.copy()
    chromosomes = data.chromosome.values
    if rows is None:
        rows = np.ones(data.shape[0], dtype=bool)
    number_imputed_total = 0
    if verbose:
        lst = tqdm(get_chromosome_coords(chromosomes).values())
    else:
        lst = get_chromosome_coords(chromosomes).values()
    for i, j in lst:
        M = methylations[i:j].T / 100
        cols = rows[i:j] # rows become cols because methyLymp works with transposed data
        if eps is None:
            M_imputed, number_imputed = methyLImpSingleChromosome(M, cols)
        else:
            M_imputed = methylations[i:j].T
            number_imputed = 0
            for col in np.where(cols)[0]:
                l = max(0, col - eps)
                r = min(j - i, col + eps)
                cols_mask = np.zeros(r - l, dtype=bool)
                cols_mask[col - l] = 1
                imputed_m, n_imputed = methyLImpSingleChromosome(M[:, l:r], cols_mask)
                number_imputed += n_imputed
                M_imputed[:, col] = 100 * imputed_m[:, col - l]
        methylations_imputed[i:j] = M_imputed.T
        number_imputed_total += number_imputed
    if verbose == 2:
        print("number of used methylations:", rows.sum() * 40)
        print(number_imputed_total, "of them imputed by methyLImp")
    return methylations_imputed
