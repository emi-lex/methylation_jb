import numpy as np
import pandas as pd
# from missforest.miss_forest import MissForest
import missForest
from tqdm import tqdm
from .binary_search import binary_search


def impute(X: np.ndarray, head=None, method="zeroes"):
    header = None
    if not head is None:
        header = '\t'.join(head)
        
    if method == "zeroes":
        res = np.concatenate((X[:, 0][:, np.newaxis], np.nan_to_num(X[:, 1:].astype(np.float64), nan=0)), axis=1)
        np.savetxt('../data/zeroes_imputed_cytosines.tsv', res, delimiter='\t', fmt='%s', header=header)
        return res
    
    if method == "1000bp":
        chromosome, position = np.transpose([name.split(".") for name in X[:, 0]])
        position = position.astype(int)

        def compare_chr(chr1, p1, chr2, p2):
            return chr1 < chr2 or (chr1 == chr2 and p1 < p2)

        def binary_search(c, p):
            left = 0
            right = len(chromosome)
            while right - left > 1:
                mid = (right + left) // 2
                if compare_chr(c, p, chromosome[mid], position[mid]):
                    right = mid
                else:
                    left = mid
            return left


        for i, j in tqdm(np.argwhere(np.isnan(X[:, 1:].astype(float)))):
            j += 1
            left_row = binary_search(chromosome[i], position[i] - 500)
            right_row = binary_search(chromosome[i], position[i] + 500)
            if chromosome[left_row] != chromosome[i] or position[left_row] != position[i] - 500:
                left_row += 1
            if right_row < left_row:
                print(right_row, left_row, i, compare_chr(chromosome[i], position[i], chromosome[right_row], position[right_row]))
            X[i, j] = X[left_row:right_row + 1, j].mean()
        res = np.concatenate((X[:, 0][:, np.newaxis], np.nan_to_num(X[:, 1:].astype(np.float64), nan=0)), axis=1)
        np.savetxt('../data/1000bp_imputed_cytosines.tsv', res, delimiter='\t', fmt='%s', header=header)
        return res
    
    if method == "missForest":
        imputer = missForest.MissForest()
        # return pd.DataFrame(X)
        X_imputed = imputer.fit_transform(pd.DataFrame(X))
        np.savetxt('../data/missForest_imputed_cytosines.tsv', X_imputed, delimiter='\t', fmt='%s', header=header)
        return X_imputed
        # X_imputed = imputer.fit_transform(pd.DataFrame(X[:, 1:]))
        # res = np.concatenate((X[:, 0], X_imputed), axis=1)
        # np.savetxt('../data/missForest_imputed_cytosines.tsv', res, delimiter='\t', fmt='%s', header=header)
        # return res


def impute_row(chromosomes_arr, positions_arr, chr, pos, data_values):
    i = binary_search(chromosomes_arr, positions_arr, chr, pos - 500, "ge")
    j = binary_search(chromosomes_arr, positions_arr, chr, pos + 500, "le")
    if i > j:
        neighbours = []
        if chromosomes_arr[i] == chr:
            neighbours.append(data_values[i][1:])
        if chromosomes_arr[j] == chr:
            neighbours.append(data_values[j][1:])
        if len(neighbours) == 0:
            return "AAA error", True
        return np.mean(neighbours, axis=0), True
    return data_values[i:j + 1][:, 1:].mean(axis=0), False