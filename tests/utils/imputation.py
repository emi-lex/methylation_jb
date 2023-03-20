import numpy as np
import pandas as pd
# from missforest.miss_forest import MissForest
from .missForest import MissForest
from tqdm import tqdm
from .binary_search import binary_search
from .preprocess import get_chromosome_coords


def impute(X: np.ndarray, head=None, method="zeroes"):
    header = None
    if not head is None:
        header = '\t'.join(head)
        
    if method == "zeroes":
        res = np.concatenate((X[:, 0][:, np.newaxis], np.nan_to_num(X[:, 1:].astype(np.float64), nan=0)), axis=1)
        np.savetxt('../data/zeroes_imputed_cytosines.tsv', res, delimiter='\t', fmt='%s', header=header)
        return res
    
    if method == "1000bp":
        chromosomes_arr, positions_arr = np.transpose([name.split(".") for name in X[:, 0]])
        positions_arr = positions_arr.astype(int)

        # def compare_chr(chr1, p1, chr2, p2):
        #     return chr1 < chr2 or (chr1 == chr2 and p1 < p2)

        # def binary_search(c, p):
        #     left = 0
        #     right = len(chromosome)
        #     while right - left > 1:
        #         mid = (right + left) // 2
        #         if compare_chr(c, p, chromosome[mid], position[mid]):
        #             right = mid
        #         else:
        #             left = mid
        #     return left


        for i, j in tqdm(np.argwhere(np.isnan(X[:, 1:].astype(float)))):
            j += 1
            left_row = binary_search(chromosomes_arr, positions_arr, chromosomes_arr[i], positions_arr[i] - 500, "ge")
            right_row = binary_search(chromosomes_arr, positions_arr, chromosomes_arr[i], positions_arr[i] + 500, "le")
            # if chromosomes_arr[left_row] != chromosomes_arr[i] or positions_arr[left_row] != positions_arr[i] - 500:
            #     left_row += 1
            # if right_row < left_row:
            #     print(right_row, left_row, i, compare_chr(chromosome[i], position[i], chromosome[right_row], position[right_row]))
            # BUG with np.mean
            # if (i == 2):
            #     print(left_row, right_row)
            X[i, j] = np.nanmean(X[left_row:right_row + 1, j])
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


# def impute_rows(data, chr_arr, pos_arr):
#     rows = []
#     for chr, pos in zip(chr_arr, pos_arr):

#         i = binary_search(chromosomes_arr, positions_arr, chr, pos - 500, "ge")
#         j = binary_search(chromosomes_arr, positions_arr, chr, pos + 500, "le")
#         if i > j:
#             neighbours = []
#             if chromosomes_arr[i] == chr:
#                 neighbours.append(data_values[i][1:])
#             if chromosomes_arr[j] == chr:
#                 neighbours.append(data_values[j][1:])
#             if len(neighbours) == 0:
#                 return "AAA error", True
#             return np.nanmean(neighbours, axis=0), True
#         rows.append(np.nanmean(data_values[i:j + 1][:, 1:], axis=0), False)


def binsearch(a, v, side="right"):
    if v < a[0]:
        return 0
    if v >= a[-1]:
        return len(a)
    return np.searchsorted(a, v, side)


def find_slice(pos, eps, positions_arr):
    left = binsearch(positions_arr, pos - eps)
    right = binsearch(positions_arr, pos + eps)
    if left == right:
        if left == 0:
            return right, right + 1
        if right == len(positions_arr):
            return left - 1, left
        return left - 1, right + 1
        
    return left, right


def impute_row(pos, eps, imputed_data, i, j):
    left, right = find_slice(pos, eps, imputed_data.position.values[i:j])
    left += i
    right += i
    return np.nanmean(imputed_data[imputed_data.columns[2:]][left:right], axis=0)


def impute_zeroes(data: pd.DataFrame):
    return data.fillna(0)


def impute_1000bp(data: pd.DataFrame, eps=500): # tested
    imputed_data = data.copy()
    chromosomes_arr = data.chromosome.values
    positions_arr = data.position.values
    values = data.values[:, 2:].astype(np.float64)
    columns = data.columns
    coords = get_chromosome_coords(chromosomes_arr)
    for row, col in tqdm(list(zip(*np.where(np.isnan(values))))):
        i, j = coords[chromosomes_arr[row]]
        left, right = find_slice(positions_arr[row], eps, positions_arr[i:j])
        left += i
        right += i
        imputed_data.loc[row, columns[col + 2]] = np.nanmean(values[left:right, col])
    return impute_zeroes(imputed_data)


# def impute_1000bp(data: pd.DataFrame, epsilon=500, print_progress=True):
#     j = 0
#     n, m = data.shape
#     imputed_data = data.copy()
#     while j < n:
#         chr = data.chromosome[j]
#         if print_progress:
#             print(chr)
#         i = j
#         while j < n and data.chromosome[j] == chr:
#             j += 1
#         # [i, j) is the semi-interval where are all samples from chromosome chr

#         left = i
#         right = i
#         for k in tqdm(range(i, j)):
#             while data.position[left] < data.position[k] - epsilon:
#                 left += 1
#             while data.position[right] <= data.position[k] + epsilon and k < j:
#                 right += 1
#             # all samples that are in chromosome chr and with (pos[k] - eps <= pos <= pos[k] + eps) are in [left, right) 

#             for c in np.where(np.isnan(data.iloc[k][2:].astype(np.float64))):
#                 col = data.columns[c + 2]
#                 if left == right:
#                     if left == i:
#                         imputed_data.loc[k, col] = data[col][right]
#                     elif right == j:
#                         imputed_data.loc[k, col] = data[col][left - 1]
#                     else:
#                         imputed_data.loc[k, col] = np.nanmean(data[col][left - 1:right + 1])
#                 else:
#                     if left >= right:
#                         print(left, right, k)
#                     imputed_data.loc[k, col] = np.nanmean(data[col][left:right])
#     return impute_zeroes(imputed_data)