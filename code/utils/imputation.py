import numpy as np
import pandas as pd
from tqdm import tqdm
from .preprocess import get_chromosome_coords


def binsearch(a, v, side="right"):
    """
    @param a: array, where to search
    @param v: value, which is needed to be found
    @param side: what side to choose if the value was not found
    @return: position of the value in array
    """
    if v < a[0]:
        return 0
    if v >= a[-1]:
        return len(a)
    return np.searchsorted(a, v, side)


def find_slice(pos, eps, positions_arr):
    """
    @param pos: position around which slice is needed
    @param eps: half of the neighborhood
    @positions array: array of the positions of the cytosines
    @return: needed slice [left, right)
    """
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
    """
    @param pos: position of the row
    @param eps: half of the neighborhood
    @param imputed_data: imputed data
    @params i, j: corresponding chromosome is contained in the interval [i:j)
    @return: imputed row
    """
    left, right = find_slice(pos, eps, imputed_data.position.values[i:j])
    left += i
    right += i
    return np.nanmean(imputed_data[imputed_data.columns[2:]][left:right], axis=0)


def impute_zeroes(data: pd.DataFrame):
    """
    @param data: data that needs to be imputed
    @return: imputed data
    """
    return data.fillna(0)


def impute_nbp(data: pd.DataFrame, eps=500, impute_positions=None, verbose=1):
    """
    @param data: data that needs to be imputed
    @param eps: half of the neighborhood
    @param impute_positions: positions that need to be imputed
    @param verbose: verbosity of the function
    @return: imputed data
    """
    imputed_data = data.copy()
    chromosomes_arr = data.chromosome.values
    positions_arr = data.position.values
    values = data.values[:, 2:].astype(np.float64)
    columns = data.columns
    coords = get_chromosome_coords(chromosomes_arr)

    if impute_positions is None:
        impute_positions = np.isnan(values)
    
    if verbose:
        lst = tqdm(list(zip(*np.where(impute_positions))))
    else:
        lst = list(zip(*np.where(impute_positions)))

    for row, col in lst:
        i, j = coords[chromosomes_arr[row]]
        left, right = find_slice(positions_arr[row], eps, positions_arr[i:j])
        left += i
        right += i
        imputed_data.loc[row, columns[col + 2]] = np.nanmean(values[left:right, col])
    return impute_zeroes(imputed_data)


def impute_people_mean(data, impute_positions=None, verbose=1):
    """
    @param data: data that needs to be imputed
    @param impute_positions: positions that need to be imputed
    @param verbose: verbosity of the function
    @return: imputed data
    """
    imputed_data = data.copy()
    values = data.values[:, 2:].astype(np.float64)
    columns = data.columns
    means = np.nanmean(values, axis=0)

    if impute_positions is None:
        impute_positions = np.isnan(values)
    
    if verbose:
        lst = tqdm(list(zip(*np.where(impute_positions))))
    else:
        lst = list(zip(*np.where(impute_positions)))

    for row, col in lst:
        imputed_data.loc[row, columns[col + 2]] = means[col]
    return imputed_data


def impute_cytosine_mean(data, impute_positions=None, verbose=1):
    """
    @param data: data that needs to be imputed
    @param impute_positions: positions that need to be imputed
    @param verbose: verbosity of the function
    @return: imputed data
    """
    imputed_data = data.copy()
    values = data.values[:, 2:].astype(np.float64)
    columns = data.columns
    means = np.nanmean(values, axis=1)

    if impute_positions is None:
        impute_positions = np.isnan(values)
    
    if verbose:
        lst = tqdm(list(zip(*np.where(impute_positions))))
    else:
        lst = list(zip(*np.where(impute_positions)))

    for row, col in lst:
        imputed_data.loc[row, columns[col + 2]] = means[row]
    return imputed_data
