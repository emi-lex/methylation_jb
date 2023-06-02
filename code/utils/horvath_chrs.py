import pandas as pd
import numpy as np
from .imputation import impute_row, find_slice
from .preprocess import get_chromosome_coords
from .methyLImp import methyLImp
from functools import reduce


def get_marker_to_chr_and_pos(map_data):
    """
    @param map_data: markers mapping information
    @return: dictionaries from markers to chromosomes and positions respectively
    """
    marker_to_chr = dict(zip(map_data.marker, map_data.chromosome)).get
    marker_to_pos = dict(zip(map_data.marker, map_data.position)).get
    return marker_to_chr, marker_to_pos


def get_features(imputed_data: pd.DataFrame, markers: np.ndarray, map_data: pd.DataFrame, eps: int=500):
    """
    @param imputed_data: imputed data
    @param markers: horvath markers
    @param map_data: markers mapping information
    @param eps: half of the neighborhood
    """
    features = []
    get_chr, get_pos = get_marker_to_chr_and_pos(map_data)
    chromosome_coords = get_chromosome_coords(imputed_data.chromosome)
    methylations = imputed_data[imputed_data.columns[2:]].values
    positions_arr = imputed_data.position.values
    for marker in markers:
        chr = get_chr(marker)
        pos = get_pos(marker)
        i, j = chromosome_coords[chr]

        ind1 = np.searchsorted(positions_arr[i:j], pos) + i
        if ind1 < len(positions_arr) and positions_arr[ind1] == pos:
            if ind1 + 1 < imputed_data.shape[0] and positions_arr[ind1 + 1] == pos + 1:
                features.append(methylations[ind1:ind1 + 2].mean(axis=0))
            else:
                features.append(methylations[ind1])
        elif ind1 < len(positions_arr) and positions_arr[ind1] == pos + 1:
            features.append(methylations[ind1])
        else:
            # both chr.pos and chr.(pos+1) don't exist in our data, so this feature needs to be imputed by n_bp
            features.append(impute_row(pos, eps, imputed_data, i, j))
    
    return np.array(features).T / 100


def get_features_with_zeroes(imputed_data: pd.DataFrame, markers: np.ndarray, map_data: pd.DataFrame, eps: int=500):
    """
    @param imputed_data: imputed data
    @param markers: horvath markers
    @param map_data: markers mapping information
    @param eps: half of the neighborhood
    """
    features = []
    get_chr, get_pos = get_marker_to_chr_and_pos(map_data)
    chromosome_coords = get_chromosome_coords(imputed_data.chromosome)
    methylations = imputed_data[imputed_data.columns[2:]].values
    positions_arr = imputed_data.position.values
    for marker in markers:
        chr = get_chr(marker)
        pos = get_pos(marker)
        i, j = chromosome_coords[chr]

        ind1 = np.searchsorted(positions_arr[i:j], pos) + i
        if ind1 < len(positions_arr) and positions_arr[ind1] == pos:
            if ind1 + 1 < imputed_data.shape[0] and positions_arr[ind1 + 1] == pos + 1:
                features.append(methylations[ind1:ind1 + 2].mean(axis=0))
            else:
                features.append(methylations[ind1])
        elif ind1 < len(positions_arr) and positions_arr[ind1] == pos + 1:
            features.append(methylations[ind1])
        else:
            features.append(np.zeros(methylations.shape[1]))
    
    return np.array(features).T / 100


def get_features_with_methiLImp(data: pd.DataFrame, markers: np.ndarray, map_data: pd.DataFrame, eps_methyLImp: int=3000, eps_nbp: int=300):
    """
    @param imputed_data: imputed data
    @param markers: horvath markers
    @param map_data: markers mapping information
    @param eps: half of the neighborhood
    """
    features = []
    get_chr, get_pos = get_marker_to_chr_and_pos(map_data)
    chromosome_coords = get_chromosome_coords(data.chromosome)
    positions_arr = data.position.values
    feature_masks = []
    for marker in markers:
        chr = get_chr(marker)
        pos = get_pos(marker)
        i, j = chromosome_coords[chr]

        ind1 = np.searchsorted(positions_arr[i:j], pos) + i

        mask = np.zeros(data.shape[0], dtype=bool)
        if ind1 < len(positions_arr) and positions_arr[ind1] == pos:
            if ind1 + 1 < data.shape[0] and positions_arr[ind1 + 1] == pos + 1:
                mask[ind1:ind1 + 2] = True
            else:
                mask[ind1] = True
        elif ind1 < len(positions_arr) and positions_arr[ind1] == pos + 1:
            mask[ind1] = True
        else:
            # both chr.pos and chr.(pos+1) don't exist in our data, so this feature needs to be imputed by n_bp
            # features.append(impute_row(pos, eps, imputed_data, i, j))
            left, right = find_slice(pos, eps_nbp, positions_arr[i:j])
            left += i
            right += i
            mask[left:right] = True
        feature_masks.append(mask)
    methylations = methyLImp(data, reduce(lambda arr1, arr2: arr1 | arr2, feature_masks))
    features = np.array([np.nanmean(methylations[mask], axis=0) for mask in feature_masks])
    
    return np.array(features).T
