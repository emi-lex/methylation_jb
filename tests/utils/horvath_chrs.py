import pandas as pd
import numpy as np
from .imputation import impute_row
from tqdm import tqdm
from .binary_search import binary_search
from .preprocess import get_chromosome_coords


# def get_chrs_from_markers(map_file, markers): # tested
#     map_table=pd.read_csv(map_file, sep='\t', header=None)
#     marker_to_chr_pos_tuple = dict(zip(map_table[4], zip(map_table[0], map_table[2] - 1)))
#     chrs = markers.apply(marker_to_chr_pos_tuple.get)
#     return chrs


# def get_imputed_features(chromosomes_arr, positions_arr, data_values, horvath_chrs):
#     features_arr = []
#     counter1, counter2 = 0, 0
#     horvath_ind_left = horvath_chrs.apply(lambda tup: binary_search(chromosomes_arr, positions_arr, *tup))
#     horvath_ind_right = horvath_chrs.apply(lambda tup: binary_search(chromosomes_arr, positions_arr, tup[0], tup[1] + 1))
#     for left_ind, right_ind, (chr, pos) in tqdm(list(zip(horvath_ind_left, horvath_ind_right, horvath_chrs))):
#         if chromosomes_arr[left_ind] == chr and positions_arr[left_ind] == pos:
#             res = data_values[int(left_ind), 1:]
#             if chromosomes_arr[right_ind] == chr and positions_arr[right_ind] == pos + 1:
#                 res = np.mean([res, data_values[int(right_ind), 1:]], axis=0)
#             features_arr.append(res)
#         elif chromosomes_arr[right_ind] == chr and positions_arr[right_ind] == pos + 1:
#             res = data_values[int(right_ind), 1:]
#             features_arr.append(res)
#         else:
#             res1, sig1 = impute_row(chromosomes_arr, positions_arr, chr, pos, data_values)
#             res2, sig2 = impute_row(chromosomes_arr, positions_arr, chr, pos + 1, data_values)
#             counter1 += sig1
#             counter2 += sig2
#             features_arr.append(np.mean([res1, res2], axis=0))
#     features_arr = np.array(features_arr)
#     return features_arr


def get_marker_to_chr_and_pos(map_data): # tested
    marker_to_chr = dict(zip(map_data.marker, map_data.chromosome)).get
    marker_to_pos = dict(zip(map_data.marker, map_data.position)).get
    return marker_to_chr, marker_to_pos


def get_features(imputed_data: pd.DataFrame, markers: np.ndarray, map_data: pd.DataFrame, eps: int=500): # tested
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
