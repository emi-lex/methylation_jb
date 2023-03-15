import pandas as pd
import numpy as np
from .imputation import impute_row
from tqdm import tqdm
from .binary_search import binary_search


def get_chrs_from_markers(map_file, markers):
    map_table=pd.read_csv(map_file, sep='\t', header=None)
    marker_to_chr_pos_tuple = dict(zip(map_table[4], zip(map_table[0], map_table[2] - 1)))
    chrs = markers.apply(marker_to_chr_pos_tuple.get)
    return chrs


def get_imputed_features(chromosomes_arr, positions_arr, data_values, horvath_chrs):
    features_arr = []
    counter1, counter2 = 0, 0
    horvath_ind_left = horvath_chrs.apply(lambda tup: binary_search(chromosomes_arr, positions_arr, *tup))
    horvath_ind_right = horvath_chrs.apply(lambda tup: binary_search(chromosomes_arr, positions_arr, tup[0], tup[1] + 1))
    for left_ind, right_ind, (chr, pos) in tqdm(list(zip(horvath_ind_left, horvath_ind_right, horvath_chrs))):
        if chromosomes_arr[left_ind] == chr and positions_arr[left_ind] == pos:
            res = data_values[int(left_ind), 1:]
            if chromosomes_arr[right_ind] == chr and positions_arr[right_ind] == pos + 1:
                res = np.mean([res, data_values[int(right_ind), 1:]], axis=0)
            features_arr.append(res)
        elif chromosomes_arr[right_ind] == chr and positions_arr[right_ind] == pos + 1:
            res = data_values[int(right_ind), 1:]
            features_arr.append(res)
        else:
            res1, sig1 = impute_row(chromosomes_arr, positions_arr, chr, pos, data_values)
            res2, sig2 = impute_row(chromosomes_arr, positions_arr, chr, pos + 1, data_values)
            counter1 += sig1
            counter2 += sig2
            features_arr.append(np.mean([res1, res2], axis=0))
    features_arr = np.array(features_arr)
    return features_arr