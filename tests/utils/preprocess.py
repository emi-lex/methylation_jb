import pandas as pd
import numpy as np


def preprocess_data(file_name: str, save: bool=True): # tested
    data = pd.read_csv(file_name, sep="\t")
    chromosomes_arr, positions_arr = np.array([s.split(".") for s in data.chrBase]).T
    positions_arr = positions_arr.astype(np.int32)
    data.drop("chrBase", axis=1, inplace=True)
    data = data.astype(np.float64)
    data.insert(0, "chromosome", chromosomes_arr)
    data.insert(1, "position", positions_arr)
    data.sort_values(["chromosome", "position"], inplace=True, ignore_index=True)
    if save:
        data.to_csv(file_name[:-4] + "_preprocessed" + file_name[-4:], sep="\t", index=False)
    return data


def preprocess_map_file(map_name: str, save: bool=True): # tested
    map_data = pd.read_csv(map_name, names=["chromosome", "position", "garbage1", "garbage2", "marker"], sep="\t")
    map_data.drop(["garbage1", "garbage2"], axis=1, inplace=True)
    map_data.position = map_data.position.astype(np.int32) + 1
    if save:
        map_data.to_csv(map_name[:-4] + "_preprocessed" + map_name[-4:], sep="\t", index=False)
    return map_data


def preprocess_horvath_file(horvath_name: str, save=True): # tested
    horvath_data = pd.read_csv(horvath_name, sep=",")
    horvath_data.drop(horvath_data.columns[2:], axis=1, inplace=True)
    horvath_data.rename(columns={"CpGmarker": "marker", "CoefficientTraining": "coefs"}, inplace=True)
    if save:
        horvath_data.to_csv(horvath_name[:-4] + "_preprocessed" + horvath_name[-4:], sep="\t", index=False)
    return horvath_data


def preprocess_all(file_name: str, map_name: str, horvath_name: str): # tested
    preprocess_data(file_name, True)
    preprocess_map_file(map_name, True)
    preprocess_horvath_file(horvath_name, True)


def get_chromosome_coords(chromosomes_arr): # tested
    i = 0
    n = len(chromosomes_arr)
    chromosome_coords = {}
    while i < n:
        start = i
        chr = chromosomes_arr[i]
        while i < n and chromosomes_arr[i] == chr:
            i += 1
        end = i
        chromosome_coords[chr] = (start, end)
    return chromosome_coords
