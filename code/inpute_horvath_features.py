import pandas as pd
import numpy as np
import utils
import sys


def __main__(filename):
    """
    This script imputes the Horvath features using the RRBS data from many people
    @param filename: name of the file (including path) with row RRBS data
    @return: imputed Horvath features
    """
    data = utils.preprocess_data(filename)
    map_data = pd.read_csv("data/hm27.hg19.manifest_preprocessed.bed", sep='\t')
    horvath_data = pd.read_csv("data/gb-2013-14-10-r115-S3_preprocessed.csv", sep='\t')

    features = utils.get_features_with_methiLImp(data, horvath_data["marker"], map_data)

    features_with_names = np.insert(features.astype(str), 0, data.columns[2:], axis=1)
    features_with_names = np.insert(features_with_names, 0, [""] + list(horvath_data["marker"]), axis=0)

    np.savetxt("data/horvath_features.tsv", features_with_names, delimiter='\t', fmt='%s')

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python3 <script_path> <file_path>")
        sys.exit(1)
    
    file_name = sys.argv[1]
    
    __main__(file_name)
