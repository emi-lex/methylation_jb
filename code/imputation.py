import numpy as np
import pandas as pd
from missforest.miss_forest import MissForest
from tqdm import tqdm

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
        for i, j in tqdm(np.argwhere(np.isnan(X[:, 1:].astype(float)))):
            j += 1
            appropriate_rows = (chromosome == chromosome[i]) & (position > position[i] - 500) & (position < position[i] + 500)
            X[i, j] = X[appropriate_rows, j].mean()
        res = np.concatenate((X[:, 0][:, np.newaxis], np.nan_to_num(X[:, 1:].astype(np.float64), nan=0)), axis=1)
        np.savetxt('../data/1000bp_imputed_cytosines.tsv', res, delimiter='\t', fmt='%s', header=header)
        return res
    if method == "missForest":
        imputer = MissForest()
        X_imputed = imputer.fit_transform(pd.DataFrame(X[:, 1:]))
        res = np.concatenate((X[:, 0], X_imputed), axis=1)
        np.savetxt('../data/missForest_imputed_cytosines.tsv', res, delimiter='\t', fmt='%s', header=header)
        return res