import numpy as np
import pandas as pd
from missforest.miss_forest import MissForest
from tqdm import tqdm

def impute(X: np.ndarray, method="zeroes"):
    if method == "zeroes":
        return np.nan_to_num(X, nan=0)
    if method == "1000bp":
        chromosome, position = np.transpose([name.split(".") for name in X[:, 0]])
        position = position.astype(int)
        for i, j in tqdm(np.argwhere(np.isnan(X[:, 1:].astype(float)))):
            j += 1
            appropriate_rows = (chromosome == chromosome[i]) & (position > position[i] - 500) & (position < position[i] + 500)
            X[i, j] = X[appropriate_rows, j].mean()
        return np.nan_to_num(X, nan=0)
    if method == "missForest":
        imputer = MissForest()
        X_imputed = imputer.fit_transform(pd.DataFrame(X[:, 1:]))
        return np.concatenate((X[:, 0], X_imputed), axis=1)