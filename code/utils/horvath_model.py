import numpy as np

def horvath_res_before_rescale(features, coefs, intercept=0.695507258):
    return (features / 100).astype(np.float64).transpose().dot(coefs.astype(np.float64)) + intercept


def apply_horvath(features, coefs, intercept=0.695507258):
    return np.vectorize(F_inv)(horvath_res_before_rescale(features, coefs, intercept=intercept))


def F_inv(x):
    if x < 0:
        return 21 * np.exp(x) - 1
    else:
        return 21 * x + 20