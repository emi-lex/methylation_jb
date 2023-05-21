import numpy as np

def horvath_res_before_rescale(features, coefs, intercept=0.695507258):
    """
    @param features: horvath features
    @param coefs: horvath coefs
    @param intercept: horvath intercept
    @return: horvath prediction
    """
    return features.dot(coefs) + intercept


def apply_horvath(features, coefs, intercept=0.695507258):
    """
    @param features: horvath features
    @param coefs: horvath coefs
    @param intercept: horvath intercept
    @return: horvath predictions
    """
    return np.vectorize(F_inv)(horvath_res_before_rescale(features, coefs, intercept=intercept))


def F_inv(x):
    """
    @param x: nonlinearized prediction
    @return: linearized prediction
    """
    if x < 0:
        return 21 * np.exp(x) - 1
    else:
        return 21 * x + 20