import numpy as np


def cov2corr(cov, return_std=False):
    """Convert covariance matrix to correlation matrix.

    Parameters
    ----------
    cov : array_like
        Covariance matrix.
    return_std : bool, optional
        If true, also return the standard deviation, by default False.

    Returns
    -------
    ndarray
        Correlation matrix.
    ndarray
        Array of standard deviations (if return_std = True).

    """
    cov = np.asanyarray(cov)
    std = np.sqrt(np.diag(cov))
    corr = cov / np.outer(std, std)
    if return_std:
        return corr, std
    else:
        return corr


def corr2cov(corr, std):
    """
    Convert correlation matrix to covariance matrix.

    Parameters
    ----------
    corr : array_like
        Correlation matrix.
    std : array_like
        Array of standard deviation.

    Returns
    -------
    cov : ndarray
        Covariance matrix.

    """
    corr = np.asanyarray(corr)
    std = np.asanyarray(std)
    cov = corr * np.outer(std, std)
    return cov
