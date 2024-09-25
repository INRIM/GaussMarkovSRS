from dataclasses import dataclass

import numpy as np
import scipy.stats as st

from gauss_markov_srs.utils import cov2unc


@dataclass
class FitResult:
    q: np.ndarray
    Cqq: np.ndarray
    weights: np.ndarray


@dataclass
class FitStats:
    n_meas: int
    n_adj: int
    n_excluded: int
    n_included: int
    n_dof: int
    self_sensitivity: np.ndarray
    residuals: np.ndarray
    relative_residuals: np.ndarray
    chi2: float
    birge: float
    birge_limit: float
    p_value: float
    n_corr_nonzero: int
    min_eigenvalue: float


def gauss_markov_fit(A, y, Cyy, mask=None):
    """Generic Gauss-Markov fit.

    Parameters
    ----------
    A : 2d array, shape (N, M)
        Model matrix.
    y : 1d array, shape (N,)
        Input data.
    Cyy : 2d array, shape (N, N)
        Covariance matrix of input data.
    mask : 1d array, shape (N,), optional
        Mask of valid data point (True=valid), by default None use all data.
    Returns
    -------
    FitResult
        `FitResult` with the following fields defined:

        q : ndarray, shape (M,)
            Adjusted values.
        cov : ndarray, shape (M, M)
            Covariance matrix of the adjusted values.
        weights : ndarray, shape (M, N)
            Matrix of weights to apply to the data to get the fit result.
    """
    A = np.asanyarray(A)
    y = np.asanyarray(y)
    Cyy = np.asanyarray(Cyy)

    if mask is None:
        mask = np.ones_like(y).astype(bool)
    else:
        mask = np.asanyarray(mask).astype(bool)

    nMeas, nAdj = A.shape

    A = A[mask]
    y = y[mask]
    Cyy = Cyy[mask][:, mask]

    invCyy = np.linalg.inv(Cyy)

    Cqq = np.linalg.inv(A.transpose() @ invCyy @ A)
    weights = Cqq @ A.transpose() @ invCyy

    q = np.dot(weights, y).T

    weights_all = np.zeros((nAdj, nMeas))
    weights_all[:, mask] = weights

    ret_q = np.squeeze(np.array(q, ndmin=1))
    ret_C = np.array(Cqq)
    ret_weights = np.squeeze(np.array(weights_all, ndmin=1))

    ret = FitResult(ret_q, ret_C, ret_weights)
    return ret


def fit_stats(A, y, Cyy, fit_result, mask=None):
    """Calculate some statistics from the fit result.

    Parameters
    ----------
    A : 2d array, shape (N, M)
        Model matrix.
    y : 1d array, shape (N,)
        Input data.
    Cyy : 2d array, shape (N, N)
        Covariance matrix of input data.
    q : ndarray, shape (M,)
        Adjusted values.
    cov : ndarray, shape (M, M)
        Covariance matrix of the adjusted values.
    weights : ndarray, shape (M, N)
        Matrix of weights to apply to the data to get the fit result.
    mask : 1d array, shape (N,), optional
        Mask of valid data point (True=valid), by default None use all data.
    Returns
    -------
    FitStats
        `FitStats` with the following fields defined:
            n_meas: int
                Total number of measurement.
            n_adj: int
                Number of adjusted constants.
            n_excluded: int
                Number of measurements excluded from the fit.
            n_included: int
                Number of measurements included in the fit.
            n_dof: int
                Number of degrees of freedom.
            self_sensitivity: np.ndarray, shape (N, )
                Self-sensitivity coefficients for the input data.
            residuals : ndarray, shape (N, )
                Residuals of the fit.
            relative_residuals : ndarray, shape (N, )
                Residuals of the fit scaled to the input uncertainty.
            chi2 : float
                Chi squared (nt reduced) of the fit
            birge : float
                Birge ratio (square root of the reduced chi squared).
            birge_limit: float
                Simple criteria for maximum acceptable Birge ratio
            p_value: float
                p_value of the chi2 (probability of chi2 < chi2observed)
            n_corr_nonzero: int
                Number of non-zero correlation coefficients.
            min_eigenvalue: float
                Square root of the minimum eigenvalue of the input covariance matrix.
    """

    A = np.asanyarray(A)
    y = np.asanyarray(y)
    Cyy = np.asanyarray(Cyy)
    yu = cov2unc(Cyy)

    if mask is None:
        mask = np.ones_like(y).astype(bool)
    else:
        mask = np.asanyarray(mask).astype(bool)

    n_meas, n_adj = A.shape

    # self-sensitivity and residuals can be calculated unmasked
    Sc = np.diag(A @ fit_result.weights)
    Q = y - A @ fit_result.q

    # after that we can mask and carry on
    A = A[mask]
    Cyy = Cyy[mask][:, mask]
    Q_mask = Q[mask]

    invCyy = np.linalg.inv(Cyy)

    n_included, n_adj = A.shape
    n_excluded = n_meas - n_included
    n_dof = n_included - n_adj

    chi2 = Q_mask.T @ invCyy @ Q_mask
    birge = (chi2 / (n_dof)) ** 0.5
    birge_limit = 1 + (2 / n_dof) ** 0.5
    p_value = st.chi2(df=n_dof).cdf(chi2)

    val, _ = np.linalg.eigh(Cyy)
    min_eigenvalue = np.amin(val) ** 0.5
    n_corr_nonzero = np.count_nonzero(Cyy - np.eye(n_included) * yu[mask] ** 2) // 2

    # prepare return arrays
    self_sensitivity = np.squeeze(np.array(Sc, ndmin=1))
    residuals = np.squeeze(np.array(Q, ndmin=1))
    relative_residuals = np.squeeze(np.array(Q / yu, ndmin=1))

    ret = FitStats(
        n_meas,
        n_adj,
        n_excluded,
        n_included,
        n_dof,
        self_sensitivity,
        residuals,
        relative_residuals,
        chi2,
        birge,
        birge_limit,
        p_value,
        n_corr_nonzero,
        min_eigenvalue,
    )

    return ret


if __name__ == "__main__":
    A = np.array([[1, 0], [1, 0.5], [1, 1]])
    y = np.array([0.0, 0.49, 1.01])
    Cyy = np.eye(3)

    ret = gauss_markov_fit(A, y, Cyy)
    stats = fit_stats(A, y, Cyy, ret)

    A = np.array([[1, 0], [1, 0.5], [1, 1], [1, 2]])
    y = np.array([0.0, 0.49, 1.01, 33.3])
    Cyy = np.eye(4)

    ret = gauss_markov_fit(A, y, Cyy, mask=[1, 1, 1, 0])
    stats = fit_stats(A, y, Cyy, ret, mask=[1, 1, 1, 0])
