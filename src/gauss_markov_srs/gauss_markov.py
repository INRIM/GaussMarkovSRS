import numpy as np


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
    tuple
        Tuple with elements:
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

    ret = (ret_q, ret_C, ret_weights)
    return ret


def fit_stats(A, y, Cyy, q, Cqq, weights, mask=None):
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
    tuple
        Tuple with elements:
        self-sensitivity : ndarray, shape (N, )
            Self-sensitivity coefficients for the input data.
        residuals : ndarray, shape (N, )
            Residuals of the fit.
        ndof : int
            Number of degrees of freedom.
        chi2 : float
            Chi squared of the fit
    """

    A = np.asanyarray(A)
    y = np.asanyarray(y)
    Cyy = np.asanyarray(Cyy)

    if mask is None:
        mask = np.ones_like(y).astype(bool)
    else:
        mask = np.asanyarray(mask).astype(bool)

    nMeas, nAdj = A.shape

    # self-sensitivity and residuals can be calculated unmasked
    Sc = np.diag(A @ weights)
    Q = y - A @ q

    # after that we can mask and carry on
    A = A[mask]
    Cyy = Cyy[mask][:, mask]
    Q_mask = Q[mask]

    invCyy = np.linalg.inv(Cyy)

    nMeas_masked, nAdj = A.shape

    ndof = nMeas_masked - nAdj
    chi2 = Q_mask.T @ invCyy @ Q_mask

    # prepare return tuple
    self_sensitivity = np.squeeze(np.array(Sc, ndmin=1))
    residuals = np.squeeze(np.array(Q, ndmin=1))
    ret = (self_sensitivity, residuals, ndof, chi2)

    return ret


if __name__ == "__main__":
    A = np.array([[1, 0], [1, 0.5], [1, 1]])
    y = np.array([0.0, 0.49, 1.01])
    Cyy = np.eye(3)

    ret = gauss_markov_fit(A, y, Cyy)
    stats = fit_stats(A, y, Cyy, *ret)

    A = np.array([[1, 0], [1, 0.5], [1, 1], [1, 2]])
    y = np.array([0.0, 0.49, 1.01, 33.3])
    Cyy = np.eye(4)

    ret = gauss_markov_fit(A, y, Cyy, mask=[1, 1, 1, 0])
    stats = fit_stats(A, y, Cyy, *ret, mask=[1, 1, 1, 0])
