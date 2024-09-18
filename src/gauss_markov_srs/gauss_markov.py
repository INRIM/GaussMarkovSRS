from dataclasses import dataclass

import numpy as np


@dataclass
class FitResult:
    q: np.ndarray
    cov: np.ndarray
    weights: np.ndarray
    self_sensitivity: np.ndarray
    residuals: np.ndarray
    ndof: int
    chi2: float


def gauss_markov_fit(A, y, Cyy):
    """Generic Gauss-Markov fit.

    Parameters
    ----------
    A : 2d array, shape (N, M)
        Model matrix.
    y : 1d array, shape (N,)
        Input data.
    Cyy : 2d array, shape (N, N)
        Covariance matrix of input data.

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
        self_sensitivity : ndarray, shape (N,)
            Self-sensitivity coefficients for the input data.
        residuals : ndarray, shape (N, )
            Residuals of the fit.
        ndof : int
            Number of degrees of freedom.
        chi2 : float
            Chi squared of the fit

    """
    A = np.matrix(A)
    Cyy = np.matrix(Cyy)

    Cqq = (A.transpose() * Cyy**-1 * A) ** -1
    weights = Cqq * A.transpose() * Cyy**-1

    q = np.dot(weights, y).T

    nm, na = A.shape
    ndof = nm - na

    Q = np.matrix(y).T - A * q
    chi2 = (Q.T * Cyy**-1 * Q)[0, 0]

    Sc = np.diag(A * weights)

    ret_q = np.squeeze(np.array(q, ndmin=1))
    ret_C = np.array(Cqq)
    ret_weights = np.squeeze(np.array(weights, ndmin=1))
    self_sensitivity = np.squeeze(np.array(Sc, ndmin=1))
    residuals = np.squeeze(np.array(Q, ndmin=1))

    ret = FitResult(
        q=ret_q,
        cov=ret_C,
        weights=ret_weights,
        self_sensitivity=self_sensitivity,
        residuals=residuals,
        ndof=ndof,
        chi2=chi2,
    )

    return ret


if __name__ == "__main__":
    A = np.array([[1, 0], [1, 0.5], [1, 1]])
    y = np.array([0.0, 0.49, 1.01])
    Cyy = np.eye(3)

    ret = gauss_markov_fit(A, y, Cyy)
