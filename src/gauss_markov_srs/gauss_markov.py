from dataclasses import dataclass

import numpy as np

# @dataclass
# class FitResult:
#     q: np.ndarray
#     cov: np.ndarray
#     weights: np.ndarray


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

    invCyy = np.linalg.inv(Cyy)

    Cqq = np.linalg.inv(A.transpose() @ invCyy @ A)
    weights = Cqq @ A.transpose() @ invCyy

    q = np.dot(weights, y).T

    # nm, na = A.shape

    ret_q = np.squeeze(np.array(q, ndmin=1))
    ret_C = np.array(Cqq)
    ret_weights = np.squeeze(np.array(weights, ndmin=1))

    ret = (ret_q, ret_C, ret_weights)
    return ret

    # self_sensitivity: np.ndarray
    # residuals: np.ndarray
    # ndof: int
    # chi2: float

    # #         self_sensitivity=self_sensitivity,
    #     residuals=residuals,
    #     ndof=ndof,
    #     chi2=chi2,

    #     self_sensitivity = np.squeeze(np.array(Sc, ndmin=1))
    # residuals = np.squeeze(np.array(Q, ndmin=1))
    # ndof = nm - na
    #  #     Q = y - A @ q
    # chi2 = Q.T @ invCyy @ Q

    # Sc = np.diag(A @ weights)
    # self_sensitivity : ndarray, shape (N,)
    #     Self-sensitivity coefficients for the input data.
    # residuals : ndarray, shape (N, )
    #     Residuals of the fit.
    # ndof : int
    #     Number of degrees of freedom.
    # chi2 : float
    #     Chi squared of the fit


if __name__ == "__main__":
    A = np.array([[1, 0], [1, 0.5], [1, 1]])
    y = np.array([0.0, 0.49, 1.01])
    Cyy = np.eye(3)

    ret = gauss_markov_fit(A, y, Cyy)
