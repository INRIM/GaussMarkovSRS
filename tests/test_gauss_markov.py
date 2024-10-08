import numpy as np
import pytest

from gauss_markov_srs.gauss_markov import gauss_markov_fit


def test_gauss_markov_fit():
    A = np.array([[1, 0], [1, 0.5], [1, 1]])
    y = np.array([0.0, 0.49, 1.01])
    Cyy = np.eye(3)

    ret = gauss_markov_fit(A, y, Cyy)

    assert ret.q == pytest.approx(np.array([-0.005, 1.01]))

    A = np.array([[1, 0], [1, 0.5], [1, 1], [1, 2]])
    y = np.array([0.0, 0.49, 1.01, 33.3])
    Cyy = np.eye(4)


def test_gauss_markov_fit_with_mask():
    A = np.array([[1, 0], [1, 0.5], [1, 1], [1, 2]])
    y = np.array([0.0, 0.49, 1.01, 33.3])
    Cyy = np.eye(4)

    ret = gauss_markov_fit(A, y, Cyy, mask=[1, 1, 1, 0])

    assert ret.q == pytest.approx(np.array([-0.005, 1.01]))
