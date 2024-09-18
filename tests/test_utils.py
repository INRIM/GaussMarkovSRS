import numpy as np
import pytest
from hypothesis import given
from hypothesis import strategies as st
from hypothesis.extra.numpy import arrays

from gauss_markov_srs.utils import corr2cov, cov2corr


@given(arrays(np.float64, (5, 5), elements=st.floats(min_value=1.0, max_value=10.0)))
def test_cov2corr_and_corr2cov_inversion(cov):
    # Ensure cov is symmetric and positive definite
    cov = np.dot(cov, cov.T)

    corr, std = cov2corr(cov, return_std=True)
    cov_reconstructed = corr2cov(corr, std)

    assert cov == pytest.approx(cov_reconstructed)
