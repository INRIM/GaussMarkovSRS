import numpy as np
import pytest
from hypothesis import given
from hypothesis import strategies as st

from gauss_markov_srs.cross_index import get_cross_index

string_strategy = st.lists(st.text(min_size=1, max_size=3), min_size=1, max_size=10, unique=True)


@given(strings=string_strategy)
def test_cross_index_inversion(strings):
    s_to_i, i_to_s = get_cross_index(strings)

    integers = s_to_i(strings)
    converted_strings = i_to_s(integers)

    assert np.array_equal(converted_strings, strings)


@given(strings=string_strategy)
def test_not_unique_strings(strings):
    with pytest.raises(ValueError):
        s_to_i, i_to_s = get_cross_index(strings * 2)


@given(strings=string_strategy)
def test_not_unique_ints(strings):
    if len(strings) > 1:
        with pytest.raises(ValueError):
            ones = np.ones_like(strings)
            s_to_i, i_to_s = get_cross_index(strings, ones)
    else:
        s_to_i, i_to_s = get_cross_index(strings)
