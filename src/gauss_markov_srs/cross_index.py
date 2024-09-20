"""
This submodule implements conversions from integer indices to string indices and vice versa.
Integer indices are easier code wise while string indices are more human readables.
"""

import numpy as np


def get_cross_index(strings, integers=None):
    """Return two vectorized functions to convert integer or string indices.

    Parameters
    ----------
    strings : array_like
        List of unique identifier.
    integers : array_like, optional
        List of unique identifiers, by default None (converts to np.arange(len(strings))).

    Returns
    -------
    callable, callable
        Vectorized function to convert from string to integers and from integers to strings.

    Raises
    ------
    ValueError
        If strings or integers are not unique or of different lengths.
    """

    unique, counts = np.unique(strings, return_counts=True)
    if (counts > 1).any():
        raise ValueError(f"Strings must be unique. Not unique strings = {unique[counts>1]}")

    if integers is None:
        integers = np.arange(len(strings))
    else:
        unique, counts = np.unique(integers, return_counts=True)
        if (counts > 1).any():
            raise ValueError(f"Integers must be unique. Not unique strings = {unique[counts>1]}")

    if len(strings) != len(integers):
        raise ValueError("String and integer indices must have the same number of elements.")

    string_to_integer = {s: i for s, i in zip(strings, integers)}
    integer_to_string = {i: s for s, i in zip(strings, integers)}
    s_to_i = np.vectorize(string_to_integer.get)
    i_to_s = np.vectorize(integer_to_string.get)

    return s_to_i, i_to_s


if __name__ == "__main__":
    strings = ["A", "B", "C"]

    s_to_i, i_to_s = get_cross_index(strings)
