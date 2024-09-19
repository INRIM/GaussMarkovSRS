"""
This submodule implements conversions from integer indices to string indices and vice versa.
Integer indices are easier code wise while string indices are more human readables.
"""

import numpy as np


def get_cross_index(strings, integers=None):
    if integers is None:
        integers = np.arange(len(strings))

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
