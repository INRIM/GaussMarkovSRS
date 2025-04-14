import re

import numpy as np

from gauss_markov_srs.cross_index import get_cross_index
from gauss_markov_srs.load import load_txt


def parse_wgfs_file(filename, reference=None):
    """Parse a file with the format shared within the WGFS returning structured arrays ready for fitting.

    Parameters
    ----------
    filename : str
        Input file
    reference : structured array, optional
        convert the names of Atom using the array as a reference, by default None
        This structured array should be in reference format and should match Id with the convention used in the input file (e.g., "nu1", "nu2", etc..).

    Returns
    -------
    ndarray, ndarray
        structured arrays of input data and input correlations
        defining the columns expected by gauss_markov_srs, with the type of each column as 'o' for unicode strings
    """
    with open(filename, "r") as f:
        lines = f.readlines()

    reference_lines = []
    input_lines = []
    corr_lines = []

    current_reading = None
    for line in lines:
        line = line.strip()
        if not line:
            continue
        elif re.match(r"^M\s+\d+$", line):
            current_reading = "reference"
            continue
        elif re.match(r"^N\s+\d+$", line):
            current_reading = "input"
            continue
        elif re.match(r"^NCorrCoeffs\s+\d+$", line):
            current_reading = "corr"
            continue
        if current_reading == "reference":
            reference_lines.append(line)
        if current_reading == "input":
            input_lines.append(line)
        if current_reading == "corr":
            corr_lines.append(line)

    dtype_input = [
        ("Id", "O"),
        ("Ref", "O"),
        ("Atom1", "O"),
        ("Atom2", "O"),
        ("Sup", "O"),
        ("Value", "O"),
        ("Unc", "O"),
        ("Note", "O"),
    ]
    input_array = np.array([parse_wgfs_input_line(line, reference) for line in input_lines], dtype=dtype_input)

    dtype_corr = [("Id1", "O"), ("Id2", "O"), ("Corr", "O")]
    corr_array = np.array([parse_wgfs_corr_line(line) for line in corr_lines], dtype=dtype_corr)

    return input_array, corr_array


def parse_wgfs_input_line(line, reference=None):
    """Parse a single line describing an input measurement from the format shared within the WGFS.

    Parameters
    ----------
    line : str
        line to be parsed
    reference : structured array, optional
        convert the names of Atom using the array as a reference, by default None
        This structured array should be in reference format and should match Id with the convention used in the input file (e.g., "nu1", "nu2", etc..).

    Returns
    -------
    tuple
        tuple of (id, ref, atom1, atom2, sup, value, unc, note)

    """
    parts = line.split(None, 4)

    id = parts[0]
    ref, note = (parts[4].strip("[]").split(",", 1) + [""])[:2]
    ref, sup = split_name_year_suffix(ref)
    ref1, ref2 = (parts[1].split("_over_", 1) + ["nu0"])[:2]
    value = parts[2]
    unc = parts[3]

    if reference is None:
        atom1, atom2 = ref1, ref2
    else:
        ref_str_to_i, i_to_ref_str = get_cross_index(reference["Id"])

        atom1 = reference["Atom"][ref_str_to_i(ref1)]
        atom2 = reference["Atom"][ref_str_to_i(ref2)]

    return (id, ref, atom1, atom2, sup, value, unc, note)


def split_name_year_suffix(s):
    """Parse the reference from the format shared within the WGFS

    Parameters
    ----------
    s : str
        reference

    Returns
    -------
    tuple
        tuple of (ref. as name+year, supplementary info)

    Examples
    --------
    "brown2022"	-> ("brown2022", "")
    "smith2024a" -> ("smith2024a","")
    "charles2023errata" ->	("charles2023", "errata")

    """
    match = re.match(r"^([a-z]+)(\d{4})([a-z]*)$", s, re.IGNORECASE)
    if not match:
        return s, ""

    name, year, suffix = match.groups()

    if suffix == "":
        return name + year, ""
    elif len(suffix) == 1:
        return name + year + suffix, ""
    else:
        return name + year, suffix


def parse_wgfs_corr_line(line):
    """Parse a line of correlations from the format shared within the WGFS.

    Parameters
    ----------
    line : str
        line to be read

    Returns
    -------
    tuple
        tuple of (Id1, id2, corr)

    Examplea
    -------
    "r(q3,q71)	0.004" -> ("q3", "q71", "0.004")

    """
    match = re.match(r"r\(([^,]+),([^)]+)\)\s+(.+)", line.strip())
    if not match:
        raise ValueError(f"Line does not match expected format: {line}")
    id1, id2, value = match.groups()
    return (id1.strip(), id2.strip(), value.strip())


if __name__ == "__main__":
    fili = "../../Data/Live/Metrologia_paper_secondary/ClockInputData2020_7.dat"
    refi = "../../Data/Reference/CCTF2021.txt"

    ref = load_txt(refi)

    a, ac = parse_wgfs_file(fili)
    b, bc = parse_wgfs_file(fili, ref)
