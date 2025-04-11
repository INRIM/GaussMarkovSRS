import numpy as np
import re


from gauss_markov_srs.cross_index import get_cross_index
from gauss_markov_srs.load import load_txt

def parse_wgfs_file(filename, reference=None):
    with open(filename, "r") as f:
        lines = f.readlines()

    reference_lines = []
    input_data_lines = []
    input_corr_lines = []

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
            input_data_lines.append(line)
        if current_reading == "corr":
            input_corr_lines.append(line)

    dtype_input = [("Id", "O"), ("Ref", "O"), ("Atom1", "O"), ("Atom2", "O"), ("Sup", "O"), ("Value", "O"), ("Unc", "O"), ("Note", "O")]
    arr2 = np.array([parse_wgfs_input(line, reference) for line in input_data_lines], dtype=dtype_input)

    return arr2

def parse_wgfs_input(line, reference = None):
    parts = line.split(None, 4)

    id = parts[0]
    ref, note = (parts[4].strip("[]").split(',' ,1)+[""])[:2]
    ref, sup = split_name_year_suffix(ref)
    ref1 ,ref2 = (parts[1].split('_over_' ,1)+["nu0"])[:2]
    value = parts[2]
    unc = parts[2]

    if reference is None:
        atom1, atom2 = ref1, ref2
    else:
        ref_str_to_i, i_to_ref_str = get_cross_index(reference["Id"])

        atom1 = reference["Atom"][ref_str_to_i(ref1)]
        atom2 = reference["Atom"][ref_str_to_i(ref1)]

    return (id, ref, atom1, atom2, sup, value, unc, note)



def split_name_year_suffix(s):
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
    

if __name__ == "__main__":
    fili = "../../Data/Live/Metrologia_paper_secondary/ClockInputData2020_7.dat"
    refi = "../../Data/Reference/CCTF2021.txt"

    ref = load_txt(refi)

    a = parse_wgfs_file(fili)
    b = parse_wgfs_file(fili, ref)
