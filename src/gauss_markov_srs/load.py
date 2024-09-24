import decimal

import numpy as np
import numpy.lib.recfunctions as rfn
import pandas as pd


def load_excel(filename, dcols=[]):
    df = pd.read_excel(filename, na_filter=False, converters={d: decimal.Decimal for d in dcols}, dtype=None)
    df.columns.str.strip()

    records = df.to_records(index=False)

    return records


def load_txt(filename, dcols=[]):
    a = np.genfromtxt(filename, dtype=None, names=True, converters={d: decimal.Decimal for d in dcols}, encoding=None)
    return a


# get a unique key for each measurement
def _key(entry):
    s = "_".join((f'{entry["Id"]}', entry["Ref"], entry["Atom1"], entry["Atom2"], entry["Sup"]))
    return s.strip("_")  # get rid of last '_' is Sup is empty


def add_long_id_column(data):
    keys = np.array([_key(d) for d in data])
    b = rfn.append_fields(x, "Long_id", keys, usemask=False)
    return b


if __name__ == "__main__":
    x = load_excel("./Data/Input/CCTF2021/inputs.xlsx", dcols=[5, 6])
    x = add_long_id_column(x)
    y = load_txt("./Data/Reference/CCTF2021.txt", dcols=[2])
