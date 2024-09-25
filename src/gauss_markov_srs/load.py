import decimal

import numpy as np
import pandas as pd


def load_excel(filename, dcols=[]):
    df = pd.read_excel(filename, na_filter=False, converters={d: decimal.Decimal for d in dcols}, dtype=None)
    df.columns.str.strip()

    records = df.to_records(index=False)

    return records


def load_txt(filename, dcols=[]):
    a = np.genfromtxt(filename, dtype=None, names=True, converters={d: decimal.Decimal for d in dcols}, encoding=None)
    return a


if __name__ == "__main__":
    x = load_excel("./Data/Input/CCTF2021/inputs.xlsx", dcols=[5, 6])
    y = load_txt("./Data/Reference/CCTF2021.txt", dcols=[2])
    z = load_txt("./Data/Input/CCTF2021/correlations.txt")
