import numpy as np
import pandas as pd


def load_excel(filename):
    """Load an excel file as a Numpy structured array.
    All fields are loaded as strings.
    It is based on pandas `read_excel`.

    Parameters
    ----------
    filename : string
        Excel file to be open.

    Returns
    -------
    ndarray
        structured array where the name of the columns is derived by the first row, and the type of each column is 'o' for unicode strings.
    """
    df = pd.read_excel(filename, na_filter=False, dtype=str)

    df.columns = df.columns.str.strip()

    records = df.to_records(index=False)

    return records


def load_txt(filename):
    """Load a space or tab separated txt file as a Numpy structured array.
    All fields are loaded as strings.
    It is based on pandas `read_csv`.

    Parameters
    ----------
    filename : string
        Text file to be open.

    Returns
    -------
    ndarray
        structured array where the name of the columns is derived by the first row, and the type of each column is 'o' for unicode strings.
    """
    # a = np.genfromtxt(filename, dtype=str, names=True)
    df = pd.read_csv(filename, na_filter=False, dtype=str, sep="\s+")

    df.columns = df.columns.str.strip("# ")

    return df.to_records(index=False)



if __name__ == "__main__":
    x = load_excel("./Data/Input/CCTF2021/inputs.xlsx")
    y = load_txt("./Data/Reference/CCTF2021.txt")
    z = load_txt("./Data/Input/CCTF2021/correlations.txt")
