import decimal

import numpy as np

from gauss_markov_srs.cross_index import get_cross_index
from gauss_markov_srs.gauss_markov import FitResult


def results_to_decimal(fit_result: FitResult, reference):
    """Cast fit results as arbitrary precision numbers.

    Parameters
    ----------

    fit_result : FitResult
        Dataclass of fit results from `gauss_markov_fit`

    Returns
    -------
    array, array
        adjusted constants and uncertainty as arrays of decimal.Decimal
    """
    # astype(Decimal) does not work
    dq = np.array([decimal.Decimal(x) for x in fit_result.q])
    dqu = np.array([decimal.Decimal(x) for x in fit_result.qu])
    ref = np.array([decimal.Decimal(x) for x in reference["nu0"][1:]])

    dres = ref * (dq + decimal.Decimal(1))
    dresu = ref * (dqu)

    return dres, dresu


# get a unique key for each measurement
def _key(entry):
    """Generate a human-readable key from an entry in the input data structured array.

    Returns
    -------
    string
        Key
    """
    s = "_".join((f'{entry["Id"]}', entry["Ref"], entry["Atom1"], entry["Atom2"], entry["Sup"]))
    return s.strip("_")  # get rid of last '_' is Sup is empty


def get_long_keys(data):
    """Generate human-readable keys from the input data structured array.

    Parameters
    ----------
    data : structured array
        Input data.

    Returns
    -------
    array of strings
        Long keys for the input data
    """
    return np.array([_key(d) for d in data])


def calc_ratio(reference, fit_result, atom1, atom2):
    """Claculate a ratio of adjusted frequencies.

    Parameters
    ----------
    reference : structured array.
        Input reference.
    fit_result : FitResult
        Dataclass of fit results from `gauss_markov_fit`
    atom1 : string
        string reference for atom1
    atom2 : string
        string reference for atom2

    Returns
    -------
    Decimal, Decimal
        value and uncertainty of the ratio Atom1 / Atom2
    """
    ref_str_to_i, i_to_ref_str = get_cross_index(reference["Atom"])

    i1 = ref_str_to_i(atom1)
    i2 = ref_str_to_i(atom2)
    x = decimal.Decimal(fit_result.q[i1 - 1] - fit_result.q[i2 - 1])
    u = (fit_result.Cqq[i1 - 1, i1 - 1] + fit_result.Cqq[i2 - 1, i2 - 1] - 2 * fit_result.Cqq[i1 - 1, i2 - 1]) ** 0.5

    nu_1 = decimal.Decimal(reference[i1]["nu0"])
    nu_2 = decimal.Decimal(reference[i2]["nu0"])

    Dr = nu_1 / nu_2 * (decimal.Decimal(1) + x)
    Du = nu_1 / nu_2 * decimal.Decimal(u)

    return Dr, Du
