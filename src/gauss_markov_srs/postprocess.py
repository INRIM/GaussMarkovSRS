import decimal

import numpy as np

from gauss_markov_srs.gauss_markov import FitResult
from gauss_markov_srs.utils import cov2unc


def results_to_decimal(fit_result: FitResult, reference):
    # astype(Decimal) does not work
    dq = np.array([decimal.Decimal(x) for x in fit_result.q])
    dqu = np.array([decimal.Decimal(x) for x in fit_result.qu])
    ref = np.array([decimal.Decimal(x) for x in reference["nu0"][1:]])

    dres = ref * (dq + decimal.Decimal(1))
    dresu = ref * (dqu)

    return dres, dresu


# get a unique key for each measurement
def key(entry):
    s = "_".join((f'{entry["Id"]}', entry["Ref"], entry["Atom1"], entry["Atom2"], entry["Sup"]))
    return s.strip("_")  # get rid of last '_' is Sup is empty


def get_keys(data):
    return np.array([key(d) for d in data])
