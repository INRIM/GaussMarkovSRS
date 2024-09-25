import decimal

import numpy as np

from gauss_markov_srs.gauss_markov import FitResult
from gauss_markov_srs.utils import cov2unc


def results_to_decimal(fit_result: FitResult, reference):
    # astype(Decimal) does not work
    dq = np.array([decimal.Decimal(x) for x in fit_result.q])
    dqu = np.array([decimal.Decimal(x) for x in fit_result.qu])

    dres = reference["nu0"][1:] * (dq + decimal.Decimal(1))
    dresu = reference["nu0"][1:] * (dqu)

    return dres, dresu
