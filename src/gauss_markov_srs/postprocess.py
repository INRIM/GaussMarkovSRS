import decimal

import numpy as np


def results_to_decimal(q, qu, reference):
    # astype(Decimal) does not work
    dq = np.array([decimal.Decimal(x) for x in q])
    dqu = np.array([decimal.Decimal(x) for x in qu])

    dres = reference["nu0"][1:] * (dq + decimal.Decimal(1))
    dresu = reference["nu0"][1:] * (dqu)

    return dres, dresu
