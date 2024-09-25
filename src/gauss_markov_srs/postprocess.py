import decimal


def results_to_decimal(q, qu, reference):
    dq = q.astype(decimal.Decimal)
    dqu = qu.astype(decimal.Decimal)

    dres = reference["nu0"][1:] * (dq + 1)
    dresu = reference["nu0"][1:] * (dqu)

    return dres, dresu
