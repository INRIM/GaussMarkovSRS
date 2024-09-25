import decimal

import numpy as np
import scipy.stats as st

from gauss_markov_srs.postprocess import results_to_decimal
from gauss_markov_srs.utils import cov2unc


def unc_to_precision(u, digits=3):
    # calculate decimal places for uncertainties with 3 digits
    exp = decimal.Decimal(10) ** (int(np.floor(np.log10(abs(float(u))))) - digits + 1)
    return exp


def decimal_to_string(d, precision):
    return f"{d.quantize(precision)}"


def pretty_print_fit_result(q, Cqq, weights, reference):
    # print adjusted constants
    # also save adjusted values as string
    # print("Transition\tvalue\t\t\tunc\trel.unc")

    qu = cov2unc(Cqq)
    dres, dresu = results_to_decimal(q, qu, reference)

    print("{}\t{: <20}\t{: <12}\t{}".format("Transition", "Value", "Unc.", "Rel. Unc."))
    for i, x, u in zip(reference["Atom"][1:], dres, dresu):
        precision = unc_to_precision(u)

        print(
            "{:8}\t{: <20}\t{: <12}\t{:.1e}".format(
                i, decimal_to_string(x, precision), decimal_to_string(u, precision), u / x
            )
        )


def pretty_print_stats(Sc, Q, nDof, chi2):
    birge = (chi2 / (nDof)) ** 0.5

    print("Birge ratio =", birge)
    print("Birge ratio limit = 1 + sqrt(2/nDof) =", 1 + (2 / nDof) ** 0.5)
    print("p value =", st.chi2(df=nDof).cdf(chi2))
    # print("Inputs > 2 sigmas =", sum(abs(Q/yu)>2))
    # print("Inputs  = ", nMeas)
    print("Inputs Sc > 1% = ", sum(Sc >= 0.01))


if __name__ == "__main__":
    x = unc_to_precision(decimal.Decimal(1.42123), 3)
    print(decimal_to_string(decimal.Decimal("1000.123456789"), x))
