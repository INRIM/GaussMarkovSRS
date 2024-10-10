import decimal

import numpy as np

from gauss_markov_srs.cross_index import get_cross_index
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


def calc_ratio(reference, fit_result, atom1, atom2):
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

    # # ratios out
    # with open(basename + ".rat", "w") as filo:
    #     filo.write("# Ratios between adjusted frequency values\n")
    #     #filo.write(whodid)
    #     fmt = "# {:8s}\t{:8s}\t" + "{: <20}\t{: <12}\t" + "{}\t" * 2 + "\n"
    #     filo.write(fmt.format("Atom A", "Atom B", "Ratio A/B", "Unc.", "r", "u"))

    #     # print ratios
    #     for j, J in enumerate(reference["Atom"][1:]):
    #         for i, I in enumerate(reference["Atom"][1:]):
    #             if i < j:
    #                 x = decimal.Decimal(fit_result.q[i] - fit_result.q[j])
    #                 u = (fit_result.Cqq[i, i] + fit_result.Cqq[j, j] - 2 * fit_result.Cqq[i, j]) ** 0.5

    #                 Dr = nu0[I] / nu0[J] * (decimal.Decimal(1) + x)
    #                 Du = nu0[I] / nu0[J] * decimal.Decimal(u)

    #                 # calculate decimal places for uncertainties with 3 digits
    #                 places = decimal.Decimal(10) ** (int(np.floor(np.log10(abs(float(Du))))) - 2)

    #                 # print("{:8}\t{:8}\t{: <20}\t{: <12}".format(I, J, Dr.quantize(places), Du.quantize(places)))

    #                 fmt = "{:8s}\t{:8s}\t" + "{: <20}\t{: <12}\t" + "{:.2e}\t" * 2 + "\n"

    #                 filo.write(fmt.format(I, J, Dr.quantize(places), Du.quantize(places), x, u))
