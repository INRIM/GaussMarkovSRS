import decimal

import numpy as np

from gauss_markov_srs.gauss_markov import FitResult, FitStats
from gauss_markov_srs.postprocess import get_keys, results_to_decimal
from gauss_markov_srs.print import decimal_to_string, unc_to_precision
from gauss_markov_srs.utils import cov2corr


def save(basename, label, data, reference, corr_data, y, yu, fit_result: FitResult, fit_stats: FitStats):
    # whodid = """# File generated on {}
    # # Script used to generate this file: {}
    # # Command line used to generate this file:
    # # {}
    # #
    # """.format(datetime.datetime.now(), sys.argv[0], ' '.join(command))

    n_adj = fit_stats.n_adj
    keys = get_keys(data)

    # adj out
    with open(basename + ".adj", "w") as filo:
        dres, dresu = results_to_decimal(fit_result, reference)

        filo.write("# Adjusted frequency values\n")
        # filo.write(whodid)
        fmt = "# {:8s}\t" + "{: <20}\t{: <12}\t" + "{}\t" * 2 + "\n"
        filo.write(fmt.format("Atom", "Frequency", "Unc.", " q", " u"))

        for a, f, fu, q, qu in zip(reference["Atom"][1:], dres, dresu, fit_result.q, fit_result.qu):
            fmt = "{:8s}\t{: <20}\t{: <12}\t{: .2e}\t{: .2e}\n"
            precision = unc_to_precision(fu)
            filo.write(fmt.format(a, decimal_to_string(f, precision), decimal_to_string(fu, precision), q, qu))

    # full adj out
    with open(basename + ".fadj", "w") as filo:
        filo.write("# Adjusted frequency values\n")
        # filo.write(whodid)
        fmt = "# {:8s}\t" + "{: <20}\t{: <12}\t" + "\n"
        filo.write(fmt.format("Atom", "Frequency", "Unc."))

        for a, f, fu in zip(reference["Atom"][1:], dres, dresu):
            fmt = "{:8s}\t{: <20}\t{: <12}\n"
            filo.write(fmt.format(a, f, fu))

    # meas out
    with open(basename + ".ins", "w") as filo:
        measout = np.column_stack(
            (y, yu, fit_stats.relative_residuals, fit_stats.self_sensitivity, fit_result.weights.T)
        )

        filo.write(
            "# Frequency deviations, uncertainties, residuals, self-sensitivity coefficients and weights used for the adjusted frequencies\n"
        )
        # filo.write(whodid)
        fmt = "# {:24s}\t" + "{}\t" * 4 + "{:8s}\t" * n_adj + "\n"
        filo.write(fmt.format("Key", "y", "u", "Q/u", "Sc", *reference["Atom"][1:]))

        for k, d in zip(keys, measout):
            fmt = "{:24s}\t" + "{: .3e}\t" * 2 + "{: .3f}\t" * 2 + "{: .3f}\t" * (len(d) - 4) + "\n"
            filo.write(fmt.format(k, *d))

    # cov and cor out
    with open(basename + ".cov", "w") as filo:
        filo.write("# Covariance of adjusted frequencies\n")
        # filo.write(whodid)

        np.savetxt(filo, fit_result.Cqq)

    with open(basename + ".cor", "w") as filo:
        filo.write("# Correlation matrix of adjusted frequencies\n")
        # filo.write(whodid)
        filo.write("# For full numerical accuracy use the .cov file instead\n#\n")

        np.savetxt(filo, cov2corr(fit_result.Cqq), fmt="% .3f\t" * n_adj)

    # # stats
    # with open(basename + ".stats", "w") as filo:
    #     filo.write("# Statistical infos\n")
    #     filo.write(whodid)

    #     fmt = "{:24s}\t{} = {}\n"
    #     filo.write(fmt.format("Number of measurements", "Nmeas", nMeas))
    #     filo.write(fmt.format("Number of adjusted frequencies", "Nadj", n_adj))
    #     filo.write(fmt.format("Degrees of freedom", "Ndof", nDof))
    #     filo.write(fmt.format("Chi squared", "Chi2", chi2))
    #     filo.write(fmt.format("Birge ratio", "sqrt(Chi2/Ndof)", birge))
    #     filo.write(fmt.format("Birge ratio limit criteria", "1 + sqrt(2/nDof)", 1 + (2 / nDof) ** 0.5))
    #     filo.write(fmt.format("p-value", "P(Chi2<Chi2obs)", st.chi2(df=nDof).cdf(chi2)))
    #     filo.write("#\n")

    #     val, vec = np.linalg.eigh(Cyy)

    #     filo.write(fmt.format("Min eigenvalue of covariance inputs", "sqrt(Min(e))", np.amin(val) ** 0.5))
    #     filo.write(fmt.format("Number of non-zero correlations", "Nnz", np.count_nonzero(cyy - np.eye(nMeas)) / 2))

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
