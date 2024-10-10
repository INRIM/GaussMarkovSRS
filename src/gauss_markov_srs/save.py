import decimal
import os

import numpy as np

from gauss_markov_srs.gauss_markov import FitResult, FitStats
from gauss_markov_srs.postprocess import calc_ratio, get_keys, results_to_decimal
from gauss_markov_srs.print import decimal_to_string, unc_to_precision
from gauss_markov_srs.utils import cov2corr


def save(dir, label, data, reference, corr_data, y, yu, fit_result: FitResult, fit_stats: FitStats):
    # whodid = """# File generated on {}
    # # Script used to generate this file: {}
    # # Command line used to generate this file:
    # # {}
    # #
    # """.format(datetime.datetime.now(), sys.argv[0], ' '.join(command))

    os.makedirs(dir, exist_ok=True)

    n_adj = fit_stats.n_adj
    keys = get_keys(data)

    label = os.path.join(dir, label)

    # adj out
    with open(label + "-adj.txt", "w") as filo:
        dres, dresu = results_to_decimal(fit_result, reference)

        filo.write("# Adjusted frequency values\n")
        # filo.write(whodid)
        fmt = "{:8s}\t" + "{: <20}\t{: <12}\t" + "{}\t" * 2 + "\n"
        filo.write(fmt.format("# Atom", "Frequency", "Unc.", " q", " u"))

        for a, f, fu, q, qu in zip(reference["Atom"][1:], dres, dresu, fit_result.q, fit_result.qu):
            fmt = "{:8s}\t{: <20}\t{: <12}\t{: .2e}\t{: .2e}\n"
            precision = unc_to_precision(fu)
            filo.write(fmt.format(a, decimal_to_string(f, precision), decimal_to_string(fu, precision), q, qu))

    # full adj out
    with open(label + "-adj-full.txt", "w") as filo:
        filo.write("# Adjusted frequency values\n")
        # filo.write(whodid)
        fmt = "{:8s}\t" + "{: <20}\t{: <12}\t" + "\n"
        filo.write(fmt.format("# Atom", "Frequency", "Unc."))

        for a, f, fu in zip(reference["Atom"][1:], dres, dresu):
            fmt = "{:8s}\t{: <20}\t{: <12}\n"
            filo.write(fmt.format(a, f, fu))

    # meas out
    with open(label + "-inputs.txt", "w") as filo:
        measout = np.column_stack(
            (y, yu, fit_stats.relative_residuals, fit_stats.self_sensitivity, fit_result.weights.T)
        )

        filo.write(
            "# Frequency deviations, uncertainties, residuals, self-sensitivity coefficients and weights used for the adjusted frequencies\n"
        )
        # filo.write(whodid)
        fmt = "{:42s}\t" + "{:8s}\t" * 2 + "{:5s}\t" * 2 + "{:8s}\t" * n_adj + "\n"
        filo.write(fmt.format("# Key", "y", "u", "Q/u", "Sc", *reference["Atom"][1:]))

        for k, d in zip(keys, measout):
            fmt = "{:42s}\t" + "{:.3e}\t" * 2 + "{:.3f}\t" * 2 + "{:.3f}   \t" * (len(d) - 4) + "\n"
            filo.write(fmt.format(k, *d))

    # cov and cor out
    with open(label + "-covariance.txt", "w") as filo:
        filo.write("# Covariance of adjusted frequencies\n")
        # filo.write(whodid)

        np.savetxt(filo, fit_result.Cqq)

    with open(label + "-correlation.txt", "w") as filo:
        filo.write("# Correlation matrix of adjusted frequencies\n")
        # filo.write(whodid)
        filo.write("# For full numerical accuracy use the covariance file instead\n#\n")

        np.savetxt(filo, cov2corr(fit_result.Cqq), fmt="% .3f\t" * n_adj)

    # stats
    with open(label + "-stats.txt", "w") as filo:
        filo.write("# Statistical infos\n")

        filo.write(f"Number of measurements = {fit_stats.n_meas}" + "\n")
        filo.write(f"Number of adjusted constants = {fit_stats.n_adj}" + "\n")
        filo.write(f"Number of measurements excluded from the fit = {fit_stats.n_excluded}" + "\n")
        filo.write(f"Number of measurements included in the fit = {fit_stats.n_included}" + "\n")
        filo.write(f"Number of degrees of freedom = {fit_stats.n_dof}" + "\n")
        filo.write(f"Chi squared = {fit_stats.chi2:.3f}" + "\n")
        filo.write(f"Birge ratio = {fit_stats.birge:.3f}" + "\n")
        filo.write(f"Birge ratio limit = 1 + sqrt(2/N_dof) = {fit_stats.birge_limit:.3f}" + "\n")
        filo.write(f"p value = prob(chi2 < chi2 observed) = {fit_stats.p_value:.3f}" + "\n")
        filo.write(f"Inputs > 2 sigmas = {np.sum(np.abs(fit_stats.relative_residuals)>2)}" + "\n")
        filo.write(f"Inputs Sc > 1% = {np.sum(fit_stats.self_sensitivity >= 0.01)}" + "\n")
        filo.write(f"Number of non-zero correlations = {fit_stats.n_corr_nonzero}" + "\n")
        filo.write(f"Square root of min eigenvalue of covariance inputs = {fit_stats.min_eigenvalue:.3e}" + "\n")

    # ratios out
    with open(label + "-ratios.txt", "w") as filo:
        filo.write("# Ratios between adjusted frequency values\n")
        # filo.write(whodid)
        fmt = "{:8s}\t{:8s}\t" + "{: <20}\t{: <12}\t" + "\n"
        filo.write(fmt.format("# Atom A", "Atom B", "Ratio A/B", "Unc."))

        # print ratios
        for i, I in enumerate(reference["Atom"][1:]):
            for j, J in enumerate(reference["Atom"][1:]):
                if i < j:
                    Dr, Du = calc_ratio(reference, fit_result, I, J)

                    precision = unc_to_precision(Du)

                    fmt = "{:8s}\t{:8s}\t" + "{: <20}\t{: <12}\t" + "\n"

                    filo.write(fmt.format(I, J, decimal_to_string(Dr, precision), decimal_to_string(Du, precision)))
