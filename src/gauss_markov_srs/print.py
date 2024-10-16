import decimal

import numpy as np
import scipy.stats as st

from gauss_markov_srs.gauss_markov import FitResult, FitStats
from gauss_markov_srs.postprocess import results_to_decimal


def unc_to_precision(u, digits=3):
    """Calculate the decimal places corresponding to a given uncertainty.

    Parameters
    ----------
    u : Decimal
        Uncertainty.
    digits : int, optional
        Numbe rof digits in the uncertainty, by default 3

    Returns
    -------
    Decimal
        Decimal precision to be passed to Decimal.quantize or `decimal_to_string`.
    """
    # calculate decimal places for uncertainties with 3 digits
    exp = decimal.Decimal(10) ** (int(np.floor(np.log10(abs(float(u))))) - digits + 1)
    return exp


def decimal_to_string(d, precision):
    """Convert a Decimal to string with a given precision.

    Parameters
    ----------
    d : Decimal
        Number to be converted.
    precision : Decimal
        Precision to be passed to Decimal.quantize

    Returns
    -------
    string
        Input number as string.
    """
    return f"{d.quantize(precision)}".lower()


def pretty_print_fit_result(fit_result: FitResult, reference):
    """Pretty print fit results.

    Parameters
    ----------
    reference : structured array.
        Input reference.
    fit_result : FitResult
        Dataclass of fit results from `gauss_markov_fit`
    """
    # print adjusted constants
    # also save adjusted values as string
    # print("Transition\tvalue\t\t\tunc\trel.unc")

    dres, dresu = results_to_decimal(fit_result, reference)

    print("{}\t{: <20}\t{: <12}\t{}".format("Transition", "Value", "Unc.", "Rel. Unc."))
    for i, x, u in zip(reference["Atom"][1:], dres, dresu):
        precision = unc_to_precision(u)

        print(
            "{:8}\t{: <20}\t{: <12}\t{:.1e}".format(
                i, decimal_to_string(x, precision), decimal_to_string(u, precision), u / x
            )
        )


def pretty_print_stats(fit_stats: FitStats):
    """Pretty print fit stats.

    Parameters
    ----------
    fit_stats : FitStats
        Dataclass of fit stats from `fit_stats`
    """
    print(f"Number of measurements = {fit_stats.n_meas}")
    print(f"Number of adjusted constants = {fit_stats.n_adj}")
    if fit_stats.n_excluded > 0:
        print(f"Number of measurements excluded from the fit = {fit_stats.n_excluded}")
        print(f"Number of measurements included in the fit = {fit_stats.n_included}")
    print(f"Number of degrees of freedom = {fit_stats.n_dof}")
    print(f"Chi squared = {fit_stats.chi2:.3f}")
    print(f"Birge ratio = {fit_stats.birge:.3f}")
    print(f"Birge ratio limit = 1 + sqrt(2/N_dof) = {fit_stats.birge_limit:.3f}")
    print(f"p value = prob(chi2 < chi2 observed) = {fit_stats.p_value:.3f}")
    print(f"Inputs > 2 sigmas = {np.sum(np.abs(fit_stats.relative_residuals)>2)}")
    print(f"Inputs Sc > 1% = {np.sum(fit_stats.self_sensitivity >= 0.01)}")
    print(f"Number of non-zero correlations = {fit_stats.n_corr_nonzero}")
    print(f"Square root of min eigenvalue of covariance inputs = {fit_stats.min_eigenvalue:.3e}")


if __name__ == "__main__":
    x = unc_to_precision(decimal.Decimal(1.42123), 3)
    print(decimal_to_string(decimal.Decimal("1000.123456789"), x))
