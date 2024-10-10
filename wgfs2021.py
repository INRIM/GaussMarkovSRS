from importlib import reload

import numpy as np
from scipy.optimize import brentq

import gauss_markov_srs as gm

gm = reload(gm)

input_data = gm.load_excel("./Data/Input/CCTF2021/inputs.xlsx")
input_corr = gm.load_txt("./Data/Input/CCTF2021/correlations.txt")
reference = gm.load_txt("./Data/Reference/CCTF2017.txt")


# input_data = gm.data_with_long_id(input_data)


y, yu = gm.get_y(input_data, reference)
M = gm.get_model_matrix(input_data, reference)
cyy = gm.get_corr_matrix(input_data, input_corr)

Cyy = gm.corr2cov(cyy, yu)


fit_results = gm.gauss_markov_fit(M, y, Cyy)
stats = gm.fit_stats(M, y, Cyy, fit_results)


gm.pretty_print_fit_result(fit_results, reference)
gm.pretty_print_stats(stats)

dir = "./Output/Adjustment-1"
figdir = "./Output/Adjustment-1/Figures"
label = "Adjustment-1"

gm.plot_residuals(figdir, label, input_data, reference, fit_results, stats)
gm.plot_correlation_matrix(figdir, label, cyy, figname="input-correlation.png")
gm.save(dir, label, reference, fit_results, stats, keys=gm.get_long_keys(input_data))


Sc = stats.self_sensitivity

mask = Sc > 0.01


fit_results2 = gm.gauss_markov_fit(M, y, Cyy, mask=mask)
stats2 = gm.fit_stats(M, y, Cyy, fit_results2, mask=mask)

print("\n====\n")
gm.pretty_print_fit_result(fit_results2, reference)
gm.pretty_print_stats(stats2)


dir = "./Output/Adjustment-2"
figdir = "./Output/Adjustment-2/Figures"
label = "Adjustment-2"

gm.plot_residuals(figdir, label, input_data, reference, fit_results2, stats2, mask=mask)
gm.plot_correlation_matrix(figdir, label, cyy, figname="input-correlation.png")
gm.save(dir, label, reference, fit_results2, stats2, mask=mask, keys=gm.get_long_keys(input_data))


# lets try a Random Effect Model
def random_effect_model(b):
    Cb = Cyy + b**2 * 1e-32 * np.eye(len(y))
    fit_results_temp = gm.gauss_markov_fit(M, y, Cb, mask=mask)
    stats_temp = gm.fit_stats(M, y, Cb, fit_results_temp, mask=mask)
    return stats_temp.chi2 / stats_temp.n_dof - 1.0


b = brentq(random_effect_model, 0, 5)
b *= 1e-16
Cb = Cyy + b**2 * np.eye(len(y))

fit_results3 = gm.gauss_markov_fit(M, y, Cb, mask=mask)
stats3 = gm.fit_stats(M, y, Cb, fit_results2, mask=mask)

print("\n====\n")
gm.pretty_print_fit_result(fit_results3, reference)
gm.pretty_print_stats(stats3)


gm.plot_fit_results(
    "./Output/", reference, [fit_results, fit_results2, fit_results3], ["Adjustment-1", "Adjustment-2", "Adjustment-3"]
)
