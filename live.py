import os
from importlib import reload

import numpy as np

import gauss_markov_srs as gm

gm = reload(gm)

reference = gm.load_txt("./Data/Reference/CCTF2021.txt")
input_data, input_corr = gm.parse_wgfs_file("./Data/Live/ClockInputData2025Reduced.dat", reference)


basedir = "./Output-Live/"
labels = []


y, yu = gm.get_y(input_data, reference)
M = gm.get_model_matrix(input_data, reference)
cyy = gm.get_corr_matrix(input_data, input_corr)

Cyy = gm.corr2cov(cyy, yu)


# Adjustment 1
label = "Adjustment-1"
labels += [label]
dir = os.path.join(basedir, label)
figdir = os.path.join(dir, "Figures")


fit_results = gm.gauss_markov_fit(M, y, Cyy)
stats = gm.fit_stats(M, y, Cyy, fit_results)


gm.pretty_print_fit_result(fit_results, reference)
gm.pretty_print_stats(stats)


gm.plot_residuals(figdir, label, input_data, reference, fit_results, stats)
gm.plot_correlation_matrix(figdir, label, cyy, figname="input-correlation.png")
gm.save(dir, label, reference, fit_results, stats, keys=gm.get_long_keys(input_data))


gm.plot_ratio_residuals(
    figdir, ["27Al+", "171Yb+E3", "171Yb", "87Sr"], label, input_data, reference, fit_results, stats
)

# Adjustment 2
# remove unpublished and problematic from ROCIT
label = "Adjustment-2"
labels += [label]
dir = os.path.join(basedir, label)
figdir = os.path.join(dir, "Figures")


unpublished = [107, 108, 110, 112, 113, 114, 130, 131, 132, 147]
rocit_not_to_be_used = [121, 124, 125, 127, 138, 157, 159, 160]


def contains_any_integer(strings, integers):
    int_strs = [str(i) for i in integers]
    return np.array([any(s_int in s for s_int in int_strs) for s in strings])


mask = ~contains_any_integer(input_data["Id"], unpublished + rocit_not_to_be_used)

fit_results2 = gm.gauss_markov_fit(M, y, Cyy, mask=mask)
stats2 = gm.fit_stats(M, y, Cyy, fit_results2, mask=mask)

print("\n====\n")
gm.pretty_print_fit_result(fit_results2, reference)
gm.pretty_print_stats(stats2)

gm.plot_residuals(figdir, label, input_data, reference, fit_results2, stats2, mask=mask)
gm.plot_correlation_matrix(figdir, label, cyy, figname="input-correlation.png")
gm.save(dir, label, reference, fit_results2, stats2, mask=mask, keys=gm.get_long_keys(input_data))

gm.plot_ratio_residuals(
    figdir,
    ["27Al+", "171Yb+E3", "171Yb", "87Sr"],
    label,
    input_data,
    reference,
    fit_results2,
    stats2,
    mask=mask,
)

# Adjustment 3
# increase unc of some ratios
label = "Adjustment-3"
labels += [label]
dir = os.path.join(basedir, label)
figdir = os.path.join(dir, "Figures")

ratio_to_increase_unc = [92, 102, 103, 104, 120, 137]
unc_mask = contains_any_integer(input_data["Id"], ratio_to_increase_unc)
yu3 = yu
yu3[unc_mask] *= 2
Cyy3 = gm.corr2cov(cyy, yu3)

fit_results3 = gm.gauss_markov_fit(M, y, Cyy3, mask=mask)
stats3 = gm.fit_stats(M, y, Cyy3, fit_results3, mask=mask)

print("\n====\n")
gm.pretty_print_fit_result(fit_results3, reference)
gm.pretty_print_stats(stats3)

gm.plot_residuals(figdir, label, input_data, reference, fit_results3, stats3, mask=mask)
gm.plot_correlation_matrix(figdir, label, cyy, figname="input-correlation.png")
gm.save(dir, label, reference, fit_results3, stats3, mask=mask, keys=gm.get_long_keys(input_data))

gm.plot_ratio_residuals(
    figdir,
    ["27Al+", "171Yb+E3", "171Yb", "87Sr"],
    label,
    input_data,
    reference,
    fit_results3,
    stats3,
    mask=mask,
)


# Adjustment 4
# reomve Sc <0.01
label = "Adjustment-4"
labels += [label]
dir = os.path.join(basedir, label)
figdir = os.path.join(dir, "Figures")

# should include adj2 exclusions
mask = stats3.self_sensitivity > 0.01

fit_results4 = gm.gauss_markov_fit(M, y, Cyy3, mask=mask)
stats4 = gm.fit_stats(M, y, Cyy3, fit_results4, mask=mask)

print("\n====\n")
gm.pretty_print_fit_result(fit_results4, reference)
gm.pretty_print_stats(stats4)

gm.plot_residuals(figdir, label, input_data, reference, fit_results4, stats4, mask=mask)
gm.plot_correlation_matrix(figdir, label, cyy, figname="input-correlation.png")
gm.save(dir, label, reference, fit_results4, stats4, mask=mask, keys=gm.get_long_keys(input_data))

gm.plot_ratio_residuals(
    figdir,
    ["27Al+", "171Yb+E3", "171Yb", "87Sr"],
    label,
    input_data,
    reference,
    fit_results4,
    stats4,
    mask=mask,
)


gm.plot_fit_results(basedir, reference, [fit_results, fit_results2, fit_results3, fit_results4], labels)
