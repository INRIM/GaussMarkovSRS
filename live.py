import os
from importlib import reload

import numpy as np

import gauss_markov_srs as gm

gm = reload(gm)

reference = gm.load_txt("./Data/Reference/CCTF2021.txt")

input_dir = "./Data/Live/August/"
base_dir = "./Output-Live/"
filis = [
    "ClockInputData2025_mod_basis_2025TAICorr.dat",
    "ClockInputData2025_mod_basis_2025TAICorr_n=153.dat",
    "ClockInputData2025_mod_basis_2025TAIROCITCorr.dat",
    "ClockInputData2025_mod_basis_2025TAIROCITCorr_n=153.dat",
]

results = []
labels = []

for file in filis:
    input_data, input_corr = gm.parse_wgfs_file(os.path.join(input_dir, file), reference)

    y, yu = gm.get_y(input_data, reference)
    M = gm.get_model_matrix(input_data, reference)
    cyy = gm.get_corr_matrix(input_data, input_corr)

    Cyy = gm.corr2cov(cyy, yu)

    # Adjustment 1
    label = file
    dir = os.path.join(base_dir, label)
    figdir = os.path.join(dir, "Figures")

    fit_results = gm.gauss_markov_fit(M, y, Cyy)
    stats = gm.fit_stats(M, y, Cyy, fit_results)

    print("\n\n")
    print(file + "\n")
    gm.pretty_print_fit_result(fit_results, reference)
    gm.pretty_print_stats(stats)

    gm.plot_residuals(figdir, label, input_data, reference, fit_results, stats)
    gm.plot_correlation_matrix(figdir, label, cyy, figname="input-correlation.png")
    gm.save(dir, label, reference, fit_results, stats, keys=gm.get_long_keys(input_data))

    gm.plot_ratio_residuals(
        figdir, ["27Al+", "171Yb+E3", "171Yb", "87Sr"], label, input_data, reference, fit_results, stats
    )

    labels += [label]
    results += [fit_results]


gm.plot_fit_results(base_dir, reference, results, filis)
