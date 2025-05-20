from importlib import reload

import numpy as np

import gauss_markov_srs as gm

gm = reload(gm)

reference = gm.load_txt("./Data/Reference/CCTF2021.txt")
input_data, input_corr = gm.parse_wgfs_file(
    "./Data/Live/ClockInputData2025Reduced.dat", reference
)

y, yu = gm.get_y(input_data, reference)
M = gm.get_model_matrix(input_data, reference)
cyy = gm.get_corr_matrix(input_data, input_corr)

Cyy = gm.corr2cov(cyy, yu)


fit_results = gm.gauss_markov_fit(M, y, Cyy)
stats = gm.fit_stats(M, y, Cyy, fit_results)


gm.pretty_print_fit_result(fit_results, reference)
gm.pretty_print_stats(stats)

dir = "./Output-Live/Adjustment-0"
figdir = "./Output-Live/Adjustment-0/Figures"
label = "Adjustment-0"

gm.plot_residuals(figdir, label, input_data, reference, fit_results, stats)
gm.plot_correlation_matrix(figdir, label, cyy, figname="input-correlation.png")
gm.save(dir, label, reference, fit_results, stats, keys=gm.get_long_keys(input_data))

gm.plot_fit_results(
    dir, reference, [fit_results, fit_results], [label, label]
)

