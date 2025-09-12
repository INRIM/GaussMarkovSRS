import os
from importlib import reload

import numpy as np
from scipy.optimize import brentq

import gauss_markov_srs as gm

gm = reload(gm)

reference = gm.load_txt("./Data/Reference/CCTF2021.txt")


input_file2021 = "./Data/Input/CCTF2021/ClockInputData2020_7.dat"
input_file = "./Data/Live/ClockInputData2025_mod_basis_n=153_TAI_05092025_AdHocCorrV1.dat"
base_dir = "./Output-Live/"


def process(input_data, y, yu, M, cyy, reference, label, mask=None, description=""):
    dir = os.path.join(base_dir, label)
    fig_dir = os.path.join(dir, "Figures")

    Cyy = gm.corr2cov(cyy, yu)
    fit_results = gm.gauss_markov_fit(M, y, Cyy, mask=mask)
    stats = gm.fit_stats(M, y, Cyy, fit_results, mask=mask)

    print("====\n")
    print(label + ": " + description + "\n")
    gm.pretty_print_fit_result(fit_results, reference)
    gm.pretty_print_stats(stats)
    gm.pretty_print_outliers(input_data, stats, mask=mask)
    print("====\n")

    gm.plot_residuals(fig_dir, label, input_data, reference, fit_results, stats, mask=mask)
    gm.plot_correlation_matrix(fig_dir, label, cyy, figname="input-correlation.png")
    gm.save(dir, label, reference, fit_results, stats, keys=gm.get_long_keys(input_data))

    gm.plot_ratio_residuals(
        fig_dir,
        ["27Al+", "171Yb+E3", "171Yb", "87Sr"],
        label,
        input_data,
        reference,
        fit_results,
        stats,
        mask=mask,
    )

    return {
        "label": label,
        "description": description,
        "results": fit_results,
        "stats": stats,
    }


def contains_any_integer(strings, integers):
    int_strs = [str(i) for i in integers]
    return np.array([any(s_int in s for s_int in int_strs) for s in strings])


adjustments = []

input_data2021, input_corr2021 = gm.parse_wgfs_file(input_file2021, reference)
input_data, input_corr = gm.parse_wgfs_file(input_file, reference)

rocit_ptb_sr3 = contains_any_integer(input_data["Id"], [117, 127, 142, 147, 151, 153])
rocit_syrte_sr2 = contains_any_integer(input_data["Id"], [118, 143, 146, 148, 149, 150, 152])
rocit = input_data["Ref"] == "Lindvall2025"
bacon = input_data["Ref"] == "Beloy2021"
ptb = input_data["Ref"] == "Dorscher2020"
riken = input_data["Ref"] == "Nemitz2016"
absolute = input_data["Atom2"] == "133Cs"

# 2021
y, yu = gm.get_y(input_data2021, reference)
M = gm.get_model_matrix(input_data2021, reference)
cyy = gm.get_corr_matrix(input_data2021, input_corr2021)
wgfs2021 = process(input_data2021, y, yu, M, cyy, reference, label="WGFS2021", description="Unexpanded")
adjustments.append(wgfs2021)

# BASE
y, yu = gm.get_y(input_data, reference)
M = gm.get_model_matrix(input_data, reference)
cyy = gm.get_corr_matrix(input_data, input_corr)
adj1 = process(input_data, y, yu, M, cyy, reference, label="Adj 1", description="All inputs")
adjustments.append(adj1)

# NO CORR
adj1nc = process(input_data, y, yu, M, np.eye(len(y)), reference, label="Adj 1nc", description="No correlations")
adjustments.append(adj1nc)

# EXCLUDING DATA
mask2 = ~rocit_ptb_sr3
adj2 = process(input_data, y, yu, M, cyy, reference, label="Adj 2", description="A1 - ROCIT PTB-Sr3", mask=mask2)
adjustments.append(adj2)


mask3 = ~(rocit_ptb_sr3 | rocit_syrte_sr2)
adj3 = process(
    input_data,
    y,
    yu,
    M,
    cyy,
    reference,
    label="Adj 3",
    description="A2 - ROCIT SYRTE-Sr2",
    mask=mask3,
)
adjustments.append(adj3)


# EXPANDING UNC DATA
mask = rocit | bacon | ptb | riken
new_yu = yu.copy()
new_yu[mask] *= 2
adj4 = process(
    input_data, y, new_yu, M, cyy, reference, label="Adj 4", description="A1 + ROCIT/BACON/PTB/RIKEN unc. x2"
)
adjustments.append(adj4)


# Only significant results
mask5 = adj3["stats"].self_sensitivity > 0.01
adj5 = process(
    input_data,
    y,
    yu,
    M,
    cyy,
    reference,
    label="Adj 5",
    description="A3 + self-sensitivity>0.01",
    mask=mask5,
)
adjustments.append(adj5)

# RANDOM EFFECT MODEL
Cyy = gm.corr2cov(cyy, yu)


def random_effect_model(b, mask):
    Cb = Cyy + b**2 * 1e-32 * np.eye(len(y))
    fit_results_temp = gm.gauss_markov_fit(M, y, Cb, mask=mask)
    stats_temp = gm.fit_stats(M, y, Cb, fit_results_temp, mask=mask)
    return stats_temp.chi2 / stats_temp.n_dof - 1.0


b = brentq(random_effect_model, 0, 5, args=None)
b *= 1e-16
Cb = Cyy + b**2 * np.eye(len(y))
cb = gm.cov2corr(Cb)
adj6 = process(
    input_data,
    y,
    gm.cov2unc(Cb),
    M,
    cb,
    reference,
    label="Adj 6",
    description=f"A1 + random effect model u = {b:.2e}",
)
adjustments.append(adj6)


b = brentq(random_effect_model, 0, 5, args=mask3)
b *= 1e-16
Cb = Cyy + b**2 * np.eye(len(y))
cb = gm.cov2corr(Cb)
adj7 = process(
    input_data,
    y,
    gm.cov2unc(Cb),
    M,
    cb,
    reference,
    label="Adj 7",
    description=f"A3 + random effect model u = {b:.2e}",
    mask=mask3,
)
adjustments.append(adj7)

# EXCLUDING DATA
mask8 = absolute
adj8 = process(input_data, y, yu, M, cyy, reference, label="Adj 8", description="No optical ratios", mask=mask8)
adjustments.append(adj8)


# SUMMARIES

# adjustments = [wgfs2021, adj1, adj1nc, adj2, adj3, adj4, adj5, adj6]

results = [_["results"] for _ in adjustments]
labels = [_["label"] + ": " + _["description"] for _ in adjustments]


gm.plot_fit_results(base_dir, reference, results, labels, figname="fit-result.png", ylim=(-60, 20))
gm.plot_fit_results(base_dir, reference, results, labels, figname="fit-result-zoom.png", ylim=(-5, 5))

print("=====\n")
for ad in adjustments:
    label = ad["label"]
    desc = ad["description"]
    stats = ad["stats"]

    string = (
        f"Nmeas = {stats.n_included},\tB = {stats.birge:.3f},\t2-sigma outliers = {np.sum(np.abs(stats.relative_residuals)>2)}"
        + "\t"
        + label
        + ": "
        + desc
    )

    print(string)


print("=====\n")
for ad in adjustments:
    label = ad["label"]
    desc = ad["description"]
    stats = ad["stats"]
