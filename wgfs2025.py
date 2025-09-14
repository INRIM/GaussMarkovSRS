import os
from importlib import reload

import numpy as np
from scipy.optimize import brentq

import gauss_markov_srs as gm

gm = reload(gm)

reference = gm.load_txt("./Data/Reference/CCTF2021.txt")


input_file2021 = "./Data/Input/CCTF2021/ClockInputData2020_7.dat"
input_file = "./Data/Live/ClockInputData2025_mod_basis_n=153_TAI_05092025_ExpandUncV1_AdHocCorrV2.dat"
base_dir = "./Output/"


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


adjustments = []

input_data2021, input_corr2021 = gm.parse_wgfs_file(input_file2021, reference)
input_data, input_corr = gm.parse_wgfs_file(input_file, reference)


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

# The input file exclude some measurements by enlarging the uncertainty by 1e4
# Here I want to receove the original unc for plotting
# as I can exclude the data with a mask
s_to_i, i_to_s = gm.get_cross_index(input_data["Id"])
excluded = s_to_i(["q51", "q117", "q127", "q142", "q147", "q151", "q153"])
yu[excluded] /= 1e4

mask1 = np.ones_like(y).astype(bool)
mask1[excluded] = False

adj1 = process(input_data, y, yu, M, cyy, reference, label="Adj 1", description="2025 Inputs", mask=mask1)
adjustments.append(adj1)


# RANDOM EFFECT MODEL
mask2 = mask1
Cyy = gm.corr2cov(cyy, yu)


def random_effect_model(b, mask):
    Cb = Cyy + b**2 * 1e-32 * np.eye(len(y))
    fit_results_temp = gm.gauss_markov_fit(M, y, Cb, mask=mask)
    stats_temp = gm.fit_stats(M, y, Cb, fit_results_temp, mask=mask)
    return stats_temp.chi2 / stats_temp.n_dof - 1.0


b = brentq(random_effect_model, 0, 5, args=mask2)
b *= 1e-16
Cb = Cyy + b**2 * np.eye(len(y))
cb = gm.cov2corr(Cb)
adj2 = process(
    input_data,
    y,
    gm.cov2unc(Cb),
    M,
    cb,
    reference,
    label="Adj 2",
    description=f"A1 + Random effect model u = {b:.2e}",
    mask=mask2,
)
adjustments.append(adj2)
print(f"u = {b:.9e}")

np.savetxt(os.path.join(base_dir, "A2-correlation-matrix.dat"), cb)
np.savetxt(os.path.join(base_dir, "A2-expanded-unc.dat"), gm.cov2unc(Cb))

# RANDOM EFFECT MODEL 2
mask3 = mask2

def random_effect_model2(b, mask):
    yub = (yu**2 + b**2 * 1e-32) ** 0.5
    Cb = gm.corr2cov(cyy, yub)
    fit_results_temp = gm.gauss_markov_fit(M, y, Cb, mask=mask)
    stats_temp = gm.fit_stats(M, y, Cb, fit_results_temp, mask=mask)
    return stats_temp.chi2 / stats_temp.n_dof - 1.0

b = brentq(random_effect_model2, 0, 5, args=mask3)
yub = (yu**2 + b**2 * 1e-32) ** 0.5
Cb = gm.corr2cov(cyy, yub)
adj3 = process(
    input_data,
    y,
    yub,
    M,
    cyy,
    reference,
    label="Adj 3",
    description=f"Test: A1 + random effect model (naive), u = {b:.2e}",
    mask=mask3,
)
adjustments.append(adj3)
print(f"u = {b:.9e}")

# Only significant results
mask4 = adj2["stats"].self_sensitivity > 0.01
adj4 = process(
    input_data,
    y,
    gm.cov2unc(Cb),
    M,
    cb,
    reference,
    label="Adj 4",
    description="Test: A2 + self-sensitivity>0.01",
    mask=mask4,
)
adjustments.append(adj4)


# Only absolute
absolute = input_data["Atom2"] == "133Cs"
mask5 = absolute & mask2
adj5 = process(input_data, y, yu, M, cyy, reference, label="Adj 5", description="Test: No optical ratios", mask=mask5)
adjustments.append(adj5)




# SUMMARIES
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
        f"Nmeas = {stats.n_included},\tB = {stats.birge:.3f} (p={stats.p_value:.3f}),\t2-sigma outliers = {np.sum(np.abs(stats.relative_residuals)>2)}"
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
