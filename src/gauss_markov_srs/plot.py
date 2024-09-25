import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Rectangle
from matplotlib.ticker import FormatStrFormatter, MultipleLocator

from gauss_markov_srs.cross_index import get_cross_index
from gauss_markov_srs.gauss_markov import FitResult, FitStats

# plot params
params = {
    "figure.figsize": (8.0, 6.0),
    "axes.labelsize": "small",
    "axes.titlesize": "medium",
    "xtick.labelsize": "small",
    "ytick.labelsize": "small",
    "legend.fontsize": "small",
    #'figure.subplot.left'    : 0.22, # the left side of the subplots of the figure
    "figure.subplot.right": 0.96,  # the right side of the subplots of the figure
    "figure.autolayout": True,
}


def plot_correlation_matrix(basename, cyy):
    plt.figure()
    plt.title(f"{basename} Correlation matrix")
    im = plt.imshow(cyy - np.eye(cyy.shape[0]), cmap="seismic", vmin=-0.8, vmax=0.8)
    plt.colorbar(im)

    plt.savefig(basename + "-corr.pdf")


def plot_residuals(basename, data, reference, fit_result: FitResult, fit_stats: FitStats, mask=None):
    ref_str_to_i, i_to_ref_str = get_cross_index(reference["Atom"])

    mpl.rcParams.update(mpl.rcParamsDefault)
    plt.style.use(params)

    if mask is None:
        fit_mask = np.ones(len(data)).astype(bool)
    else:
        fit_mask = np.asanyarray(mask).astype(bool)

    # to convert data for plotting, I need a conversion factor that are the adjusted constants + 0 for Cs
    k = np.append([0.0], fit_result.q)
    ku = np.append([0.0], fit_result.qu)

    for entry in reference[1:]:
        ref_id, ref_atom, ref_value = entry

        mask1 = data["Atom1"] == ref_atom
        mask2 = data["Atom2"] == ref_atom

        plot_mask = mask1 | mask2
        N = sum(plot_mask)

        width, height = 8, max(6, N * 0.16 + 0.5)

        plt.figure(figsize=(width, height))
        plt.title(ref_atom)

        ref_unc = ku[ref_str_to_i(ref_atom)]
        plt.axvline(
            0,
            color="gray",
            linewidth=1,
        )
        plt.axvspan(-ref_unc, ref_unc, alpha=0.25, color="gray")

        X = fit_stats.residuals[plot_mask]
        U = fit_stats.yu[plot_mask]
        ref1 = data["Atom1"][plot_mask]
        ref2 = data["Atom2"][plot_mask]
        input_used = fit_mask[plot_mask]

        # some data may require a change of sign
        mask3 = ref2 == ref_atom
        X[mask3] = -X[mask3]
        ref = ref2
        ref[mask3] = ref1[mask3]

        # y labels
        Y = np.arange(N)

        labels = data["Ref"][plot_mask].astype(str)
        markers = ["o", "d", "s", "v", "p", "^", "X", "*", "<", ">", "P", "|", "2", "h", "."]

        mark = 0
        for i, entry2 in enumerate(reference):
            ref_id2, ref_atom2, ref_value2 = entry2

            plot_mask2 = ref == ref_atom2

            if plot_mask2.any():
                plt.errorbar(
                    x=X[plot_mask2],
                    y=Y[plot_mask2],
                    xerr=(U[plot_mask2] ** 2 + ku[i] ** 2) ** 0.5,
                    label="vs {}".format(ref_atom2),
                    fmt=markers[mark],
                    color=f"C{mark}",
                )

                not_used = ~input_used[plot_mask2]
                if not_used.any():
                    plt.errorbar(
                        x=X[plot_mask2][not_used],
                        y=Y[plot_mask2][not_used],
                        fmt=markers[mark],
                        color=f"C{mark}",
                        markerfacecolor="white",
                    )

                mark += 1

        plt.legend(loc=0)

        plt.ylim(-0.5, N - 0.5)
        plt.grid(which="both")

        # adjust too big xlim
        xl = plt.xlim()
        plt.xlim(np.clip(xl, -40 * ref_unc, 40 * ref_unc))

        plt.yticks(np.arange(len(labels)), list(labels))
        plt.xlabel("$y$")

        plt.tight_layout()
        figname = basename + "-fig{}-{}.png".format(ref_id, ref_atom)
        plt.savefig(figname)
        plt.close()


# 	# to convert data for plotting, I need a conversion factor that are the adjusted constants + 0 for Cs
# 	k = np.append([0.], q)
# 	ku = np.append([0.], qu)


# 		# y labels
# 		Y = np.arange(N)
# 		#labels = ["{}({})".format(a,b) if b  else "{}".format(a) for a, b in zip(data['Id'][mask], data['Sup'][mask]) ]
# 		labels = data['Id'][mask]
# 		markers = ['o', 'd', 's', 'v', 'p', '^', 'X', '*',  '<', '>', 'P', '|', '2', 'h']

# 		# I have to convert the data for plotting
# 		# note the change of sign
# 		k = 0
# 		for i, atom2 in enumerate(Atoms):
# 			mask = ref == atom2

# 			if mask.any():
# #				if i == 0:
# #					fmt = 'o'
# #				else:
# #					fmt = 'd'

# 				plt.errorbar(x = X[mask], y = Y[mask], xerr = (U[mask]**2 + ku[i]**2)**0.5, label='vs {}'.format(atom2), fmt=markers[k])
# 				k+=1
