#-----------------------------------------------------
# Make plots from matplotlib using data exported by
# DNSS.jl
# Soham M 05/2022
#-----------------------------------------------------

import numpy as np
import glob
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

import matplotlib
mpl.rcParams.update({
        "font.size": 34.0,
        "axes.titlesize": 34.0,
        "axes.labelsize": 34.0,
        "xtick.labelsize": 34.0,
        "ytick.labelsize": 34.0,
        "legend.fontsize": 34.0,
        "figure.figsize": (25, 10),
        "figure.dpi": 300,
        "savefig.dpi": 300,
        "text.usetex": True
})

def plot_solution():
    fguv = glob.glob("../data/minkowski/constraints/minkowski_guv*")
    fgrr = glob.glob("../data/minkowski/constraints/minkowski_grr*")

    guvmax = np.amax(list(map(lambda x: np.amax(np.load(x)["w"]), fguv)))
    guvmin = np.amin(list(map(lambda x: np.amin(np.load(x)["w"]), fguv)))
    grrmax = np.amax(list(map(lambda x: np.amax(np.load(x)["w"]), fgrr)))
    grrmin = np.amin(list(map(lambda x: np.amin(np.load(x)["w"]), fgrr)))

    guvlevels = np.linspace(guvmin, guvmax, 40)
    grrlevels = np.linspace(grrmin, grrmax, 40)

    fig, (ax1, ax2) = plt.subplots(1, 2, sharey=False, sharex=True)

    for (_guv, _grr) in zip(fguv, fgrr):
        guv = np.load(_guv)
        grr = np.load(_grr)

        A1 = ax1.contourf(guv["v"], guv["u"], guv["w"], vmax=np.amax(guvlevels), vmin=np.amin(guvlevels), levels=guvlevels)
        A2 = ax2.contourf(grr["v"], grr["u"], grr["w"], vmax=np.amax(grrlevels), vmin=np.amin(grrlevels), levels=grrlevels)


    ax1.tick_params(axis='both', which='major', size=10)
    ax1.set_xlabel(r"$v$")
    ax1.set_ylabel(r"$u$")
    fig.colorbar(A1, ax=ax1)
    fig.colorbar(A2, ax=ax2)
    plt.tight_layout()
    fig.savefig("minkowski_constraints.pdf")

    return 0

plot_solution()
