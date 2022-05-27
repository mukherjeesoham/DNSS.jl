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

# Plot convergence
def plot_convergence():
    deltapsi_0p1_pconv = np.load("../data/scalarwave_SG/scalarwave_SG_k0p1_off_axis_pconv")
    deltapsi_0p3_pconv = np.load("../data/scalarwave_SG/scalarwave_SG_k0p3_off_axis_pconv")
    deltapsi_0p5_pconv = np.load("../data/scalarwave_SG/scalarwave_SG_k0p5_off_axis_pconv")
    
    deltapsi_0p1_hconv = np.load("../data/scalarwave_SG/scalarwave_SG_k0p1_off_axis_hconv")
    deltapsi_0p3_hconv = np.load("../data/scalarwave_SG/scalarwave_SG_k0p3_off_axis_hconv")
    deltapsi_0p5_hconv = np.load("../data/scalarwave_SG/scalarwave_SG_k0p5_off_axis_hconv")

    fig, (ax1, ax2) = plt.subplots(1, 2, sharey=False, sharex=False)
    ax1.semilogy(deltapsi_0p1_pconv["n"][:-2], deltapsi_0p1_pconv["e"][:-2], "-o", label= r"$k = 0.1$")
    ax1.semilogy(deltapsi_0p3_pconv["n"][:-1], deltapsi_0p3_pconv["e"][:-1], "-o", label= r"$k = 0.3$")
    ax1.semilogy(deltapsi_0p5_pconv["n"],      deltapsi_0p5_pconv["e"],      "-o", label= r"$k = 0.5$")

    ax2.plot(np.log2(deltapsi_0p1_hconv["n"][1:]), 2 * np.true_divide(deltapsi_0p1_hconv["e"][:-1], deltapsi_0p1_hconv["e"][1:])[0:], "-o", label= r"$k = 0.1$")
    ax2.plot(np.log2(deltapsi_0p3_hconv["n"][1:]), 2 * np.true_divide(deltapsi_0p3_hconv["e"][:-1], deltapsi_0p3_hconv["e"][1:])[0:], "-o", label= r"$k = 0.3$")
    ax2.plot(np.log2(deltapsi_0p5_hconv["n"][1:]), 2 * np.true_divide(deltapsi_0p5_hconv["e"][:-1], deltapsi_0p5_hconv["e"][1:])[0:], "-o", label= r"$k = 0.5$")

    ax1.set_xlabel(r"$p$")
    ax1.set_ylabel(r"$\|\mathcal{C}_p\|_2$")
    ax1.legend(frameon=False)

    ax2.set_xlabel(r"$h$")
    ax2.set_ylabel(r"$\|\mathcal{C}_{h}\|_2 ~/~ \|\mathcal{C}_{h+1}\|_2$")
    ax2.legend(frameon=False)


    ax2.set_xlim(1.5, 7.5)
    ax1.tick_params(axis='both', which='major', size=10)
    ax2.tick_params(axis='both', which='major', size=10)

    fig.subplots_adjust(hspace=0.1)
    plt.tight_layout()
    fig.savefig("scalarwave_SG_convergence.pdf")
    return 0

# Plot the solution and the associated error
def plot_solution():
    uf = glob.glob("../data/scalarwave_SG/scalarwave_SG_psi*")
    vf = glob.glob("../data/scalarwave_SG/scalarwave_SG_constraints*")
    
    umax = np.amax(list(map(lambda x: np.amax(np.load(x)["w"]), uf)))
    umin = np.amax(list(map(lambda x: np.amin(np.load(x)["w"]), uf)))

    vmax = np.amax(list(map(lambda x: np.amax(np.load(x)["w"]), vf)))
    vmin = np.amax(list(map(lambda x: np.amin(np.load(x)["w"]), vf)))
 
    ulevels = np.linspace(umin, umax, 40)
    vlevels = np.linspace(vmin, vmax, 40)
    
    fig, ax1 = plt.subplots(1, 1, sharey=False, sharex=True)

    for _u in fpsi:
        psi = np.load(fname)
        A1 = ax1.contourf(psi["v"], psi["u"], psi["w"], 
                vmax=np.amax(psi_levels), vmin=np.amin(psi_levels), levels=psi_levels)


    ax1.tick_params(axis='both', which='major', size=10)
    # ax2.tick_params(axis='both', which='major', size=10)

    ax1.set_xlabel(r"$u$")
    ax1.set_ylabel(r"$v$")

    # ax2.set_xlabel(r"$u$")
    # ax2.set_ylabel(r"$v$")

    fig.colorbar(A1, ax=ax1)
    # fig.colorbar(A2, ax=ax2)
    fig.subplots_adjust(hspace=0.1)
    plt.tight_layout()
    fig.savefig("minkowski-2D-solution.pdf")

    return 0

plot_convergence()
