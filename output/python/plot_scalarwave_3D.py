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
    deltapsi_pconv = np.load("../data/scalar_wave_spherical_symmetry/scalar_waave_spherical_symmetry_off_axis_pconv")
    deltapsi_hconv = np.load("../data/scalar_wave_spherical_symmetry/scalar_waave_spherical_symmetry_off_axis_hconv")

    fig, (ax1, ax2) = plt.subplots(1, 2, sharey=False, sharex=False)
    ax1.semilogy(deltapsi_pconv["n"], deltapsi_pconv["e"], "-o", label= r"$k = 4$")

    ax2.plot(np.log2(deltapsi_hconv["n"][1:]), np.true_divide(deltapsi_hconv["e"][:-1], deltapsi_hconv["e"][1:]), "-o", label= r"$k = 4$")

    ax1.set_xlabel(r"$p$")
    ax1.set_ylabel(r"$\|\mathcal{E}_p\|_2$")
    ax1.legend(frameon=False)

    ax2.set_xlabel(r"$h$")
    ax2.set_ylabel(r"$\|\mathcal{E}_{h}\|_2 ~/~ \|\mathcal{E}_{h+1}\|_2$")
    ax2.legend(frameon=False)

    ax1.tick_params(axis='both', which='major', size=10)
    ax2.tick_params(axis='both', which='major', size=10)

    fig.subplots_adjust(hspace=0.1)
    plt.tight_layout()
    fig.savefig("scalar_wave-3D-off-axis-convergence.pdf")
    return 0

# Plot the solution and the associated error
def plot_solution():
    fpsi = glob.glob("../data/scalarwave_3D/k_04_np_18_nh_1/scalarwave_3D_psi*")
    fdeltapsi = glob.glob("../data/scalarwave_3D/k_04_np_18_nh_1/scalarwave_3D_deltapsi*")    

    psi_maximum = np.amax(list(map(lambda x: np.amax(np.load(x)["w"]), fpsi)))
    psi_minimum = np.amax(list(map(lambda x: np.amin(np.load(x)["w"]), fpsi)))

    deltapsi_maximum = np.amax(list(map(lambda x: np.amax(np.load(x)["w"]), fdeltapsi)))
    deltapsi_minimum = np.amax(list(map(lambda x: np.amin(np.load(x)["w"]), fdeltapsi)))
 
    psi_levels = np.linspace(psi_minimum - 1.2, psi_maximum + 0.2, 40)
    deltapsi_levels = np.linspace(deltapsi_minimum - 0.1, deltapsi_maximum + 0.1, 40)
    
    fig, (ax1, ax2) = plt.subplots(1, 2, sharey=False, sharex=False)

    for fname in fdeltapsi:
        deltapsi = np.load(fname)
        A1 = ax1.contourf(deltapsi["v"], deltapsi["u"], deltapsi["w"], 
                vmin=np.amin(deltapsi_levels) , vmax=np.amax(deltapsi_levels), levels=deltapsi_levels)

    for fname in fpsi:
        psi = np.load(fname)
        A2 = ax2.semilogy(np.amax(np.abs(psi["wlm"]), axis=1), "m-o")

    # ax1.tick_params(axis='both', which='major', size=10)
    ax1.tick_params(axis='both', which='major', size=10)
    # ax1.xaxis.set_major_locator(MaxNLocator(integer=True))
    # ax1.yaxis.set_major_locator(MaxNLocator(integer=True))

    ax2.set_xlabel(r"$p_u$")
    ax2.set_ylabel(r"$\max_{p_v}|\Psi_{{l_u}, {l_v}}|$")

    ax1.set_xlabel(r"$u$")
    ax1.set_ylabel(r"$v$")

    fig.colorbar(A1, ax=ax1)
    # fig.colorbar(A2, ax=ax2)
    fig.subplots_adjust(hspace=0.1)
    plt.tight_layout()
    fig.savefig("scalarwave_3D_solution.pdf")


    return 0

# plot_convergence()
plot_solution()
