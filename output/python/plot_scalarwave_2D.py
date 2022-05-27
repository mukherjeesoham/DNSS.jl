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
    deltapsi_0p1_pconv = np.load("../data/scalarwave_2D/sigma_0p1/minkowski_2D_sigma_0p1_deltapsi_pconv")
    deltapsi_0p3_pconv = np.load("../data/scalarwave_2D/sigma_0p3/minkowski_2D_sigma_0p3_deltapsi_pconv")
    deltapsi_0p5_pconv = np.load("../data/scalarwave_2D/sigma_0p5/minkowski_2D_sigma_0p5_deltapsi_pconv")
    deltapsi_0p7_pconv = np.load("../data/scalarwave_2D/sigma_0p7/minkowski_2D_sigma_0p7_deltapsi_pconv")
    deltapsi_0p9_pconv = np.load("../data/scalarwave_2D/sigma_0p9/minkowski_2D_sigma_0p9_deltapsi_pconv")
    
    deltapsi_0p1_hconv = np.load("../data/scalarwave_2D/sigma_0p1/minkowski_2D_sigma_0p1_deltapsi_hconv")
    deltapsi_0p3_hconv = np.load("../data/scalarwave_2D/sigma_0p3/minkowski_2D_sigma_0p3_deltapsi_hconv")
    deltapsi_0p5_hconv = np.load("../data/scalarwave_2D/sigma_0p5/minkowski_2D_sigma_0p5_deltapsi_hconv")
    deltapsi_0p7_hconv = np.load("../data/scalarwave_2D/sigma_0p7/minkowski_2D_sigma_0p7_deltapsi_hconv")
    deltapsi_0p9_hconv = np.load("../data/scalarwave_2D/sigma_0p9/minkowski_2D_sigma_0p9_deltapsi_hconv")

    fig, (ax1, ax2) = plt.subplots(1, 2, sharey=False, sharex=False)
    # ax1.semilogy(deltapsi_0p1_pconv["n"], deltapsi_0p1_pconv["e"], "-o", label= r"$\sigma = 0.1$")
    ax1.semilogy(deltapsi_0p3_pconv["n"], deltapsi_0p3_pconv["e"], "-o", label= r"$\sigma = 0.3$")
    ax1.semilogy(deltapsi_0p5_pconv["n"], deltapsi_0p5_pconv["e"], "-o", label= r"$\sigma = 0.5$")
    ax1.semilogy(deltapsi_0p7_pconv["n"], deltapsi_0p7_pconv["e"], "-o", label= r"$\sigma = 0.7$")
    # ax1.semilogy(deltapsi_0p9_pconv["n"], deltapsi_0p9_pconv["e"], "-o", label= r"$\sigma = 0.9$")

    # ax2.plot(np.log2(deltapsi_0p1_hconv["n"][1:]), np.true_divide(deltapsi_0p1_hconv["e"][:-1], deltapsi_0p1_hconv["e"][1:]), "-o", label= r"$\sigma = 0.1$")
    ax2.plot(np.log2(deltapsi_0p3_hconv["n"][2:]), np.true_divide(deltapsi_0p3_hconv["e"][:-1], deltapsi_0p3_hconv["e"][1:])[1:], "-o", label= r"$\sigma = 0.3$")
    ax2.plot(np.log2(deltapsi_0p5_hconv["n"][2:]), np.true_divide(deltapsi_0p5_hconv["e"][:-1], deltapsi_0p5_hconv["e"][1:])[1:], "-o", label= r"$\sigma = 0.5$")
    ax2.plot(np.log2(deltapsi_0p7_hconv["n"][2:]), np.true_divide(deltapsi_0p7_hconv["e"][:-1], deltapsi_0p7_hconv["e"][1:])[1:], "-o", label= r"$\sigma = 0.7$")
    # ax2.plot(np.log2(deltapsi_0p9_hconv["n"][1:]), np.true_divide(deltapsi_0p9_hconv["e"][:-1], deltapsi_0p9_hconv["e"][1:]), "-o", label= r"$\sigma = 0.9$")

    ax1.set_xlabel(r"$p$")
    ax1.set_ylabel(r"$\|\mathcal{E}_p\|_2$")
    ax1.legend(frameon=False)

    ax2.set_xlabel(r"$h$")
    ax2.set_ylabel(r"$\|\mathcal{E}_{h}\|_2 ~/~ \|\mathcal{E}_{h+1}\|_2$")
    ax2.legend(frameon=False)


    ax2.set_xlim(1.5, 7.5)
    ax1.tick_params(axis='both', which='major', size=10)
    ax2.tick_params(axis='both', which='major', size=10)

    fig.subplots_adjust(hspace=0.1)
    plt.tight_layout()
    fig.savefig("scalarwave_2D_convergence.pdf")
    return 0

# Plot the solution and the associated error
def plot_solution():
    fpsi = glob.glob("../data/minkowski_2D/sigma_0p5_np_30_nh_4/minkowski_2D_sigma_0p5_psi_*")
    fpsi_source = glob.glob("../data/minkowski_2D/sigma_0p5_np_30_nh_4_source/minkowski_2D_sigma_0p5_psi_*")
    
    psi_maximum = np.amax(list(map(lambda x: np.amax(np.load(x)["w"]), fpsi)))
    psi_minimum = np.amax(list(map(lambda x: np.amin(np.load(x)["w"]), fpsi)))

    psi_source_maximum = np.amax(list(map(lambda x: np.amax(np.load(x)["w"]), fpsi_source)))
    psi_source_minimum = np.amax(list(map(lambda x: np.amin(np.load(x)["w"]), fpsi_source)))
 
    psi_levels = np.linspace(psi_minimum - 1.2, psi_maximum + 0.2, 40)
    psi_source_levels = np.linspace(psi_source_minimum - 0.01, psi_source_maximum + 0.01, 40)
    
    fig, ax1 = plt.subplots(1, 1, sharey=False, sharex=True)

    for fname in fpsi:
        psi = np.load(fname)
        A1 = ax1.contourf(psi["v"], psi["u"], psi["w"], 
                vmax=np.amax(psi_levels), vmin=np.amin(psi_levels), levels=psi_levels)

    # for fname in fpsi_source:
        # print(fname)
        # psi_source = np.load(fname)
        # # A2 = ax2.contourf(psi_source["v"], psi_source["u"], psi_source["w"], 
                # # vmin=np.amin(psi_source_levels) , vmax=np.amax(psi_source_levels), levels=psi_source_levels)
        # A2 = ax2.contourf(psi_source["v"], psi_source["u"], psi_source["w"], levels=psi_source_levels) 
                # # vmin=np.amin(psi_source_levels) , vmax=np.amax(psi_source_levels), levels=psi_source_levels)

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
# plot_solution()
