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
import warnings
warnings.filterwarnings("ignore")

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

def plot_convergence():
    deltapsi_0p1_pconv = np.load("../data/schwarzschild/schwarzschild_M_0p1_pconv")
    deltapsi_0p3_pconv = np.load("../data/schwarzschild/schwarzschild_M_1_pconv")
    deltapsi_0p5_pconv = np.load("../data/schwarzschild/schwarzschild_M_8_pconv")
    
    deltapsi_0p1_hconv = np.load("../data/schwarzschild/schwarzschild_M_0p1_hconv")
    deltapsi_0p3_hconv = np.load("../data/schwarzschild/schwarzschild_M_1_hconv")
    deltapsi_0p5_hconv = np.load("../data/schwarzschild/schwarzschild_M_8_hconv")

    fig, (ax1, ax2) = plt.subplots(1, 2, sharey=False, sharex=False)
    ax1.semilogy(deltapsi_0p1_pconv["n"], deltapsi_0p1_pconv["e"], "-o", label= r"$M = 0.125$")
    ax1.semilogy(deltapsi_0p3_pconv["n"], deltapsi_0p3_pconv["e"], "-o", label= r"$M = 1.0$")
    ax1.semilogy(deltapsi_0p5_pconv["n"], deltapsi_0p5_pconv["e"], "-o", label= r"$M = 8.0$")

    ax2.plot(np.log2(deltapsi_0p5_hconv["n"][1:]), np.true_divide(deltapsi_0p5_hconv["e"][:-1], deltapsi_0p5_hconv["e"][1:])[0:], "--o", markersize=20.0, label= r"$M = 8.0$")

    ax2.plot(np.log2(deltapsi_0p1_hconv["n"][1:]), np.true_divide(deltapsi_0p1_hconv["e"][:-1], deltapsi_0p1_hconv["e"][1:])[0:], "--o", markersize=10.0, label= r"$M = 0.125$")
    ax2.plot(np.log2(deltapsi_0p3_hconv["n"][1:]), np.true_divide(deltapsi_0p3_hconv["e"][:-1], deltapsi_0p3_hconv["e"][1:])[0:], "k--o", markersize=5.0, label= r"$M = 1.0$")
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
    fig.savefig("schwarzschild_convergence.pdf")
    return 0
def plot_solution():
    unames = glob.glob("../data/schwarzschild/schwarzschild_eta_*")
    umax    = np.nanmax(list(map(lambda x: np.nanmax(np.load(x)["w"]), unames)))
    umin    = np.nanmin(list(map(lambda x: np.nanmin(np.load(x)["w"]), unames)))
    ulevels = np.linspace(umin, umax, 40)
    print(ulevels)

    fig, (ax1, ax2) = plt.subplots(1, 2, sharey=False, sharex=False)

    for fname in unames:
        u = np.load(fname)
        A1 = ax1.contourf(u["v"], u["u"], np.nan_to_num(u["w"]),
                vmax=np.amax(ulevels), vmin=np.amin(ulevels), levels=ulevels)
        # ax1.plot(u["v"], 1 / u["v"], "--k")

    qnames = glob.glob("../data/schwarzschild/schwarzschild_constraints*")
    qmax    = np.nanmax(list(map(lambda x: np.nanmax(np.load(x)["w"]), qnames)))
    qmin    = np.nanmin(list(map(lambda x: np.nanmin(np.load(x)["w"]), qnames)))
    print(qmax, qmin)
    # qlevels = np.log10(np.abs(np.linspace(qmin, qmax, 40)))
    qlevels = np.arange(-10, 1, 0.6) 
    

    for fname in qnames:
        q = np.load(fname)
        lp = np.log10(np.abs(np.nan_to_num(q["w"])))
        A2 = ax2.contourf(q["v"], q["u"], lp,  vmax=np.amax(qlevels), vmin=np.amin(qlevels), levels=qlevels)
        # ax2.plot(q["v"], 1 / q["v"], "--k")

    ax1.tick_params(axis='both', which='major', size=10)
    ax1.set_xlabel(r"$v$")
    ax1.set_ylabel(r"$u$")
    ax1.set_ylim(-3, 0)
    ax1.set_xlim(2, 5)
    ax1.set_title("$r(u,v)$", pad=20)
    ax2.set_title("$\log_{10} \mathcal{|C|}$", pad=20)

    ax2.tick_params(axis='both', which='major', size=10)
    ax2.set_xlabel(r"$v$")
    ax2.set_ylabel(r"$u$")
    ax2.set_ylim(-3, 0)
    ax2.set_xlim(2, 5)

    fig.colorbar(A1, ax=ax1)
    fig.colorbar(A2, ax=ax2)
    plt.tight_layout()
    fig.savefig("schwarzschild-solution.pdf")

    return 0

plot_convergence()
