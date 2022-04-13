#--------------------------------------------------------------------
# DNSS.jl 
# Soham 03-2022
# Plotting routines for 1D, 2D and Arrays. 
# TODO: Plotting still royally sucks in Julia. I wouldn't 
# waste time with this. 
#--------------------------------------------------------------------

using PyPlot, LaTeXStrings, LinearAlgebra


"""
    Plot fields in 1D space
"""
function PyPlot. plot(u::Field{S}, path) where {S}
    rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
    rcParams["font.size"] = 8
    rcParams["text.usetex"] = true
    rcParams["axes.titlesize"] = 8.0
    rcParams["axes.labelsize"] = 10.0
    rcParams["xtick.labelsize"] = 8.0
    rcParams["ytick.labelsize"] = 8.0
    rcParams["legend.fontsize"] =  8.0
    rcParams["figure.figsize"] = (10, 10)
    rcParams["figure.dpi"] = 300
    rcParams["savefig.dpi"] = 300
    w, h = plt[:figaspect](10)
    figure(figsize=(w,h))

    x = Field(u.space, x->x)
    plot(x.value, u.value, "-") 

    # xlabel(L"$v$")
    # ylabel(L"$a(u_0, v)$")
    tight_layout()
    savefig("$path")
    close()
end

"""
    Plot fields in 2D space
"""
function PyPlot. contourf(w::Field{ProductSpace{S1, S2}}, levels; wmin=nothing, wmax=nothing) where {S1, S2} 
    u  = Field(w.space.S1, u->u)
    v  = Field(w.space.S2, v->v)
    contourf(v.value, u.value, w.value, levels; vmin=wmin, vmax=wmax, cmap="Blues")
end

"""
    Plot fields in 2D array space
"""
function PyPlot. contourf(AoF::T, levels::Int, path::String) where {T<:Array{Field, 2}}

    rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
    rcParams["font.size"] = 8
    rcParams["text.usetex"] = true
    rcParams["axes.titlesize"] = 8.0
    rcParams["axes.labelsize"] = 10.0
    rcParams["xtick.labelsize"] = 8.0
    rcParams["ytick.labelsize"] = 8.0
    rcParams["legend.fontsize"] =  8.0
    rcParams["figure.figsize"] = (10, 10)
    rcParams["figure.dpi"] = 300
    rcParams["savefig.dpi"] = 300
    # w, h = plt[:figaspect](10)
    # figure(figsize=(w,h))

    globalmin = minimum(AoF)
    globalmax = maximum(AoF)
    globallevels = collect(range(globalmin, stop=globalmax, length=levels))
    for index in CartesianIndices(AoF)
        if eltype(AoF[index].value) <: Number
            contourf(AoF[index], globallevels, wmin=globalmin, wmax=globalmax)
        end
    end
    colorbar()
    xlabel(L"$v$")
    ylabel(L"$u$")
    tight_layout()
    savefig("$path")
    close("all")
end

"""
Plot p convergence
"""
function plotpconv(n_::Vector, l_::Vector, path::String)
    rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
    rcParams["font.size"] = 8
    rcParams["text.usetex"] = true
    rcParams["axes.titlesize"] = 8.0
    rcParams["axes.labelsize"] = 10.0
    rcParams["xtick.labelsize"] = 8.0
    rcParams["ytick.labelsize"] = 8.0
    rcParams["legend.fontsize"] =  8.0
    rcParams["figure.figsize"] = (10, 10)
    rcParams["figure.dpi"] = 300
    rcParams["savefig.dpi"] = 300

    w, h = plt[:figaspect](10)
    figure(figsize=(w,h))

    semilogy(n_, l_, "o--", linewidth=0.5,  markersize=1.0)
    xlabel(L"$N_p$")
    xticks(collect(12:2:24))
    ylabel(L"$\|\mathcal{C}\|_2$")
    tight_layout()
    savefig("$path")
    close()
end

"""
Plot h convergence
"""
function plothconv(n_::Vector, l_::Vector, path::String)
    rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
    rcParams["font.size"] = 8
    rcParams["text.usetex"] = true
    rcParams["axes.titlesize"] = 8.0
    rcParams["axes.labelsize"] = 10.0
    rcParams["xtick.labelsize"] = 8.0
    rcParams["ytick.labelsize"] = 8.0
    rcParams["legend.fontsize"] =  8.0
    rcParams["figure.figsize"] = (10, 10)
    rcParams["figure.dpi"] = 300
    rcParams["savefig.dpi"] = 300
    ll_ = similar(l_) 
    w, h = plt[:figaspect](10)
    figure(figsize=(w,h))

    for index in 2:length(l_) 
        ll_[index] = l_[index - 1] / l_[index]
    end

    plot(log.(2, n_[2:end]), ll_[2:end], "o--", linewidth=0.5, markersize=1.0)
    xlabel(L"$\log_2(N_h)$")
    ylabel(L"$\Delta_{N_h} / \Delta_{N_{h+1}} $")
    tight_layout()
    savefig("$path")
    close()
end

function plotmodes(u::Field{ProductSpace{ChebyshevGL{Tag1, N1, T},
                                         ChebyshevGL{Tag2, N2, T}}}, path::String) where {Tag1, Tag2, N1, N2, T}
    rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
    rcParams["font.size"] = 8
    rcParams["text.usetex"] = true
    rcParams["axes.titlesize"] = 8.0
    rcParams["axes.labelsize"] = 10.0
    rcParams["xtick.labelsize"] = 8.0
    rcParams["ytick.labelsize"] = 8.0
    rcParams["legend.fontsize"] =  8.0
    rcParams["figure.figsize"] = (10, 10)
    rcParams["figure.dpi"] = 300
    rcParams["savefig.dpi"] = 300
    w, h = plt[:figaspect](10)
    figure(figsize=(w,h))

    # Compute the maximum along the secondary diagonals
    ull = basistransform(u) 
    B   = rotr90(ull.value,3)
    (lu, lv) = size(B)
    @assert lu == lv
    ul = [maximum(abs.(B[diagind(B, ind)])) for ind in -lu+1:lv-1] 

    semilogy(ul, "o--", linewidth=0.5, markersize=1.0)
    xlabel(L"$l_u + l_v$")
    ylabel(L"$\mathrm{max}(|c_{uv}|)$")
    tight_layout()
    savefig("$path")
    close()
end
