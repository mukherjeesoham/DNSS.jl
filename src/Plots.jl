#--------------------------------------------------------------------
# DNSS.jl 
# Soham 03-2022
# Plotting routines for 1D, 2D and Arrays. 
# TODO: Plotting still royally sucks in Julia. I wouldn't 
# waste time with this. 
#--------------------------------------------------------------------

rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
rcParams["font.size"] = 15
rcParams["text.usetex"] = true
rcParams["axes.titlesize"] = 18.0
rcParams["axes.labelsize"] = 15.0
rcParams["xtick.labelsize"] = 15.0
rcParams["ytick.labelsize"] = 15.0
rcParams["legend.fontsize"] =  15.0
rcParams["figure.figsize"] = (12, 10)
rcParams["figure.dpi"] = 300
rcParams["savefig.dpi"] = 300

"""
    Plot fields in 1D space
"""
function PyPlot. plot(u::Field{S}; plotstyle="-o", label="") where {S}
    x = Field(u.space, x->x)
    plot(x.value, u.value, plotstyle, label=label) 
end

"""
    Plot fields in 2D space
"""
function PyPlot. contourf(w::Field{ProductSpace{S1, S2}}, levels; wmin=nothing, wmax=nothing) where {S1, S2} 
    u  = Field(w.space.S1, u->u)
    v  = Field(w.space.S2, v->v)
    contourf(v.value, u.value, w.value, levels; vmin=wmin, vmax=wmax)
end

"""
    Plot fields in 2D array space
"""
function PyPlot. contourf(AoF::T, levels::Int, path::String) where {T<:Array{Field, 2}}
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
    savefig("$path")
    close()
end

"""
Plot p convergence
"""
function plotpconv(n_::Vector, l_::Vector, path::String)
    semilogy(n_, l_, "o--")
    xlabel(L"$p$")
    ylabel(L"$L_{2}(u - u_{0}) $")
    savefig("$path")
    close()
end

"""
Plot h convergence
"""
function plothconv(n_::Vector, l_::Vector, path::String)
    for index in 2:length(l_) 
        plot(n_[index], l_[index - 1] / l_[index], "o--")
    end
    xlabel(L"$2^h$")
    ylabel(L"$L_{2}(u - u_{0}) $")
    savefig("$path")
    close()
end

