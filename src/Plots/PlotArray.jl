#--------------------------------------------------------------------
# DNSS.jl
# Soham 01-2019
# Plotting routines for Arrays
#--------------------------------------------------------------------

using PyPlot
export contourf, scatter, plotaxis

function PyPlot. contourf(AoF::T, levels::Int) where {T<:Array{Field, 2}}
    globalmin = minimum(AoF)
    globalmax = maximum(AoF)
    globallevels = collect(range(globalmin, stop=globalmax, length=levels))
    for index in CartesianIndices(AoF)
        if eltype(AoF[index].value) <: Number
            contourf(AoF[index], globallevels, wmin=globalmin, wmax=globalmax)
        end
    end
    colorbar()
    xlabel("v")
    ylabel("u")
end

function PyPlot. scatter(AoS::T) where {T<:Array{Union{ProductSpace, NTuple{N, Field}}}} where {N}
    for index in CartesianIndices(AoS)
        scatter(AoS[index])
    end
    xlabel("v")
    ylabel("u")
end

function plotaxis(params::Parameters)
    v = collect(range(params.vbounds[1], stop=params=params.vbounds[2], length=3))
    plot(v, v, "k--", linewidth=0.5)
end

