#--------------------------------------------------------------------
# DNSS.jl 
# Soham 03-2022
# Plotting routines for 1D, 2D and Arrays. 
# TODO: Plotting still royally sucks in Julia. I wouldn't 
# waste time with this. 
#--------------------------------------------------------------------

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
