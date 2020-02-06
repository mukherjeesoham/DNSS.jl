#--------------------------------------------------------------------
# DNSS.jl
# Soham 01-2019
# Plotting routines for 2D
#--------------------------------------------------------------------

using PyPlot
export contourf, scatter

function PyPlot. contourf(w::Field{ProductSpace{S1, S2}}, levels; wmin=nothing, wmax=nothing) where {S1, S2} 
    u  = Field(w.space.S1, u->u)
    v  = Field(w.space.S2, v->v)
    contourf(v.value, u.value, w.value, levels; vmin=wmin, vmax=wmax)
end

function PyPlot. scatter(PS::ProductSpace{S1, S2}) where {S1, S2}
    u = Field(PS, (u,v)->u)
    v = Field(PS, (u,v)->v)
    for index in CartesianIndices(v.value)
        scatter(v.value[index], u.value[index], s=1)
    end
end

