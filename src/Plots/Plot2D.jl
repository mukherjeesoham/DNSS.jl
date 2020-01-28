#--------------------------------------------------------------------
# DNSS.jl
# Soham 01-2019
# Plotting routines for 2D
#--------------------------------------------------------------------
using PyPlot
export plot, pcolormesh, contourf

function PyPlot. pcolormesh(f::Field{ProductSpace{S1, S2}}) where {S1, S2} 
    # FIXME: Check this against known functions. 
    u  = Field(f.space.S1, u->u)
    v  = Field(f.space.S2, v->v)
    wu = [integral(f.space.S1, i) for i in 1:size(f.space.S1)]
    wv = [integral(f.space.S2, j) for j in 1:size(f.space.S2)]
    cu = u.value .- (wu/2) 
    cv = v.value .- (wv/2) 
    append!(cu, u.value[end] + (wu[end]/2))
    append!(cv, v.value[end] + (wv[end]/2))
    pcolormesh(cv, cu, f.value, snap=true)
    xlabel("v")
    ylabel("u")
    colorbar()
end

function PyPlot. contourf(f::Field{ProductSpace{S1, S2}}, levels::Union{Array{Number, 1}, Int}) where {S1, S2} 
    u  = Field(f.space.S1, u->u)
    v  = Field(f.space.S2, v->v)
    contourf(v.value, u.value, f.value, levels)
    contourlines = contour(v.value, u.value, f.value, levels, colors="k")
    clabel(contourlines, inline=1, fontsize=5, colors="k")
    xlabel("v")
    ylabel("u")
    colorbar()
end

