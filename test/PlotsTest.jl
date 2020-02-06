#--------------------------------------------------------------------
# DNSS.jl
# Soham 08-2019
# Experiment with staggered grids
#--------------------------------------------------------------------
using PyPlot

parameters = Parameters((8, 8), (34, 34), (-2.0, 2.0), (-2.0, 2.0))

# Minkowski spacetime
newparameters = stagger(parameters, 1e-6)
AoS = setup(newparameters)
AoF = Field(AoS, (u,v)->(v-u))
contourf(AoF, 10)
plotaxis(newparameters)
axis("square")
savefig("Minkowski.pdf")
close()


# Schwarzschild spacetime
AoS = setup(parameters)

function excision(PS::ProductSpace{S1, S2})::Bool where {S1, S2}
    r = Field(PS, (u,v)->v*u)
    v = Field(PS, (u,v)->v)
    return any(value(r) .>= 1) || any(value(v) .< eps(eltype(value(r)))) 
end

M = 1.0
AoF = Array{Field, 2}(undef, size(AoS))
for index in CartesianIndices(AoS)
    if excision(AoS[index]) == true
        AoF[index] = Field(AoS[index], Array{Missing,2}(missing, size(AoS[index])))
    else
        r = Field(AoS[index], (u,v)->v*u)
        AoF[index] = Field(AoS[index], (u,v)->find_r_of_UV(u,v,M))
    end
end

function plotsingularity(params::Parameters)
    v = collect(range(params.vbounds[1], stop=params.vbounds[2], length=10*params.npoints[2]))
    plot(v, 1 ./v, "k--", linewidth=0.8)
    xlim(params.vbounds)
    ylim(params.ubounds)
end

contourf(AoF, 10)
plotsingularity(parameters)
plotaxis(newparameters)
savefig("Schwarzschild.pdf")
