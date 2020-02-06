#--------------------------------------------------------------------
# DNSS.jl
# Soham 08-2019
# Simulate Minkowski spacetime 
#--------------------------------------------------------------------

MPLBACKEND = "WXAgg"
using PyPlot

parameters = Parameters((12, 12), (10, 10), (0.0, 4.0), (0.0, 4.0))
tree = distribute(parameters, excision, computeUboundary, computeVboundary)

AoF = extract(tree, 2)

globalmin = minimum(AoF)
globalmax = maximum(AoF)
levels = collect(range(globalmin, stop=globalmax, length=10))

for index in CartesianIndices(AoF)
    if eltype(AoF[index].value) <: Number
        contourf(AoF[index], levels)
    end
    xlim(parameters.vbounds)
    ylim(parameters.ubounds)
    plot(parameters.vbounds, parameters.ubounds)
    colourbar()
end

show()

