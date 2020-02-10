#--------------------------------------------------------------------
# DNSS.jl
# Soham 08-2019
# Simulate Minkowski spacetime 
#--------------------------------------------------------------------

using PyPlot
parameters = stagger(Parameters((12, 12), (1, 1), (-2.0, 2.0), (-2.0, 2.0)), 1e-2)
AoT = distribute(parameters, excision, computeUboundary, computeVboundary)
AoF = extract(AoT, 2)

contourf(AoF, 10)
plotaxis(parameters)
savefig("MinkowskiEvolveR.pdf")

L2C = maximum([maximum(L2.(constraints(T...))) for T in AoT])
@show L2C 
