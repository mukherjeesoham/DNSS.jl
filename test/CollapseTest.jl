#--------------------------------------------------------------------
# DNSS.jl
# Soham 08-2019
# Simulate Minkowski spacetime 
#--------------------------------------------------------------------

using PyPlot
parameters = stagger(Parameters((12, 12), (4, 4), (-2.0, 2.0), (-2.0, 2.0)), 1e-2)
AoT = distribute(parameters, excision, computeUboundary, computeVboundary)

AoF = extract(AoT, 1)
contourf(AoF, 20)
plotaxis(parameters)
axis("square")
tight_layout()
savefig("CollapseA.pdf")
close()

AoF = extract(AoT, 2)
contourf(AoF, 20)
plotaxis(parameters)
axis("square")
tight_layout()
savefig("CollapseR.pdf")
close()

AoF = extract(AoT, 3)
contourf(AoF, 20)
plotaxis(parameters)
axis("square")
tight_layout()
savefig("CollapseF.pdf")
close()

L2C = maximum([maximum(L2.(constraints(T...))) for T in AoT])
@show L2C 
