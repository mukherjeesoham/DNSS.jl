#--------------------------------------------------------------------
# DNSS.jl
# Soham 08-2019
# Simulate Minkowski spacetime 
#--------------------------------------------------------------------

parameters = Parameters((12, 12), (8, 8), (0.0, 4.0), (0.0, 4.0))
tree = distribute(parameters, excision, computeUboundary, computeVboundary)
