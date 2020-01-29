#--------------------------------------------------------------------
# DNSS.jl
# Soham 01-2020
#--------------------------------------------------------------------

module DNSS

include("./Types/AbstractTypes.jl")
include("./Spectral/Basis/BasisTypes.jl")
include("./Spectral/Basis/ChebyshevGL.jl")
include("./Spectral/Spaces/1Dspace.jl")
include("./Spectral/Spaces/2Dspace.jl")
include("./Spectral/Spaces/AnySpace.jl")
# include("./Plots/Plot1D.jl")
# include("./Plots/Plot2D.jl")

include("./Utilities/AxiSymmetry.jl")
include("./Utilities/BoundaryUtils.jl")
include("./Utilities/NonlinearSolver.jl")

include("./Distribute/Distribute.jl")
include("./Distribute/Futures.jl")

include("./Physics/InitialData.jl")
include("./Physics/Residuals.jl")
include("./Physics/Utilities/Diagnostics.jl")
include("./Physics/Utilities/DoubleNullCoordinates.jl")

include("./Physics/Spacetimes/Minkowski.jl")

end


