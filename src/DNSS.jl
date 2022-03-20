#--------------------------------------------------------------------
# DNSS.jl
# Soham 03-2022
#--------------------------------------------------------------------

module DNSS

    using PyPlot, NLsolve, Random, LaTeXStrings

    export distribute
    export Manifold, Space, ProductSpace, SingleSpaces,
           Field, Operator, U, V, Parameters
    export plot, contourf, plotpconv, plothconv
    export extract, range, Field, stagger
    export reshapeFromTuple, reshapeToTuple, enforcebc!
    export setup, distribute
    export Cardinal, ChebyshevGL, Chebyshev

    include("AbstractTypes.jl")
    include("Spectral/ChebyshevGL.jl")
    include("Spectral/1Dspace.jl")
    include("Spectral/2Dspace.jl")
    include("Spectral/AnySpace.jl")
    include("Spectral/BasisTransform1D.jl")
    include("Spectral/BasisTransform2D.jl")
    include("ArraySpace.jl")
    include("BoundaryUtils.jl")
    include("SolverUtils.jl")
    include("Distribute.jl")
    include("Plots.jl")
end


