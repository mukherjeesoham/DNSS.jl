#--------------------------------------------------------------------
# DNSS.jl
# Soham 03-2022
#--------------------------------------------------------------------

module DNSS

    using PyPlot, NLsolve, Random, LaTeXStrings, LinearAlgebra

    export distribute
    export Manifold, Space, ProductSpace, SingleSpaces,
           Field, Operator, U, V, Parameters
    export plot, contourf, plotpconv, plothconv
    export extract, range, Field, adjust, rmse
    export reshapeFromTuple, reshapeToTuple, enforcebc!
    export setup, distribute
    export Cardinal, ChebyshevGL, Chebyshev
    export value, space, linsolve, norm
    export basistransform, project, prolongate, restrict 

    include("AbstractTypes.jl")
    include("Spectral/ChebyshevGL.jl")
    include("Spectral/1Dspace.jl")
    include("Spectral/2Dspace.jl")
    include("Spectral/BasisTransform.jl")
    include("AnySpace.jl")
    include("ArraySpace.jl")
    include("BoundaryUtils.jl")
    include("SolverUtils.jl")
    include("Distribute.jl")
    include("Plots.jl")
end


