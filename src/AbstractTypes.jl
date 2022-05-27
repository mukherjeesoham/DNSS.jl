#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 08-2018
# Define the Abstract Types
#--------------------------------------------------------------------
abstract type Manifold{Tag} end
abstract type Space{Tag, N} <: Manifold{Tag} end

struct U end
struct V end
struct Singularity end

"""
    Struct to hold a field
    in space S, dimension D and eltype T.
"""
struct Field{S, D, T}
    space::S
    value::AbstractArray{T, D}
end

"""
    Struct to hold a field
    in space S, dimension D and eltype T.
"""
struct Operator{S, D, T}
    space::S
    value::AbstractArray{T, D}
end

struct ProductSpace{T1, T2}
    S1::T1
    S2::T2
end

mutable struct Parameters{T}
    npoints::NTuple{2, Int}
    npatchs::NTuple{2, Int}
    ubounds::NTuple{2, T}
    vbounds::NTuple{2, T}
    nfields::Int
end


struct ChebyshevGL{Tag, N, T} <: Space{Tag, N} 
    min::T
    max::T
    ChebyshevGL{Tag, N, T}(min, max) where {Tag, N, T} = max > min ? new{Tag, N, T}(min, max) : error("bounds are out of order")
end

struct Chebyshev{Tag, N, T} <: Space{Tag, N} 
    min::T
    max::T
    Chebyshev{Tag, N, T}(min, max) where {Tag, N, T} = max > min ? new{Tag, N, T}(min, max) : error("bounds are out of order")
end

Cardinal{Tag, N, T} = Union{ChebyshevGL{Tag, N, T}, 
                              Chebyshev{Tag, N, T}} 

