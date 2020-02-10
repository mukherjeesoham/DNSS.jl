#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 08-2018
# Add basis functions
#--------------------------------------------------------------------

export Cardinal, PointSpace
export ChebyshevGL, LegendreGL, FourierGEP, FourierCEP, Chebyshev

struct ChebyshevGL{Tag, N, T} <: Space{Tag, N} 
    min::T
    max::T
    ChebyshevGL{Tag, N, T}(min, max) where {Tag, N, T} = max > min ? new{Tag, N, T}(min, max) : error("bounds are out of order")
end

struct LegendreGL{Tag, N, T} <: Space{Tag, N} 
    min::T
    max::T
    LegendreGL{Tag, N, T}(min, max) where {Tag, N, T} = max > min ? new{Tag, N, T}(min, max) : error("bounds are out of order")
end

struct FourierGEP{Tag, N, T} <: Space{Tag, N} 
    min::T
    max::T
    FourierGEP{Tag, N, T}(min, max) where {Tag, N, T} = new{Tag, N, T}(0, 2*pi) 
end

struct FourierCEP{Tag, N, T} <: Space{Tag, N} 
    min::T
    max::T
    FourierCEP{Tag, N, T}(min, max) where {Tag, N, T} = new{Tag, N, T}(0, pi) 
end

struct Chebyshev{Tag, N, T} <: Space{Tag, N} 
    min::T
    max::T
    Chebyshev{Tag, N, T}(min, max) where {Tag, N, T} = max > min ? new{Tag, N, T}(min, max) : error("bounds are out of order")
end

# FIXME: Chebyshev is in Cardinal
Cardinal{Tag, N, T} = Union{ChebyshevGL{Tag, N, T}, 
                            LegendreGL{Tag, N, T}, 
                            FourierGEP{Tag, N, T}, 
                            FourierCEP{Tag, N, T}, 
                            Chebyshev{Tag, N, T}} 
