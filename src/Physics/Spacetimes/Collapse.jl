#--------------------------------------------------------------------
# DNSS.jl
# Soham 09-2019
# Set initial data for Minkowski spacetime with a scalar field
#--------------------------------------------------------------------

export computeUboundary, computeVboundary, initialguess, excision

function pulse(u::T , v::T)::T where {T<:Number}
    p = 1e-6
    return p*exp(-u^2/0.1) + p*exp(-v^2/0.1)
end

function computeUboundary(PS::ProductSpace{S1, S2})::NTuple{3, Field{S2}} where {S1, S2}
    a0 = Field(PS, (u,v)->1) 
    η0 = Field(PS, (u,v)->(v-u)/2) 
    ϕ0 = Field(PS, (u,v)->pulse(u,v))
    return computeUboundary((a0, η0, ϕ0))
end

function computeVboundary(PS::ProductSpace{S1, S2})::NTuple{3, Field{S1}} where {S1, S2}
    a0 = Field(PS, (u,v)->1) 
    η0 = Field(PS, (u,v)->(v-u)/2) 
    ϕ0 = Field(PS, (u,v)->pulse(u,v))
    return computeVboundary((a0, η0, ϕ0))
end

function initialguess(PS::ProductSpace{S1, S2})::NTuple{3, Field{ProductSpace{S1, S2}}} where {S1, S2}
    a0 = Field(PS, (u,v)->1) 
    η0 = Field(PS, (u,v)->(v-u)/2) 
    ϕ0 = Field(PS, (u,v)->pulse(u,v))
    return (a0, η0, ϕ0)
end

function excision(PS::ProductSpace{S1, S2})::Bool where {S1, S2}
    return false
end
