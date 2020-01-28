#--------------------------------------------------------------------
# DNSS.jl
# Soham 09-2019
# Set initial data for Minkowski spacetime
#--------------------------------------------------------------------

function pulse(p::T, u::T , v::T)::T where {T<:Number}
    if v >= 4 && v <= 6
         return p*exp(-(v-5)^2) - 1/exp(1)
    else
        return 0
    end
end

function computeUboundary(PS::ProductSpace{S1, S2})::NTuple{3, Field{S2}} where {S1, S2}
    p0 = 1.0
    a0 = Field(PS, (u,v)->1) 
    η0 = Field(PS, (u,v)->(v-u)/2) 
    ϕ0 = Field(PS, (u,v)->pulse(p0,u,v))
    return computeVboundary((a0, η0, ϕ0))
end

function computeVboundary(PS::ProductSpace{S1, S2})::NTuple{3, Field{S1}} where {S1, S2}
    p0 = 1.0
    a0 = Field(PS, (u,v)->1) 
    η0 = Field(PS, (u,v)->(v-u)/2) 
    ϕ0 = Field(PS, (u,v)->pulse(p0,u,v))
    return computeUboundary((a0, η0, ϕ0))
end

