#--------------------------------------------------------------------
# DNSS.jl
# Soham 09-2019
# Set initial data for Schwarzschild spacetime
# NOTE: Explicitly set mass of Schwarzschild. 
# You can choose to either use the initial data solver, 
# or set the initial data using analytic expressions
#--------------------------------------------------------------------

export computeUboundary, computeVboundary, excision
export setSchwarzschild

function computeUboundary(PS::ProductSpace{S1, S2})::NTuple{3, Field{S2}} where {S1, S2}
    M  = 1.0
    η0 = Field(PS, (u,v)->find_r_of_UV(u,v,M))
    ϕ0 = Field(PS, (u,v)->0)
    f0 = ((16*M^3)/η0)*exp(-η0/2M)
    a0 = -sqrt(2*f0) 
    return computeUboundary((a0, η0, ϕ0))
end

function computeVboundary(PS::ProductSpace{S1, S2})::NTuple{3, Field{S1}} where {S1, S2}
    M  = 1.0
    η0 = Field(PS, (u,v)->find_r_of_UV(u,v,M))
    ϕ0 = Field(PS, (u,v)->0)
    f0 = ((16*M^3)/η0)*exp(-η0/2M)
    a0 = -sqrt(2*f0) 
    return computeVboundary((a0, η0, ϕ0))
end

function excision(PS::ProductSpace{S1, S2})::Bool where {S1, S2}
    r = Field(PS, (u,v)->v-u)
    v = Field(PS, (u,v)->v)
    return any(value(r) .>= eptype(value(r))(1)) || any(value(v) .< eps(eltype(value(r)))) 
end

function initialguess(PS::ProductSpace{S1, S2})::NTuple{3, Field{ProductSpace{S1, S2}}} where {S1, S2}
    a = Field(PS, (u,v)->1)
    η = Field(PS, (u,v)->(v-u)/2)
    ϕ = Field(PS, (u,v)->0)
    return (a, η, ϕ)
end

function extractUboundary(PS::ProductSpace{S1, S2}, boundarytype::Symbol)::NTuple{3, Field{S2}} where {S1, S2}
    η0 = Field(PS, (u,v)->missing)
    ϕ0 = Field(PS, (u,v)->missing)
    a0 = Field(PS, (u,v)->missing)
    return extractUboundary((a, η, ϕ), boundarytype)
end

function extractVboundary(PS::ProductSpace{S1, S2}, boundarytype::Symbol)::NTuple{3, Field{S1}} where {S1, S2}
    η0 = Field(PS, (u,v)->missing)
    ϕ0 = Field(PS, (u,v)->missing)
    a0 = Field(PS, (u,v)->missing)
    return extractVboundary((a, η, ϕ), boundarytype)
end
function setSchwarzschild(PS::ProductSpace{S1, S2})::NTuple{3, Field{ProductSpace{S1, S2}}} where {S1, S2}
    M  = 1.0
    η0 = Field(PS, (u,v)->find_r_of_UV(u,v,M))
    ϕ0 = Field(PS, (u,v)->0)
    f0 = ((16*M^3)/η0)*exp(-η0/2M)
    a0 = -sqrt(2*f0) 
    return (a0, η0, ϕ0)
end
