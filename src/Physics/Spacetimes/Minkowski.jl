#--------------------------------------------------------------------
# DNSS.jl
# Soham 09-2019
# Set initial data for Minkowski spacetime
#--------------------------------------------------------------------

import Base.Threads.@spawn
export computeUboundary, computeVboundary, excision
export extractUboundary, extractVboundary
export FextractUboundary, FextractVboundary

function computeUboundary(PS::ProductSpace{S1, S2})::NTuple{3, Field{S2}} where {S1, S2}
    a = Field(PS, (u,v)->1)
    η = Field(PS, (u,v)->(v-u)/2)
    ϕ = Field(PS, (u,v)->0)
    return extractUboundary((a, η, ϕ), :incoming)
end

function computeVboundary(PS::ProductSpace{S1, S2})::NTuple{3, Field{S1}} where {S1, S2}
    a = Field(PS, (u,v)->1)
    η = Field(PS, (u,v)->(v-u)/2)
    ϕ = Field(PS, (u,v)->0)
    return extractVboundary((a, η, ϕ), :incoming)
end

function excision(PS::ProductSpace{S1, S2})::Bool where {S1, S2}
    r = Field(PS, (u,v)->v-u)
    return any(r.value .<= eps(eltype(r.value)))
end

function initialguess(PS::ProductSpace{S1, S2})::NTuple{3, Field{ProductSpace{S1, S2}}} where {S1, S2}
    a = Field(PS, (u,v)->1)
    η = Field(PS, (u,v)->(v-u)/2)
    ϕ = Field(PS, (u,v)->0)
    return (a + η, 1 + η, ϕ)
end

function extractUboundary(PS::ProductSpace{S1, S2}, boundarytype::Symbol)::NTuple{3, Field{S2}} where {S1, S2}
    a = Field(PS, (u,v)->1)
    η = Field(PS, (u,v)->(v-u)/2)
    ϕ = Field(PS, (u,v)->0)
    return extractUboundary((a, η, ϕ), boundarytype)
end

function extractVboundary(PS::ProductSpace{S1, S2}, boundarytype::Symbol)::NTuple{3, Field{S1}} where {S1, S2}
    a = Field(PS, (u,v)->1)
    η = Field(PS, (u,v)->(v-u)/2)
    ϕ = Field(PS, (u,v)->0)
    return extractVboundary((a, η, ϕ), boundarytype)
end
