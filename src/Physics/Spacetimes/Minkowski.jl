#--------------------------------------------------------------------
# DNSS.jl
# Soham 09-2019
# Set initial data for Minkowski spacetime
#--------------------------------------------------------------------

function computeUboundary(PS::ProductSpace{S1, S2})::NTuple{3, Field{S2}} where {S1, S2}
    a = Field(PS. (u,v)->1)
    η = Field(PS, (u,v)->(v-u)/2)
    ϕ = Field(PS, (u,v)->0)
    return extractUboundary(a, η, ϕ)
end

function computeVboundary(PS::ProductSpace{S1, S2})::NTuple{3, Field{S1}} where {S1, S2}
    a = Field(PS. (u,v)->1)
    η = Field(PS, (u,v)->(v-u)/2)
    ϕ = Field(PS, (u,v)->0)
    return extractVboundary(a, η, ϕ)
end
