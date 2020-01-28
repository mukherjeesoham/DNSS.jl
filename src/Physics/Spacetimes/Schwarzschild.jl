#--------------------------------------------------------------------
# DNSS.jl
# Soham 09-2019
# Set initial data for Schwarzschild spacetime
# NOTE: Explicitly set mass of Schwarzschild. 
# You can choose to either use the initial data solver, 
# or set the initial data using analytic expressions
#--------------------------------------------------------------------

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
