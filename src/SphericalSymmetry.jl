#--------------------------------------------------------------------
# DNSS.jl
# Soham 04-2022
# Functions for spherical symmetry
#--------------------------------------------------------------------

"""
 Compute initial data for η, given a and ψ
"""
function computeη(U::NTuple{3,Field{S}})::Field{S} where {S<:Space{Tag}} where {Tag}
    PS = first(U).space
    D = derivative(PS)
    I = identity(PS)
    (a, η, ϕ) = U
    B = incomingboundary(PS) + outgoingboundary(PS)
    A = D*D - (2/a)*(D*a)*D + 4pi*(D*ϕ)^2*I
    η = linsolve((I - B)*A + B, B*η)
    return η
end

"""
    Compute the Ricci scalar
"""
function R(U::NTuple{3, Field{S}})::Field{S} where {S}
    (a, η, ϕ) = U
    DU, DV = derivative(first(U).space)
    return - (8 * (DU * ϕ) * (DV * ϕ)) / a^2  
end

"""
    Compute the Hawking mass
"""
function M(U::NTuple{3, Field{S}})::Field{S} where {S}
    (a, η, ϕ) = U
    DU, DV = derivative(first(U).space)
    return (η / 2) * (1 + (4 * (DU * η ) * (DV * η)) / η^2) 
end

"""
    Compute the expansion
    TODO: Check if the sign is right.
"""
function expansion(U::NTuple{3, Field{S}})::Field{S} where {S}
    (a, η, ϕ) = U
    DU, DV = derivative(first(U).space)
    return (DV * η) 
end

function constraints(U::NTuple{3, Field})::NTuple{1, Field}
    PS = first(U).space
    DU, DV = derivative(PS)
    (a, η, ϕ) = U
    C1 = DU*(DU*η) - (2/a)*(DU*a)*(DU*η) + (4pi*η)*(DU*ϕ)^2
    C2 = DV*(DV*η) - (2/a)*(DV*a)*(DV*η) + (4pi*η)*(DV*ϕ)^2
    return (sqrt(C1*C1 + C2*C2), )
end

function lineconstraints(U::NTuple{3, Field})::Field
    PS = first(U).space
    DU = derivative(PS)
    (a, η, ϕ) = U
    return DU*(DU*η) - (2/a)*(DU*a)*(DU*η) + (4pi*η)*(DU*ϕ)^2
end

function residuals(U::NTuple{3, Field})::NTuple{1, Field}
    PS = first(U).space
    DU, DV = derivative(PS)
    (a, η, ϕ) = U
    F1 = DU*(DV*a) - (1/a)*(DU*a)*(DV*a) + (a/η)*(DU*(DV*η)) + (4pi*a)*(DU*ϕ)*(DV*ϕ)
    F2 = DU*(DV*η) + (1/η)*(DU*η)*(DV*η) + (1/4)*(1/η)*(a^2)
    F3 = DU*(DV*ϕ) + (1/η)*(DU*η)*(DV*ϕ) + (1/η)*(DV*η)*(DU*ϕ)
    return (sqrt(F1*F1 + F2*F2 + F3*F3), )
end
