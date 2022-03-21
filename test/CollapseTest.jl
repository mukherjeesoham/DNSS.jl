#--------------------------------------------------------------------
# DNSS.jl
# Soham 03-2022
# Simulate a self-gravitating scalar field
#--------------------------------------------------------------------
# TODO: Try solving for a instead of η
# TODO: Check how you computed the error for the Minkowski case. You might 
# need factors of N and h in the error.

using NLsolve, Printf

#----------------------------------------------
# Compute initial data
#----------------------------------------------
function phi(u::T , v::T)::T where {T<:Number}
    p = 1e-6
    return p*exp(-u^2/0.1) + p*exp(-v^2/0.1)
end

function computeη(a::Field{S}, η::Field{S}, ϕ::Field{S})::Field{S} where {S<:Space{Tag}} where {Tag}
    D = derivative(a.space)
    I = identity(a.space)
    B = incomingboundary(a.space) + outgoingboundary(a.space)
    A = D*D - (2/a)*(D*a)*D + 4pi*(D*ϕ)^2*I
    η = linsolve((I - B)*A + B, B*η)
    @assert norm(D*(D*η) - (2/a)*(D*a)*(D*η) + (4pi*η)*(D*ϕ)^2) < 1e-10
    return η
end

function computeUboundary(PS::ProductSpace{S1, S2})::NTuple{3, Field{S2}} where {S1, S2}
    a0 = extractUboundary(Field(PS, (u,v)->1), :incoming) 
    η0 = extractUboundary(Field(PS, (u,v)->(v-u)/2), :incoming)
    ϕ0 = extractUboundary(Field(PS, (u,v)->phi(u,v)), :incoming)
    return (a0, computeη(a0, η0, ϕ0), ϕ0)
end

function computeVboundary(PS::ProductSpace{S1, S2})::NTuple{3, Field{S1}} where {S1, S2}
    a0 = extractVboundary(Field(PS, (u,v)->1), :incoming) 
    η0 = extractVboundary(Field(PS, (u,v)->(v-u)/2), :incoming)
    ϕ0 = extractVboundary(Field(PS, (u,v)->phi(u,v)), :incoming)
    return (a0, computeη(a0, η0, ϕ0), ϕ0)
end

#----------------------------------------------
# Compute function on each patch 
#----------------------------------------------
function computePatch(PS::ProductSpace{S1, S2}, ubnd::NTuple{3, Field{S2}}, vbnd::NTuple{3, Field{S1}})::NTuple{3, Field{ProductSpace{S1, S2}}} where {S1, S2}
    B = incomingboundary(PS)
    I = identity(PS)
    DU, DV = derivative(PS)
    (abnd, ηbnd, ϕbnd) = combineUVboundary(ubnd, vbnd, :incoming)

    function residual(U::NTuple{3, Field{S}})::NTuple{3, Field{S}} where {S}
        (a, η, ϕ) = U
        F1 = DU*(DV*a) - (1/a)*(DU*a)*(DV*a) + (a/η)*(DU*(DV*η)) + (4pi*a)*(DU*ϕ)*(DV*ϕ)
        F2 = DU*(DV*η) + (1/η)*(DU*η)*(DV*η) + (1/4)*(1/η)*(a^2)
        F3 = DU*(DV*ϕ) + (1/η)*(DU*η)*(DV*ϕ) + (1/η)*(DV*η)*(DU*ϕ)

        # TODO: This is SLOW. Speed this up by replacing the residual in place. 
        F1 = (I - B) * F1 + B*(a - abnd)
        F2 = (I - B) * F2 + B*(η - ηbnd)
        F3 = (I - B) * F3 + B*(ϕ - ϕbnd)
        return (F1, F2, F3)
    end

    function initialguess(PS::ProductSpace{S1, S2})::NTuple{3, Field{ProductSpace{S1, S2}}} where {S1, S2}
        a0 = Field(PS, (u,v)->1) 
        η0 = Field(PS, (u,v)->1) 
        ϕ0 = Field(PS, (u,v)->1)
        return (a0, η0, ϕ0)
    end

    function f!(fvec::Array{T,1}, x::Array{T,1}) where {T}
        fvec[:] = reshapeFromTuple(residual(reshapeToTuple(PS, 3, x)))
    end

    # Call the non-linear solver
    U0 = reshapeFromTuple(initialguess(PS))
    U = nlsolve(f!, U0; method=:trust_region, autodiff=:forward, show_trace=false, ftol=1e-9, iterations=100)
    return reshapeToTuple(PS, 3, U.zero)
end

#----------------------------------------------
# Constraint equations and convergence 
#----------------------------------------------
function constraints(U::NTuple{3, Field{S}})::Field{S} where {S}
    (a, η, ϕ) = U
    DU, DV = derivative(a.space)
    C1 = DU*(DU*η) - (2/a)*(DU*a)*(DU*η) + (4pi*η)*(DU*ϕ)^2
    C2 = DV*(DV*η) - (2/a)*(DV*a)*(DV*η) + (4pi*η)*(DV*ϕ)^2
    return sqrt(C1*C1 + C2*C2)
end

function constraints(AoT::Matrix{Union{NTuple{3, Field}, ProductSpace}})::Matrix{Field} 
    return [constraints(AoT[index]) for index in CartesianIndices(AoT)]
end

function rmse(AoT::Matrix{Union{NTuple{3, Field}, ProductSpace}})::Real
    return norm(norm.(constraints(AoT)))
end

function pconv(min, max)
    println("Testing p convergence")
    n_ = collect(min:max)
    l_ = zeros(size(n_))
    for index in CartesianIndices(n_) 
        n = n_[index]
        l_[index] = rmse(distribute(Parameters((n, n), (1,1), urange, vrange, nfields), computePatch, computeUboundary, computeVboundary))
        @printf("  n = %i, rmse = %e\n", n_[index], l_[index])
    end
    plotpconv(n_, l_, "../output/vaidya-constraints-pconv.pdf")
end

function hconv(min, max)
    println("Testing h convergence")
    n_ = collect(min:max)
    l_ = zeros(size(n_))
    for index in CartesianIndices(n_) 
        np = 2^n_[index]
        n_[index] = np
        l_[index] = rmse(distribute(Parameters((4, 4), (np,np), urange, vrange, nfields), computePatch, computeUboundary, computeVboundary))
        @printf("  n = %2i, rmse = %e\n", n_[index], l_[index])
    end
    plothconv(n_, l_, "../output/minkowski-constraints-hconv.pdf")
end
#-----------------------------------------
# Setup simulation grid parameters and
# call distribute.
#-----------------------------------------
npoints  = (4, 4)
npatches = (2, 2) 
urange   = (-1.0, 1.0)
vrange   = (-1.0, 1.0) 
nfields  = 3
params   = adjust(Parameters(npoints, npatches, urange, vrange, nfields), 1e-2)

#-----------------------------------------
# Cover the grid
#-----------------------------------------
U = distribute(params, computePatch, computeUboundary, computeVboundary)
@show rmse(U)

#-----------------------------------------
# Plot the solutions and constraints
#-----------------------------------------
contourf(extract(U, 3), 20, "../output/vaidya-psi.pdf")
contourf(constraints(U), 20, "../output/vaidya-constraints-vector-norm.pdf")

#-----------------------------------------
# Check and plot convergence
#-----------------------------------------
pconv(4, 8)
hconv(0, 5)

