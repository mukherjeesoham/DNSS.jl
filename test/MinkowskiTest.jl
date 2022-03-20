#--------------------------------------------------------------------
# DNSS.jl
# Soham 08-2019
# Test scalar field on Minkowski spacetime
#--------------------------------------------------------------------
# FIXME: Why are we converging instantly? 
# TODO: Try other functions. Since the convergence is suspicious. 
# TODO: Check with the linear solver
# TODO: Check how we were computing the error there. 

using NLsolve, LinearAlgebra, Random
using PyPlot, LaTeXStrings, Printf

#-----------------------------------------
# Set up initial conditions on u0 and v0
#-----------------------------------------
function psibar(u::Real, v::Real)::Real
    return exp(-u^2) * sin(pi*u) + exp(-v^2) * sin(pi*v)
end

function ubnd(PS::ProductSpace{S1, S2})::NTuple{1, Field{S2}} where {S1, S2}
    return (extractUboundary(Field(PS, psibar), :incoming), )
end

function vbnd(PS::ProductSpace{S1, S2})::NTuple{1, Field{S1}} where {S1, S2}
    return (extractVboundary(Field(PS, psibar), :incoming), )
end

#------------------------------------------------------------
# Set up computation on each patch using a non-linear solver
#------------------------------------------------------------
function compute(PS::ProductSpace{S1, S2}, ubnd::NTuple{1,Field{S2}}, vbnd::NTuple{1,Field{S1}})::NTuple{1, Field{ProductSpace{S1, S2}}} where {S1, S2}

    # Define all the operators we use local to the patch. 
    B = incomingboundary(PS)
    I = identity(PS)
    DU, DV = derivative(PS)
    psibnd = combineUVboundary(ubnd, vbnd, :incoming)

    # Compute the residual. This is slow but works. Sets the residual in the
    # interior as well as the boundary.
    function residual(psi::NTuple{1, Field{S}})::NTuple{1, Field{S}} where {S}
        psi = first(psi)
        F = (I - B) * (DU*(DV*psi)) + B * (psi - first(psibnd))
        return (F,)
    end

    # Residual to pass to the non-linear solver in NLsolve
    function f!(fvec::Array{T,1}, x::Array{T,1}) where {T}
        fvec[:]  = reshapeFromTuple(residual(reshapeToTuple(PS, 1, x)))
    end

    # Specify an initial guess
    function initialguess(PS::ProductSpace{S1, S2})::NTuple{1, Field{ProductSpace{S1, S2}}} where {S1, S2}
        return (Field(PS, (u,v)->1), ) 
    end

    # Compute the initial guess to pass to the solver
    U0 = reshapeFromTuple(initialguess(PS))

    # Call the non-linear solver
    U = nlsolve(f!, U0; method=:trust_region, autodiff=:forward, show_trace=false, ftol=1e-8, iterations=100)
    return reshapeToTuple(PS, 1, U.zero)
end

#-----------------------------------------
# Compare with analytic solution
#-----------------------------------------
function L2(psi::Field{S})::Real where {S}
    return norm(psi - Field(psi.space, psibar))
end

function L2(AoF::Matrix{Field})::Real where {S}
    return norm([L2(AoF[index]) for index in CartesianIndices(AoF)])
end

#-----------------------------------------
# h and p convergence helper 
#-----------------------------------------
function pconv(min, max)
    println("Testing p convergence")
    n_ = collect(min:max)
    l_ = zeros(size(n_))
    for index in CartesianIndices(n_) 
        n = n_[index]
        l_[index]  =  L2(extract(distribute(Parameters((n, n), (1,1), urange, vrange, nfields), compute, ubnd, vbnd), 1))
        @printf("  n = %i, L2 error = %e\n", n_[index], l_[index])
    end
    plotpconv(n_, l_)
end

function hconv(min, max)
    println("Testing h convergence")
    println("Testing p convergence")
    n_ = collect(min:max)
    l_ = zeros(size(n_))
    for index in CartesianIndices(n_) 
        np = 2^n_[index]
        l_[index] =  L2(extract(distribute(Parameters((4, 4), (np, np), urange, vrange, nfields), compute, ubnd, vbnd), 1))
        n_[index] =  np
        @printf("  n = %2i, L2 error = %e\n", n_[index], l_[index])
    end
    plotpconv(n_, l_)
end

#-----------------------------------------
# Setup simulation grid parameters 
#-----------------------------------------
npoints = (12, 12)
npatches = (4, 4) 
urange = (-1.0, 1.0)
vrange = (-1.0, 1.0) 
nfields = 1
Grid = Parameters(npoints, npatches, urange, vrange, nfields)

#-----------------------------------------
# Compute the scalar field over the grid
#-----------------------------------------
ϕ_ = distribute(Grid, compute, ubnd, vbnd)
ϕ  = extract(ϕ_, 1)
@test L2(ϕ) < 1e-12

#-----------------------------------------
# Plot the scalar field
#-----------------------------------------
contourf(ϕ, 20)

#-----------------------------------------
# Now test h-p convergence
#-----------------------------------------
pconv(2, 8)
hconv(0, 5)
