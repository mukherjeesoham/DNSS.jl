#--------------------------------------------------------------------
# DNSS.jl
# Soham 08-2019
# Test scalar field on Minkowski spacetime
#--------------------------------------------------------------------
# TODO: Evaluate the residual at twice the number of grid points.
# TODO: Change the function 

using Test, NLsolve, LinearAlgebra, Random
using LaTeXStrings, Printf, PyPlot

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
function nonlin(PS::ProductSpace{S1, S2}, ubnd::NTuple{1,Field{S2}}, vbnd::NTuple{1,Field{S1}})::NTuple{1, Field{ProductSpace{S1, S2}}} where {S1, S2}

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

#------------------------------------------------------------
# Set up computation on each patch using a linear solver
#------------------------------------------------------------
function lin(PS::ProductSpace{S1, S2}, ubnd::NTuple{1,Field{S2}}, vbnd::NTuple{1,Field{S1}})::NTuple{1, Field{ProductSpace{S1, S2}}} where {S1, S2}

    # Define all the operators we use local to the patch. 
    B = incomingboundary(PS)
    DU, DV = derivative(PS)
    bnd = combineUVboundary(ubnd, vbnd, :incoming)
    L   = 2*DU*DV 
    
    # Using the linearity of the system to solve the system
    # L u  = 0 # Evolution equations
    # B u =  b # Boundary conditions
    # (L + B) u = b
    U = linsolve(L + B,  first(bnd)) 
    return (U, )
end

#-----------------------------------------
# Compare with analytic solution
#-----------------------------------------

function deltapsi(U::NTuple{1, Field})::NTuple{1, Field} where {S}
    # Compute the error at twice the number of points
    psi = project(first(U), prolongate(first(U).space))
    psiexact = Field(psi.space, psibar)
    return (psi - psiexact, ) 
end

#-----------------------------------------
# h and p convergence functions 
#-----------------------------------------
function pconv(min, max)
    println("Testing p convergence")
    n_ = collect(min:max)
    l_ = zeros(size(n_))
    for index in CartesianIndices(n_) 
        n = n_[index]
        l_[index]  =  rmse(extract(map(deltapsi, distribute(Parameters((n, n), (1,1), urange, vrange, nfields), compute, ubnd, vbnd)), 1))
        @printf("  n = %i,  rmse = %e\n", n_[index], l_[index])
    end
    plotpconv(n_, l_, "../output/minkowski-psi-pconv.pdf")
end

function hconv(min, max)
    println("Testing h convergence")
    n_ = collect(min:max)
    l_ = zeros(size(n_))
    for index in CartesianIndices(n_) 
        np = 2^n_[index]
        n_[index] =  np
        l_[index] =  rmse(extract(map(deltapsi, distribute(Parameters((4, 4), (np, np), urange, vrange, nfields), compute, ubnd, vbnd)), 1))
        l0 = index[1] > 1 ? l_[index[1] - 1] : 1
        @printf("  n = %3i, rmse = %e\n", n_[index], l0 / l_[index])
    end
    plothconv(n_, l_, "../output/minkowski-psi-hconv.pdf")
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

# Choose whether to use linear or non-linear solver.
compute = nonlin

#-----------------------------------------
# Compute the scalar field over the grid
#-----------------------------------------
ϕ_ = distribute(Grid, compute, ubnd, vbnd)
@test rmse(extract(map(deltapsi, ϕ_), 1)) < 1e-9

#-----------------------------------------
# Plot the scalar field and the error
#-----------------------------------------
contourf(extract(ϕ_, 1), 20, "../output/minkowski-psi.pdf")
contourf(extract(map(deltapsi, ϕ_), 1), 20, "../output/minkowski-psi-error.pdf")

#-----------------------------------------
# Now test h-p convergence
#-----------------------------------------
pconv(4, 30)
hconv(0, 7)
