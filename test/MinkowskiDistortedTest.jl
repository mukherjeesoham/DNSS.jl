#--------------------------------------------------------------------
# DNSS.jl
# Soham 08-2019
# Test scalar field on Minkowski spacetime
#--------------------------------------------------------------------

using Test, NLsolve, LinearAlgebra, Random
using LaTeXStrings, Printf, PyPlot

#-------------------------------------------
# Set up initial conditions 
#-------------------------------------------
function twist(uv::NTuple{2, Real})::NTuple{2, Real}
    u, v = uv
    Ω =  (pi/8) * cospi(u/2) * cospi(v/2)
    ũ =  u*cos(Ω) + v*sin(Ω)
    ṽ = -u*sin(Ω) + v*cos(Ω)
    return (ũ, ṽ)
end

function psibar(u::Real, v::Real)::Real
    (u,v) = twist((u,v))
    return exp(-u^2) * sin(pi*u) + exp(-v^2) * sin(pi*v)
end

function ubnd(PS::ProductSpace{S1, S2})::NTuple{1, Field{S2}} where {S1, S2}
    return (extractUboundary(Field(PS, psibar), :incoming), )
end

function vbnd(PS::ProductSpace{S1, S2})::NTuple{1, Field{S1}} where {S1, S2}
    return (extractVboundary(Field(PS, psibar), :incoming), )
end

function Power(u::Field{S}, p::Int)::Field{S} where {S}
    return u^p
end

#------------------------------------------------------------
# Set up computation on each patch using a non-linear solver
#------------------------------------------------------------
function compute(PS::ProductSpace{S1, S2}, ubnd::NTuple{1,Field{S2}}, vbnd::NTuple{1,Field{S1}})::NTuple{1, Field{ProductSpace{S1, S2}}} where {S1, S2}

    # Define all the operators we use local to the patch. 
    B = incomingboundary(PS)
    I = identity(PS)
    DU, DV = derivative(PS)
    Psibnd = combineUVboundary(ubnd, vbnd, :incoming)

    P = Field(PS, (u,v)->twist((u,v))[1])
    Q = Field(PS, (u,v)->twist((u,v))[2])
    
      duP = DU * P
      dvP = DV * P
    duduP = DU * (DU * P)
    dvdvP = DV * (DV * P)
    dudvP = DU * (DV * P)
    
      duQ = DU * Q
      dvQ = DV * Q
    duduQ = DU * (DU * Q)
    dvdvQ = DV * (DV * Q)
    dudvQ = DU * (DV * Q)

    function residual(psi::NTuple{1, Field{S}})::NTuple{1, Field{S}} where {S}
        Psi = first(psi)
        duPsi = DU * Psi
        dvPsi = DV * Psi
        duduPsi = DU * (DU * Psi)
        dvdvPsi = DV * (DV * Psi)
        dudvPsi = DU * (DV * Psi)
    
        # Compute the expression using Mathematica.  FIXME: The residual is not
        # zero even for an exact solution. Fix the computation in
        # the Mathematica notebook.

        F = ((2*(dvP*(dudvPsi*Power(duQ,3)*Power(dvP,2) + dudvQ*duPsi*duQ*Power(dvP,2)*(-duQ + dvP) + dudvP*Power(duQ,4)*dvPsi - dudvP*Power(duQ,3)*dvP*dvPsi -
          dudvP*duPsi*Power(duQ,3)*dvQ + 2*dudvP*duPsi*Power(duQ,2)*dvP*dvQ + duduQ*duPsi*duQ*Power(dvP,2)*dvQ - duduPsi*Power(duQ,2)*Power(dvP,2)*dvQ -
          2*duduQ*duPsi*Power(dvP,3)*dvQ + duduPsi*duQ*Power(dvP,3)*dvQ - duduP*Power(duQ,3)*dvPsi*dvQ + 2*duduP*Power(duQ,2)*dvP*dvPsi*dvQ -
          duduP*duQ*Power(dvP,2)*dvPsi*dvQ + duduP*duPsi*Power(duQ,2)*Power(dvQ,2) - 3*duduP*duPsi*duQ*dvP*Power(dvQ,2) + 2*duduP*duPsi*Power(dvP,2)*Power(dvQ,2)) -
          Power(duP,3)*dvQ*(-(dvdvQ*dvP*dvPsi) + dvQ*(dvdvPsi*dvP - dudvQ*dvPsi + dvdvP*dvPsi + dudvPsi*dvQ) + duQ*(dvdvQ*dvPsi - dvdvPsi*dvQ)) +
          Power(duP,2)*(Power(duQ,2)*dvQ*(-4*dudvQ*dvPsi + 3*dvdvP*dvPsi + 2*dudvPsi*dvQ) - Power(duQ,2)*dvdvQ*(3*dvP*dvPsi + duPsi*dvQ) +
          duQ*dvdvQ*dvP*(dvP*dvPsi + 2*duPsi*dvQ) + Power(duQ,3)*(2*dvdvQ*dvPsi - dvdvPsi*dvQ) -
          duQ*dvQ*(-(dvdvPsi*Power(dvP,2)) - 5*dudvQ*dvP*dvPsi + dvdvP*dvP*dvPsi - 2*dudvQ*duPsi*dvQ + 2*duPsi*dvdvP*dvQ + dudvPsi*dvP*dvQ - 2*duduQ*dvPsi*dvQ +
             2*dudvP*dvPsi*dvQ + duduPsi*Power(dvQ,2)) + dvP*dvQ*(-3*dudvQ*dvP*dvPsi + dvQ*(2*dudvPsi*dvP - 2*duduQ*dvPsi + 2*dudvP*dvPsi + duduPsi*dvQ)) -
          duPsi*dvQ*(dvdvQ*Power(dvP,2) + dvQ*(2*dudvQ*dvP - 2*dvdvP*dvP + duduQ*dvQ - dudvP*dvQ))) +
            duP*(Power(duQ,4)*(dvdvPsi*dvP - 2*dvdvP*dvPsi) + Power(duQ,3)*
           (-(duPsi*dvdvQ*dvP) - dvdvPsi*Power(dvP,2) + dvdvP*dvP*dvPsi + 2*duPsi*dvdvP*dvQ - 2*dudvPsi*dvP*dvQ + 3*dudvP*dvPsi*dvQ) -
          dvP*dvQ*(-3*dudvQ*duPsi*Power(dvP,2) - duduQ*dvP*(2*dvP*dvPsi + 3*duPsi*dvQ) +
             dvQ*(4*dudvP*duPsi*dvP + duduPsi*Power(dvP,2) + duduP*dvP*dvPsi + duduP*duPsi*dvQ)) +
          Power(duQ,2)*(2*duPsi*dvdvQ*Power(dvP,2) + 2*dudvQ*Power(dvP,2)*dvPsi + duPsi*dvQ*(2*dudvQ*dvP - 2*dvdvP*dvP - 3*dudvP*dvQ) +
             dvQ*(dudvPsi*Power(dvP,2) - 5*dudvP*dvP*dvPsi + duduPsi*dvP*dvQ - duduP*dvPsi*dvQ)) -
          duQ*(duPsi*dvdvQ*Power(dvP,3) + dudvQ*Power(dvP,3)*dvPsi + duPsi*dvQ*(5*dudvQ*Power(dvP,2) + dvQ*(duduQ*dvP - 5*dudvP*dvP - duduP*dvQ)) +
               2*dvP*dvQ*(dudvPsi*Power(dvP,2) + dvPsi*(duduQ*dvP - dudvP*dvP - duduP*dvQ))))))/Power(duQ*dvP - duP*dvQ,3))

        F = (I - B) * F + B * (Psi - first(Psibnd))
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
    U = nlsolve(f!, U0; method=:trust_region, autodiff=:forward, show_trace=true, ftol=1e-8, iterations=100)
    return reshapeToTuple(PS, 1, U.zero)
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
    plotpconv(n_, l_, "../output/minkowski-distorted-psi-pconv.pdf")
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
    plothconv(n_, l_, "../output/minkowski-distorted-psi-hconv.pdf")
end

#-----------------------------------------
# Setup simulation grid parameters 
#-----------------------------------------
npoints = (20, 20)
npatches = (1, 1) 
urange = (-1.0, 1.0)
vrange = (-1.0, 1.0) 
nfields = 1
Grid = Parameters(npoints, npatches, urange, vrange, nfields)

#-----------------------------------------
# Compute the scalar field over the grid
#-----------------------------------------
ϕ_ = distribute(Grid, compute, ubnd, vbnd)
# @test rmse(extract(map(deltapsi, ϕ_), 1)) < 1e-9

#-----------------------------------------
# Plot the scalar field and the error
#-----------------------------------------
contourf(extract(ϕ_, 1), 20, "../output/minkowski-distorted-psi.pdf")
# contourf(extract(map(deltapsi, ϕ_), 1), 20, "../output/minkowski-distorted-psi-error.pdf")

#-----------------------------------------
# Now test h-p convergence
#-----------------------------------------
# pconv(4, 30)
# hconv(0, 7)
