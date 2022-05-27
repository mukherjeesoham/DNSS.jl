#--------------------------------------------------------------------
# DNSS.jl
# Soham 08-2019
# Test scalar field on Minkowski spacetime
#--------------------------------------------------------------------

using Test, NLsolve, LinearAlgebra, Random, Printf, ForwardDiff
# using LaTeXStrings, Printf, PyPlot

#-----------------------------------------
# Set up initial conditions on u0 and v0
#-----------------------------------------
function psibar(u::Real, v::Real)::Real
    return  exp(-u^2 / σ^2) * sin(pi*u) + exp(-v^2 / σ^2) * sin(pi*v)
end

function computeUboundary(PS::ProductSpace{S1, S2})::NTuple{1, Field{S2}} where {S1, S2}
    return (extractUboundary(Field(PS, psibar), :incoming), )
end

function computeVboundary(PS::ProductSpace{S1, S2})::NTuple{1, Field{S1}} where {S1, S2}
    return (extractVboundary(Field(PS, psibar), :incoming), )
end

#------------------------------------------------------------
# Source term or the change of coordinate functions 
#------------------------------------------------------------
function twist(uv::Vector{T})::Vector{T} where {T<:Real}
    u, v = uv
    ω =  (pi/8) * cospi(u/2) * cospi(v/2)
    ũ =  u*cos(ω) + v*sin(ω)
    ṽ = -u*sin(ω) + v*cos(ω)
    return [ũ, ṽ]
end

function dxdX(u::Real, v::Real, component::Symbol)
    if component == :dUdŨ
        return inv(ForwardDiff.jacobian(twist, [u,v]))[1,1]
    elseif component == :dUdṼ
        return inv(ForwardDiff.jacobian(twist, [u,v]))[1,2]
    elseif component == :dVdŨ
        return inv(ForwardDiff.jacobian(twist, [u,v]))[2,1]
    elseif component == :dVdṼ
        return inv(ForwardDiff.jacobian(twist, [u,v]))[2,2]
    end
end

function source(u::Real, v::Real)::Real
    t = u + v
    r = v - u
    return 20 * sinpi(t) * sinpi(r) *  exp(-r^2 / 0.1)
end

#------------------------------------------------------------
# Set up computation on each patch using a non-linear solver
#------------------------------------------------------------
function computePatch(PS::ProductSpace{S1, S2}, ubnd::NTuple{1,Field{S2}}, vbnd::NTuple{1,Field{S1}})::NTuple{1, Field{ProductSpace{S1, S2}}} where {S1, S2}
    if solver == :nonlin
        # Define all the operators we use local to the patch. 
        B = incomingboundary(PS)
        I = identity(PS)
        DU, DV = derivative(PS)
        psibnd = combineUVboundary(ubnd, vbnd, :incoming)

        dUdŨ   = Field(PS, (u,v)-> dxdX(u,v, :dUdŨ))
        dVdŨ   = Field(PS, (u,v)-> dxdX(u,v, :dVdŨ))
        dUdṼ   = Field(PS, (u,v)-> dxdX(u,v, :dUdṼ))
        dVdṼ   = Field(PS, (u,v)-> dxdX(u,v, :dVdṼ))

        # Compute the residual. This is slow but works. Sets the residual in the
        # interior as well as the boundary.
        function residual(psi::NTuple{1, Field{S}})::NTuple{1, Field{S}} where {S}
            psi   = first(psi)
            # FIXME: The twisted coordinates do not work
            dpsidŨ    = dUdŨ * (DU * psi) + dVdŨ * (DV * psi) 
            ddpsidŨdṼ = dUdṼ * (DU * dpsidŨ) + dVdṼ * (DV * dpsidŨ) 
            F = (I - B) * ddpsidŨdṼ  + B * (psi - first(psibnd))
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

    elseif solver == :lin
        # Define all the operators we use local to the patch. 
        B = incomingboundary(PS)
        DU, DV = derivative(PS)
        bnd = combineUVboundary(ubnd, vbnd, :incoming)
        L   = 2*DU*DV 
        s   = 0 * Field(PS, source)
        b   = first(bnd)
        
        # Using the linearity of the system to solve the system
        # L u  = s # Wave equation with a source 
        # B u =  b # Boundary conditions
        # (L + B) u = b
        U = linsolve(L + B,  s + b) 
        return (U, )
    end
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
function pconv(range::StepRange, nh::Int=1)
    println("Testing p convergence")
    n_ = collect(range)
    l_ = zeros(size(n_))
    for index in CartesianIndices(n_) 
        np = n_[index]
        params = Parameters((np, np), (nh,nh), urange, vrange, 1)
        q = distribute(params, computePatch, computeUboundary, computeVboundary)
        l_[index] = rmse(extract(map(deltapsi, q), 1))
        @printf("  n = %i, rmse = %e\n", n_[index], l_[index])
    end
    return (n_, l_)
end

function hconv(range::StepRange, np::Int=4)
    println("Testing h convergence")
    n_ = collect(range)
    l_ = zeros(size(n_))
    for index in CartesianIndices(n_) 
        nh = n_[index] = 2^n_[index]
        params = Parameters((np, np), (nh, nh), urange, vrange, 1)
        l_[index] = rmse(extract(map(deltapsi, distribute(params, computePatch, computeUboundary, computeVboundary)), 1))
        l0 = index[1] > 1 ? l_[index[1] - 1] : 1
        @printf("  n = %3i, rmse[2h] / rmse[h] = %e\n", n_[index], l0 / l_[index])
    end
    return (n_, l_)
end

#-----------------------------------------
# Setup simulation grid parameters 
#-----------------------------------------
# TODO: Check if we find convergence for different sizes of domains. What
# happens if the domain size increases? Is the solution unstable for larger 
# domain sizes? What we find is that one needs more spacetime elements for
# larger domain sizes to see p-convergence. For a single or fewer elements, 
# the p-convergence is either non-existent, or slow. 
# FIXME: In any case, we saturate around 1e-10 with all the errors. Is this due
# to the floor set by the condition number of the second derivative matrices?
# What happens if we choose to work with BigFloat?
npoints  = (20, 20)
npatches = (1, 1) 
urange   = (-1.0, 1.0)
vrange   = (-1.0, 1.0)
nfields  = 1
Grid     = Parameters(npoints, npatches, urange, vrange, nfields)
solver   = :lin # Use a linear solver instead of the non-linear solver
σ        = 0.5

#-----------------------------------------
# Compute the scalar field over the grid
#-----------------------------------------
# ϕ_ = distribute(Grid, computePatch, computeUboundary, computeVboundary)
# save("output/data/minkowski_2D/sigma_0p5_np_30_nh_4_source/minkowski_2D_sigma_0p5_psi", extract(ϕ_, 1))
# save("output/data/minkowski_2D/sigma_0p5_np_30_nh_4_source/minkowski_2D_sigma_0p5_deltapsi", extract(map(deltapsi,ϕ_), 1))

#-----------------------------------------
# Now test h-p convergence of the error
# in the solution
#-----------------------------------------
np, ep = pconv(2:2:40, 12)
nh, eh = hconv(0:1:9, 4)

#-----------------------------------------
# Save all the data for solution and convergence 
#-----------------------------------------
# save("output/data/minkowski_2D_sigma_0p9_psi", extract(ϕ_, 1))
# save("output/data/minkowski_2D_sigma_0p9_deltapsi_pconv", np, ep)
# save("output/data/minkowski_2D_sigma_0p9_deltapsi_hconv", nh, eh)
