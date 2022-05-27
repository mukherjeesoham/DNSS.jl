#--------------------------------------------------------------------
# DNSS.jl
# Soham 03-2022
# Simulate the Minkowski spacetime
#--------------------------------------------------------------------
using NLsolve, Printf

function minkowski(PS::ProductSpace{S1, S2})::NTuple{2, Field{ProductSpace{S1, S2}}} where {S1, S2}
    a = Field(PS, (u,v)->2) 
    η = Field(PS, (u,v)->v-u) 
    return (a, η)
end

function computeUboundary(PS::ProductSpace{S1, S2})::NTuple{2, Field{S2}} where {S1, S2}
    a0 = extractUboundary(Field(PS, (u,v)->2), :incoming) 
    η0 = extractUboundary(Field(PS, (u,v)->(v-u)), :incoming)
    return (a0, η0)
end

function computeVboundary(PS::ProductSpace{S1, S2})::NTuple{2, Field{S1}} where {S1, S2}
    a0 = extractVboundary(Field(PS, (u,v)->2), :incoming) 
    η0 = extractVboundary(Field(PS, (u,v)->(v-u)), :incoming)
    return (a0, η0)
end

function initialguess(PS::ProductSpace{S1, S2})::NTuple{2, Field{ProductSpace{S1, S2}}} where {S1, S2}
    (a0, η0) = minkowski(PS)
    n1       = Field(PS, (u,v)-> 1e-3 * rand())
    n2       = Field(PS, (u,v)-> 1e-3 * rand())
    return (a0 + n1, η0 + n2) 
end

function computePatch(PS::ProductSpace{S1, S2}, ubnd::NTuple{2, Field{S2}}, vbnd::NTuple{2, Field{S1}})::NTuple{2, Field{ProductSpace{S1, S2}}} where {S1, S2}
    B = incomingboundary(PS)
    I = identity(PS)
    DU, DV = derivative(PS)
    ∂U = combineUVboundary(ubnd, vbnd, :incoming)

    function constr(U::NTuple{2, Field{S}})::NTuple{2, Field{S}} where {S}
        (a, η) = U
        C1 = DU*(DU*η) - (2/a)*(DU*a)*(DU*η)
        C2 = DV*(DV*η) - (2/a)*(DV*a)*(DV*η)
        C  = sqrt(C1*C1 + C2*C2)
        return (C, C)
    end

    function offaxis(U::NTuple{2, Field{S}})::NTuple{2, Field{S}} where {S}
        (a, η) = U
        F1 = DU*(DV*a) - (1/a)*(DU*a)*(DV*a) + (a/η)*(DU*(DV*η))
        F2 = DU*(DV*η) + (1/η)*(DU*η)*(DV*η) + (1/4)*(1/η)*(a^2)
        return (F1, F2) 
    end

    function onaxis(U::NTuple{2, Field{S}})::NTuple{2, Field{S}} where {S}
        (a, η) = U
        A1 = DU*(DV*a) - (1/a)*(DU*a)*(DV*a)
        A2 = DU*(DV*η)
        return (A1, A2)
    end

    function boundary(U::NTuple{2, Field{S}})::NTuple{2, Field{S}} where {S}
        return U - ∂U 
    end

    function residual(U::NTuple{2, Field{S}})::NTuple{2, Field{S}} where {S}
        F = offaxis(U) + constr(U)
        enforceregularity!(F, onaxis(U))
        enforcebc!(F, boundary(U)) 
        return F
    end

    function f!(fvec::Array{T,1}, x::Array{T,1}) where {T}
        fvec[:] = reshapeFromTuple(residual(reshapeToTuple(PS, 2, x)))
    end

    r = initialguess(PS)
    s = nlsolve(f!, reshapeFromTuple(r); method=solver, autodiff=:forward, show_trace=debug >= 3, ftol=1e-8, iterations=80)
    k = reshapeToTuple(PS, 2, s.zero)

    if debug >= 2
        println(" [Before solve] norm(constraints) = ", norm(constr(r)))
        println(" Converged?                       = ", converged(s))
        println(" [After solve] norm(constraints)  = ", norm(constr(k)))
    end

    return k
end

function constraints(U::NTuple{2, Field{S}})::NTuple{1, Field{S}} where {S}
    (a, η) = U
    DU, DV = derivative(a.space)
    C1 = DU*(DU*η) - (2/a)*(DU*a)*(DU*η)
    C2 = DV*(DV*η) - (2/a)*(DV*a)*(DV*η)
    C  = sqrt(C1*C1 + C2*C2)
    return (C, ) 
end

function pconv(range::StepRange, nh::Int=1)
    println("Testing p convergence")
    n_ = collect(range)
    l_ = zeros(size(n_))
    for index in CartesianIndices(n_) 
        np = n_[index]
        params = Parameters((np, np), (nh,nh), urange, vrange, 2)
        q = distribute(params, computePatch, computeUboundary, computeVboundary)
        l_[index] = rmse(extract(map(constraints, q), 1))
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
        params = Parameters((np, np), (nh, nh), urange, vrange, 2)
        l_[index] = rmse(extract(map(constraints, distribute(params, computePatch, computeUboundary, computeVboundary)), 1))
        l0 = index[1] > 1 ? l_[index[1] - 1] : 1
        @printf("  n = %3i, rmse[2h] / rmse[h] = %e\n", n_[index], l0 / l_[index])
    end
    return (n_, l_)
end

#-----------------------------------------
# Setup simulation grid parameters, 
#-----------------------------------------
npoints  = (8, 8)
npatches = (4, 4) 
urange   = (0.0, 1.0)  
vrange   = (2.0, 3.0)  
nfields  = 2
params   = Parameters(npoints, npatches, urange, vrange, nfields)
solver   = :trust_region 
debug    = 1

#-----------------------------------------
# Check solution solution
#-----------------------------------------
# NOTE: We found the solution to be okayshily stable
# with the modified regularity conditions and after including the constraints. 
# TODO: Check for convergence
# S = distribute(params, computePatch, computeUboundary, computeVboundary, debug)
# save("output/data/minkowski/constraints/minkowski_guv", extract(S, 1))
# save("output/data/minkowski/constraints/minkowski_grr", extract(S, 2))

#-----------------------------------------
# Now test h-p convergence of the error
#-----------------------------------------
np, ep = pconv(2:2:12, 3)
nh, eh = hconv(0:1:7, 4)
save("output/data/minkowski/minkowski_constraints_pconv", np, ep)
save("output/data/minkowski/minkowski_constraints_hconv", nh, eh)

# Notes
# [Adding noise]
# 1/ Adding noise to guv results in the non-linear solver failing to converge.
# 2/ Adding noise to psi allows the non-linear solver to converge. 
# 3/ However the higher modes saturate around 1e-5, which explains why the constraint
#    violations are so large (i.e., once you take the second derivatives of the basis functions
#    their values at the boundary are orders of magnitude higher than in the interior.)
#    Filtering, at least in the Minkowski case is tricky, since most modes saturate around 1e-4. Filtering 
#    the top-half is already quite aggresive but that doesn't work either. If we set all the modes below
#    say 1e-3 to zero, we get the constraints to the satisfied well. 
