#--------------------------------------------------------------------
# DNSS.jl
# Soham 03-2022
# Simulate the Minkowski spacetime
#--------------------------------------------------------------------

using NLsolve, Printf

# TODO: Add some noise to the initial and boundary conditions
# TODO: Instead of addding it to guv, what happens when we add a little
# noise to ψ?
noise = 1e-4 * rand()

function computeUboundary(PS::ProductSpace{S1, S2})::NTuple{3, Field{S2}} where {S1, S2}
    a0 = extractUboundary(Field(PS, (u,v)->2), :incoming) 
    η0 = extractUboundary(Field(PS, (u,v)->(v-u)), :incoming)
    ϕ0 = extractUboundary(Field(PS, (u,v)->noise), :incoming)
    return (a0, η0, ϕ0)
end

function computeVboundary(PS::ProductSpace{S1, S2})::NTuple{3, Field{S1}} where {S1, S2}
    a0 = extractVboundary(Field(PS, (u,v)->2 + noise), :incoming) 
    η0 = extractVboundary(Field(PS, (u,v)->(v-u)), :incoming)
    ϕ0 = extractVboundary(Field(PS, (u,v)->noise), :incoming)
    return (a0, η0, ϕ0)
end

function initialguess(PS::ProductSpace{S1, S2})::NTuple{3, Field{ProductSpace{S1, S2}}} where {S1, S2}
    a0 = Field(PS, (u,v)->2) 
    η0 = Field(PS, (u,v)->v-u) 
    ϕ0 = Field(PS, (u,v)->noise)
    return (a0, η0, ϕ0)
end

function computePatch(PS::ProductSpace{S1, S2}, ubnd::NTuple{3, Field{S2}}, vbnd::NTuple{3, Field{S1}})::NTuple{3, Field{ProductSpace{S1, S2}}} where {S1, S2}
    B = incomingboundary(PS)
    I = identity(PS)
    DU, DV = derivative(PS)
    ∂U = combineUVboundary(ubnd, vbnd, :incoming)
    btol = norm(map(norm ∘ lineconstraints, (ubnd, vbnd)))

    function constr(U::NTuple{3, Field{S}})::NTuple{3, Field{S}} where {S}
        (a, η, ϕ) = U
        C1 = DU*(DU*η) - (2/a)*(DU*a)*(DU*η) + (4pi*η)*(DU*ϕ)^2
        C2 = DV*(DV*η) - (2/a)*(DV*a)*(DV*η) + (4pi*η)*(DV*ϕ)^2
        C  = sqrt(C1*C1 + C2*C2)
        return (C, C, C)
    end

    function offaxis(U::NTuple{3, Field{S}})::NTuple{3, Field{S}} where {S}
        (a, η, ϕ) = U
        F1 = DU*(DV*a) - (1/a)*(DU*a)*(DV*a) + (a/η)*(DU*(DV*η)) + (4pi*a)*(DU*ϕ)*(DV*ϕ)
        F2 = DU*(DV*η) + (1/η)*(DU*η)*(DV*η) + (1/4)*(1/η)*(a^2)
        F3 = DU*(DV*ϕ) + (1/η)*(DU*η)*(DV*ϕ) + (1/η)*(DV*η)*(DU*ϕ)
        return (F1, F2, F3) 
    end

    function onaxis(U::NTuple{3, Field{S}})::NTuple{3, Field{S}} where {S}
        (a, η, ϕ) = U
        A1 = (DU*η) * (DV*η) + (1/4) * a^2
        A2 = DU*η + DV*η
        A3 = DU*ϕ - DV*ϕ
        return (A1, A2, A3)
    end

    function boundary(U::NTuple{3, Field{S}})::NTuple{3, Field{S}} where {S}
        return U - ∂U 
    end

    function residual(U::NTuple{3, Field{S}})::NTuple{3, Field{S}} where {S}
        F = offaxis(U) + λ * constr(U)
        enforceregularity!(F, onaxis(U))
        enforcebc!(F, boundary(U)) 
        return F
    end

    function f!(fvec::Array{T,1}, x::Array{T,1}) where {T}
        fvec[:] = reshapeFromTuple(residual(reshapeToTuple(PS, 3, x)))
    end

    s = nlsolve(f!, reshapeFromTuple(initialguess(PS)); method=solver, autodiff=:forward, show_trace=debug >= 3, ftol=max(stol, btol), iterations=40)
    k = reshapeToTuple(PS, 3, s.zero)

    if debug >= 1
        println("    Converged?            = ", converged(s))
        println("    maximum(Ricci Scalar) = ", maximum(R(k)))
        println("    maximum(Hawking Mass) = ", maximum(M(k)))
        println("    norm(constraints)     = ", norm(constraints(k)))
        println("    minimum(expansion)    = ", minimum(expansion(k)))
    end

    return k
end

function pconv(range::StepRange, nh::Int=1)
    println("Testing p convergence")
    n_ = collect(r)
    l_ = zeros(size(n_))
    for index in CartesianIndices(n_) 
        np = n_[index]
        params = Parameters((np, np), (nh,nh), urange, vrange, 3)
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
        params = Parameters((np, np), (nh, nh), urange, vrange, 3)
        l_[index] = rmse(extract(map(constraints, distribute(params, computePatch, computeUboundary, computeVboundary)), 1))
        l0 = index[1] > 1 ? l_[index[1] - 1] : 1
        @printf("  n = %3i, rmse[2h] / rmse[h] = %e\n", n_[index], l0 / l_[index])
    end
    return (n_, l_)
end

# Implement some local filtering
function aggresivefilter(u::Field{S}) where {S}
    # XXX: Setting an threshold of 1e-5
    threshold = 1e-3
    ulm = basistransform(u)
    for index in CartesianIndices(u.value)
        if ulm.value[index] <= threshold
            ulm.value[index] = 1e-15
        end
    end
    return basistransform(ulm)
end

function aggresivefilter(U::NTuple{N,Field{S}}) where {N,S}
    return map(aggresivefilter, U)
end

#-----------------------------------------
# Setup simulation grid parameters, 
#-----------------------------------------
# FIXME: Does this work for n = 2?
npoints  = (18, 18)
npatches = (1, 1) 
urange   = (0.0, 1.0)  
vrange   = (0.0, 1.0)  
params   = Parameters(npoints, npatches, urange, vrange, 3)
solver   = :trust_region 
λ        = 0.0 
itol     = 1e-9
stol     = 1e-9
debug    = 4
path     = "."

# Check solution solution
S = distribute(params, computePatch, computeUboundary, computeVboundary, debug)
# S = map(filtertophalf, S)
# S = map(aggresivefilter, S)
@show rmse(extract(map(constraints, S), 1))
# contourf(extract(S, 1), 20, "$path/collase-a.pdf")
# contourf(extract(S, 2), 20, "$path/collase-eta.pdf")
# contourf(extract(S, 3), 20, "$path/collase-psi.pdf")
# contourf(extract(map(constraints, S), 1), 20, "$path/collapse-constraints.pdf")
# contourf(extract(map(residuals, S), 1), 20, "$path/collapse-residuals.pdf")

# Check convergence
# (np, lp) = pconv(4:2:16, 1)
# (nh, lh) = hconv(0:1:4, 2)
# plotpconv(np, lp, "$path/collapse-constraints-pconv.pdf")
# plotpconv(nh, lh, "$path/collapse-constraints-hconv.pdf")

# al = extract(S, 1)[1,1]
# ηl = extract(S, 2)[1,1]
# ψl = extract(S, 3)[1,1]
# plotmodes(al, "$path/guv_l.pdf")
# plotmodes(ηl, "$path/grr_l.pdf")
# plotmodes(ψl, "$path/psi_l.pdf")

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
#    TODO: Also understand the cancellations that happen for Minkowski for the case of 4 points but not
#    for the more general Vaidya spacetime. 
