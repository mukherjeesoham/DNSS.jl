#--------------------------------------------------------------------
# DNSS.jl
# Soham 03-2022
# Simulate the collapse of a self-gravitating scalar field
# in Vaidya spacetime
#--------------------------------------------------------------------

using NLsolve, Printf

function phi(x::T)::T where {T<:Number}
    (ul, ur) = (-0.5, 0.5)
    if abs(x) <= 0.5 
        return x * (1-(x/ur))^4  * (1 - (x/ul))^4 
    else
        return 0.0
    end
end

function phi(u::T , v::T)::T where {T<:Number}
    A = 1    
    return A * (phi(u) + phi(v))
end

# FIXME: Are the boundary computations working as expected. i.e. are they converging?
function computeUboundary(PS::ProductSpace{S1, S2})::NTuple{3, Field{S2}} where {S1, S2}
    a0 = extractUboundary(Field(PS, (u,v)->2), :incoming) 
    η0 = extractUboundary(Field(PS, (u,v)->(v-u)), :incoming)
    ϕ0 = extractUboundary(Field(PS, (u,v)->phi(u,v)), :incoming)
    return (a0, computeη((a0, η0, ϕ0)), ϕ0)
end

function computeVboundary(PS::ProductSpace{S1, S2})::NTuple{3, Field{S1}} where {S1, S2}
    a0 = extractVboundary(Field(PS, (u,v)->2), :incoming) 
    η0 = extractVboundary(Field(PS, (u,v)->(v-u)), :incoming)
    ϕ0 = extractVboundary(Field(PS, (u,v)->phi(u,v)), :incoming)
    return (a0, computeη((a0, η0, ϕ0)), ϕ0)
end

function initialguess(PS::ProductSpace{S1, S2})::NTuple{3, Field{ProductSpace{S1, S2}}} where {S1, S2}
    a0 = Field(PS, (u,v)->2) 
    η0 = Field(PS, (u,v)->v-u) 
    ϕ0 = Field(PS, phi)
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

    println("    norm(constraints)     = ", norm(constraints(initialguess(PS))))
    println("    norm(residuals)       = ", norm(residuals(initialguess(PS))))

    btol = 1e-9
    s = nlsolve(f!, reshapeFromTuple(initialguess(PS)); method=solver, autodiff=:forward, show_trace=debug >= 3, ftol=max(stol, btol), iterations=40)
    k = reshapeToTuple(PS, 3, s.zero)

    if debug >= 1
        println("    Converged?            = ", converged(s))
        println("    maximum(Ricci Scalar) = ", maximum(R(k)))
        println("    maximum(Hawking Mass) = ", maximum(M(k)))
        println("    norm(constraints)     = ", norm(constraints(k)))
        println("    norm(residuals)       = ", norm(residuals(k)))
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

#-----------------------------------------
# Setup simulation grid parameters, 
#-----------------------------------------
npoints  = (18, 18)
npatches = (1, 1) 
urange   = (-1.0, 1.0) #(0.0, 1.0)  
vrange   = (-1.0, 1.0) #(0.0, 1.0)  
params   = Parameters(npoints, npatches, urange, vrange, 3)
solver   = :trust_region 
λ        = 0.0 
itol     = 1e-9
stol     = 1e-9
debug    = 4
path     = "."

S = distribute(params, computePatch, computeUboundary, computeVboundary, debug)
# @show rmse(extract(map(constraints, S), 1))
# @show rmse(extract(map(residuals, S), 1))

# contourf(extract(S, 1), 20, "$path/collase-a.pdf")
# contourf(extract(S, 2), 20, "$path/collase-eta.pdf")
# contourf(extract(S, 3), 20, "$path/collase-psi.pdf")
# contourf(extract(map(constraints, S), 1), 20, "$path/collapse-constraints.pdf")
# contourf(extract(map(residuals, S), 1), 20, "$path/collapse-residuals.pdf")

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
