#--------------------------------------------------------------------
# DNSS.jl
# Soham 03-2022
# Simulate the collapse of a self-gravitating scalar field
# in Vaidya spacetime
# FIXME: Compute initial data
#--------------------------------------------------------------------

using NLsolve, Printf

function psi(u::T, v::T)::T where {T<:Number}
    t = v + u
    r = v - u
    if r == 0.0 
        return cos(k*t) 
    else
        return cos(k*t) * (sin(k*r) / (k*r))
    end
end

function computeUboundary(PS::ProductSpace{S1, S2})::NTuple{3, Field{S2}} where {S1, S2}
    a0 = extractUboundary(Field(PS, (u,v)->2), :incoming) 
    η0 = extractUboundary(Field(PS, (u,v)->(v-u)), :incoming)
    ϕ0 = extractUboundary(Field(PS, psi), :incoming)
    ηs = computeη((a0, η0, ϕ0))
    # @show norm(lineconstraints((a0, ηs, ϕ0))) 
    # @assert norm(lineconstraints((a0, ηs, ϕ0))) < 1e-9
    return (a0, ηs, ϕ0)
end

function computeVboundary(PS::ProductSpace{S1, S2})::NTuple{3, Field{S1}} where {S1, S2}
    a0 = extractVboundary(Field(PS, (u,v)->2), :incoming) 
    η0 = extractVboundary(Field(PS, (u,v)->(v-u)), :incoming)
    ϕ0 = extractVboundary(Field(PS, psi), :incoming)
    ηs = computeη((a0, η0, ϕ0))
    # @show norm(lineconstraints((a0, ηs, ϕ0))) 
    # @assert norm(lineconstraints((a0, ηs, ϕ0))) < 1e-9
    return (a0, ηs, ϕ0)
end

function initialguess(PS::ProductSpace{S1, S2})::NTuple{3, Field{ProductSpace{S1, S2}}} where {S1, S2}
    a0 = Field(PS, (u,v)->2) 
    η0 = Field(PS, (u,v)->v-u) 
    ϕ0 = Field(PS, psi)
    return (a0, η0, ϕ0)
end

function computePatch(PS::ProductSpace{S1, S2}, ubnd::NTuple{3, Field{S2}}, 
                                                vbnd::NTuple{3, Field{S1}})::NTuple{3, Field{ProductSpace{S1, S2}}} where {S1, S2}
    B = incomingboundary(PS)
    I = identity(PS)
    DU, DV = derivative(PS)
    ∂U = combineUVboundary(ubnd, vbnd, :incoming)

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
        A2 = DU*(DV*η)
        A3 = DU*(DV*ϕ)
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

    s = nlsolve(f!, reshapeFromTuple(initialguess(PS)); method=:trust_region, 
                autodiff=:forward, show_trace=debug >= 3, ftol=1e-9, iterations=40)
    k = reshapeToTuple(PS, 3, s.zero)

    if debug >= 2
        println(" [Before solve] norm(constraints) = ", norm(constraints(k)))
        println("                Converged?        = ", converged(s))
        println(" [After solve]  norm(constraints) = ", norm(constraints(k)))
    end

    return k
end

function pconv(range::StepRange, nh::Int=1)
    println("Testing p convergence")
    n_ = collect(range)
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

function constraint(U::NTuple{3, Field{S}})::NTuple{1, Field{S}} where {S}
    (a, η, ϕ) = U
    PS = first(U).space
    DU, DV = derivative(PS)
    C1 = DU*(DU*η) - (2/a)*(DU*a)*(DU*η) + (4pi*η)*(DU*ϕ)^2
    C2 = DV*(DV*η) - (2/a)*(DV*a)*(DV*η) + (4pi*η)*(DV*ϕ)^2
    C  = sqrt(C1*C1 + C2*C2)
    return (C,) 
end

#-----------------------------------------
# Setup simulation grid parameters, 
#-----------------------------------------
npoints  = (4, 4)
npatches = (4, 4) 
urange   = (2.0, 4.0)  
vrange   = (5.0, 7.0)  
params   = Parameters(npoints, npatches, urange, vrange, 3)
λ        = 0.0 
debug    = 1
k        = 0.1

S = distribute(params, computePatch, computeUboundary, computeVboundary, debug)
save("output/data/scalarwave_SG/scalarwave_SG_psi", extract(S, 3))
save("output/data/scalarwave_SG/scalarwave_SG_constraints", extract(map(constraint, S), 1))

# np, ep = pconv(2:2:14, 2)
nh, eh = hconv(1:1:5, 4)
# save("output/data/scalarwave_SG/scalarwave_SG_k0p7_off_axis_pconv", np, ep)
# save("output/data/scalarwave_SG/scalarwave_SG_k0p7_off_axis_hconv", nh, eh)
