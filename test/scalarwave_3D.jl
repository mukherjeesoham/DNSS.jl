#--------------------------------------------------------------------
# DNSS.jl
# Soham 03-2022
# Simulate a scalar field in fixed Minkowski spacetime
#--------------------------------------------------------------------
using NLsolve, Printf, DoubleFloats

function psi(u::Real, v::Real)::Real
    t = v + u
    r = v - u
    if r == 0
        return cos(k * t) # Explicitly regularized for r = 0 
    else
        return cos(k * t) * (sin(k* r)  / (k * r))
    end
end

function computeUboundary(PS::ProductSpace{S1, S2})::NTuple{1, Field{S2}} where {S1, S2}
    ϕ0 = extractUboundary(Field(PS, psi), :incoming)
    return (ϕ0, )
end

function computeVboundary(PS::ProductSpace{S1, S2})::NTuple{1, Field{S1}} where {S1, S2}
    ϕ0 = extractVboundary(Field(PS, psi), :incoming)
    return (ϕ0, )
end

function initialguess(PS::ProductSpace{S1, S2})::NTuple{1, Field{ProductSpace{S1, S2}}} where {S1, S2}
    ϕ0 = Field(PS, (u,v)->1)
    return (ϕ0, )
end

function computePatch(PS::ProductSpace{S1, S2}, ubnd::NTuple{1, Field{S2}}, vbnd::NTuple{1, Field{S1}})::NTuple{1, Field{ProductSpace{S1, S2}}} where {S1, S2}
    B = incomingboundary(PS)
    I = identity(PS)
    DU, DV = derivative(PS)
    ∂U = combineUVboundary(ubnd, vbnd, :incoming)

    function offaxis(U::NTuple{1, Field{S}})::NTuple{1, Field{S}} where {S}
        (ϕ, ) = U
        η  = Field(PS, (u,v)->v-u)
        a  = Field(PS, (u,v)->2)
        F3 = DU*(DV*ϕ) + (1/η)*(DU*η)*(DV*ϕ) + (1/η)*(DV*η)*(DU*ϕ)
        return (F3, ) 
    end

    function onaxis(U::NTuple{1, Field{S}})::NTuple{1, Field{S}} where {S}
        (ϕ, ) = U
        A3 = DU*(DV*ϕ)
        return (A3, )
    end

    function boundary(U::NTuple{1, Field{S}})::NTuple{1, Field{S}} where {S}
        return U - ∂U 
    end

    function residual(U::NTuple{1, Field{S}})::NTuple{1, Field{S}} where {S}
        F = offaxis(U)
        enforceregularity!(F, onaxis(U))
        enforcebc!(F, boundary(U)) 
        return F
    end

    function f!(fvec::Array{T,1}, x::Array{T,1}) where {T}
        fvec[:] = reshapeFromTuple(residual(reshapeToTuple(PS, 1, x)))
    end

    s = nlsolve(f!, reshapeFromTuple(initialguess(PS)); method=solver, autodiff=:forward, 
                show_trace=debug >= 3, ftol=1e-12, iterations=100)
    k = reshapeToTuple(PS, 1, s.zero)

    if debug >= 2
        println(" Converged? = ", converged(s))
    end
    return k
end

function deltapsi(U::NTuple{1, Field})::NTuple{1, Field} where {S}
    psinum   = project(first(U), prolongate(first(U).space))
    psiexact = Field(psinum.space, psi)
    return (psinum - psiexact, )
end

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
# Setup simulation grid parameters, 
#-----------------------------------------
npoints  = (18, 18)
npatches = (1, 1) 
urange   = (0.0, 4.0) 
vrange   = (0.0, 4.0)  
params   = Parameters(npoints, npatches, urange, vrange, 1)
solver   = :trust_region 
debug    = 3
k        = 0.5

# Check the solution
# It is unstable, and the instability is most promiment near the axis and grows with time. 
# However, these are with partiy conditions on the axis. If we impose the other condition on the axis, the solution
# looks more stable.
# TODO: Check mode fall-off?
S = distribute(params, computePatch, computeUboundary, computeVboundary, debug)
save("output/data/scalarwave_3D/k_04_np_18_nh_1/scalarwave_3D_psi", extract(S, 1))
save("output/data/scalarwave_3D/k_04_np_18_nh_1/scalarwave_3D_deltapsi", extract(map(deltapsi, S), 1))

#-----------------------------------------
# Now test h-p convergence of the error
#-----------------------------------------
# np, ep = pconv(2:2:26, 12)
# nh, eh = hconv(0:1:9, 4)

#-----------------------------------------
# Save all the data for solution and convergence 
#-----------------------------------------
# save("output/data/scalar_wave_spherical_symmetry/scalar_waave_spherical_symmetry_off_axis_pconv", np, ep)
# save("output/data/scalar_wave_spherical_symmetry/scalar_waave_spherical_symmetry_off_axis_hconv", nh, eh)


