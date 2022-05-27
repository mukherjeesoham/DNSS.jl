#--------------------------------------------------------------------
# DNSS.jl
# Soham 03-2022
# Simulate Schwarzschild spacetime 
#--------------------------------------------------------------------
using NLsolve, Printf, Roots

function find_t_of_UV(U::T, V::T, M::Number)::T where {T<:Number}
    @assert V > 0   # ensure you're in region I or II
    @assert U*V < T(1) # ensure you don't hit the singularity
    if U*V == 0     # r = 2M
        t = V       # enforce uniqueness
    elseif U > 0    # r < 2M
        t = -2M*log(U/V)
    elseif U < 0    # r > 2M
        t = -2M*log(-U/V)
    else
        error("Domain error")
    end
    return t
end

function find_r_of_UV(U::T, V::T, M::Number)::T where {T<:Number}
    @assert V > T(0)       # ensure you're in region I or II
    @assert U*V < T(1)     # ensure you don't hit the singularity
    if U*V == T(0)         # r = 2M
        r = 2M
    else                # r < 2M or r > 2M
        f(r) = (r/2M - 1)*exp(r/2M) + U*V
        try
            # r    = find_zero(f, 2M)
            # NOTE: Introduce bracketing
            if U*V > T(0)
                r = find_zero(f, (eps(T), 2*M))
            else
                # WARNING: Would fail if r > 1000*M
                r = find_zero(f, (2M, 1000*M))
            end
        catch
            @show U, V
            exit()
        end
    end
    @assert r > 0
    return r
end

function find_U_of_tr(t::T, r::T, M::Number)::T where {T<:Number}
    @assert r > 0
    if r == 2M
        U = 0
    else
        rstar = r + 2M*log(abs((r/2M)-1))
        u = t - rstar
        r > 2M ?  U = -exp(-u/4M) : U = +exp(-u/4M)
    end
    return U
end

function find_V_of_tr(t::T, r::T, M::Number)::T where {T<:Number}
    @assert r > 0
    if r == 2M
        V = t           # make the inverse mapping unique.
    else
        rstar = r + 2M*log(abs((r/2M)-1))
        v = t + rstar
        V = exp(v/4M)
    end
    return V
end

function schwarzschild(PS::ProductSpace{S1, S2})::NTuple{2, Field{ProductSpace{S1, S2}}} where {S1, S2}
    η0 = Field(PS, (u,v)->find_r_of_UV(u,v,M))
    ϕ0 = Field(PS, (u,v)->0)
    f0 = ((16*M^3)/η0)*exp(-η0/2M)
    a0 = -sqrt(2*f0) 
    return (a0, η0)
end

function singularity(PS::ProductSpace{S1, S2})::NTuple{2, Field{ProductSpace{S1, S2}}} where {S1, S2}
    a = Field(PS, (u,v)->NaN)
    η = Field(PS, (u,v)->NaN)
    return (a, η)
end

function excision(PS::ProductSpace{S1, S2})::Bool where {S1, S2}
    uv = Field(PS, (u,v)->u*v)
    v  = Field(PS, (u,v)->v)
    return (maximum(uv) >= 1.0) || (minimum(v) <= 0.0) 
end

function excision(bnd::NTuple{2, Field{S}})::Bool where {S}
    return any(map(isnan, bnd))
end

function computeUboundary(PS::ProductSpace{S1, S2})::NTuple{2, Field{S2}} where {S1, S2}
    if excision(PS) == true 
        return map(x->extractUboundary(x, :incoming), singularity(PS))
    else
        return map(x->extractUboundary(x, :incoming), schwarzschild(PS))
    end
end

function computeVboundary(PS::ProductSpace{S1, S2})::NTuple{2, Field{S1}} where {S1, S2}
    if excision(PS) == true 
        return map(x->extractVboundary(x, :incoming), singularity(PS))
    else
        return map(x->extractVboundary(x, :incoming), schwarzschild(PS))
    end
end

function computePatch(PS::ProductSpace{S1, S2}, ubnd::NTuple{2, Field{S2}}, 
                                                vbnd::NTuple{2, Field{S1}})::NTuple{2, Field{ProductSpace{S1, S2}}} where {S1, S2}
    if excision(PS) || excision(ubnd) || excision(vbnd)
        return singularity(PS) 
    else
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
            C1 = DU*(DU*η) - (2/a)*(DU*a)*(DU*η)
            C2 = DV*(DV*η) - (2/a)*(DV*a)*(DV*η)
            return (F1 , F2) 
        end

        function boundary(U::NTuple{2, Field{S}})::NTuple{2, Field{S}} where {S}
            return U - ∂U 
        end

        function residual(U::NTuple{2, Field{S}})::NTuple{2, Field{S}} where {S}
            F = offaxis(U) + λ * constr(U)
            enforcebc!(F, boundary(U)) 
            return F
        end

        function initialguess(PS::ProductSpace{S1, S2})::NTuple{2, Field{ProductSpace{S1, S2}}} where {S1, S2}
            a = Field(PS, (u,v)->2)
            η = Field(PS, (u,v)->(v-u))
            return (a, η)
        end

        function f!(fvec::Array{T,1}, x::Array{T,1}) where {T}
            fvec[:] = reshapeFromTuple(residual(reshapeToTuple(PS, 2, x)))
        end

        s = nlsolve(f!, reshapeFromTuple(initialguess(PS)); method=:trust_region, autodiff=:forward, 
                    show_trace=debug >= 3, ftol=1e-9, iterations=100)
        k = reshapeToTuple(PS, 2, s.zero)

        if debug > 1
            println("    Converged?            = ", converged(s))
        end

        if converged(s) == false
            return singularity(PS)
        else
            return k
        end
    end
end

function pconv(range::StepRange, nh::Int=1)
    println("Testing p convergence")
    n_ = collect(range)
    l_ = zeros(size(n_))
    for index in CartesianIndices(n_) 
        np = n_[index]
        params = Parameters((np, np), (nh,nh), urange, vrange, 2)
        q = distribute(params, computePatch, computeUboundary, computeVboundary)
        l_[index] = rmse(extract(map(constraint, q), 1))
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
        l_[index] = rmse(extract(map(constraint, distribute(params, computePatch, computeUboundary, computeVboundary)), 1))
        l0 = index[1] > 1 ? l_[index[1] - 1] : 1
        @printf("  n = %3i, rmse[2h] / rmse[h] = %e\n", n_[index], l0 / l_[index])
    end
    return (n_, l_)
end

function constraint(U::NTuple{2, Field{S}})::NTuple{1, Field{S}} where {S}
    (a, η) = U
    PS = first(U).space
    DU, DV = derivative(PS)
    C1 = DU*(DU*η) - (2/a)*(DU*a)*(DU*η)
    C2 = DV*(DV*η) - (2/a)*(DV*a)*(DV*η)
    C  = sqrt(C1*C1 + C2*C2)
    return (C,) 
end

function error(U::NTuple{2, Field{S}})::NTuple{2, Field{S}} where {S}
    (a, η) = U
    PS = first(U).space
    if excision(PS) == false
        (as, ηs) = schwarzschild(PS)
    else
        (as, ηs) = singularity(PS)
    end
    return (a - as, η - ηs)
end

#-----------------------------------------
# Setup simulation grid parameters 
#-----------------------------------------
npoints  = (10, 10)
npatches = (12, 12) 
M        = 0.125 
debug    = 1
λ        = 0.0

# urange   = (-3.0, 3.0)
# vrange   = ( 2.0, 8.0)
urange   = (-9.0, -3.0) 
vrange   = ( 2.0,  3.0) 

params   = Parameters(npoints, npatches, urange, vrange, 2)
s        = distribute(params, computePatch, computeUboundary, computeVboundary, debug)
@show minimum(extract(s, 2))
save("output/data/schwarzschild/schwarzschild_M_0p1_a", extract(s, 1))
save("output/data/schwarzschild/schwarzschild_M_0p1_eta", extract(s, 2))
save("output/data/schwarzschild/schwarzschild_M_0p1_constraints", extract(map(constraint, s), 1))

#-----------------------------------------
# Check convergence away from the singularity
# and then near the singularity
# NOTE: We only get convergence when we don't include
# the constraints
#-----------------------------------------
np, ep = pconv(2:2:12, 12)
nh, eh = hconv(1:1:7, 4)
save("output/data/schwarzschild/schwarzschild_M_0p1_pconv", np, ep)
save("output/data/schwarzschild/schwarzschild_M_0p1_hconv", nh, eh)
