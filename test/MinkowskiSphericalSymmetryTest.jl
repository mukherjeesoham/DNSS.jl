#--------------------------------------------------------------------
# DNSS.jl
# Soham 03-2022
# Simulate the collapse of a self-gravitating scalar field
#--------------------------------------------------------------------
# TODO: Investigating stability

using NLsolve, Printf

function computeUboundary(PS::ProductSpace{S1, S2})::NTuple{3, Field{S2}} where {S1, S2}
    a0 = extractUboundary(Field(PS, (u,v)->2),   :incoming)
    η0 = extractUboundary(Field(PS, (u,v)->v-u), :incoming) 
    ϕ0 = extractUboundary(Field(PS, (u,v)-> 1e-8 * rand()),   :incoming)
    return (a0, η0, ϕ0)
end

function computeVboundary(PS::ProductSpace{S1, S2})::NTuple{3, Field{S1}} where {S1, S2}
    a0 = extractVboundary(Field(PS, (u,v)->2),   :incoming)
    η0 = extractVboundary(Field(PS, (u,v)->v-u), :incoming) 
    ϕ0 = extractVboundary(Field(PS, (u,v)->0)  , :incoming)
    return (a0, η0, ϕ0)
end

#----------------------------------------------
# Compute function on each patch 
#----------------------------------------------
function computePatch(PS::ProductSpace{S1, S2}, ubnd::NTuple{3, Field{S2}}, vbnd::NTuple{3, Field{S1}})::NTuple{3, Field{ProductSpace{S1, S2}}} where {S1, S2}
    B = incomingboundary(PS)
    I = identity(PS)
    DU, DV = derivative(PS)

    function enforceregularity!(F::NTuple{3, Field{S}}, A::NTuple{3, Field{S}})::NTuple{3, Field{S}} where {S}
	    r = Field(first(F).space, (u,v)->v-u) 
        for index in CartesianIndices(r.value)
            if r.value[index] == 0.0
                F[1].value[index] = A[1].value[index]
                F[2].value[index] = A[2].value[index]
                F[3].value[index] = A[3].value[index]
            end
        end
        return F 
    end

    function offaxis(U::NTuple{3, Field{S}})::NTuple{3, Field{S}} where {S}
        (a, η, ϕ) = U
        F1 = DU*(DV*a) - (1/a)*(DU*a)*(DV*a) + (a/η)*(DU*(DV*η)) + (4pi*a)*(DU*ϕ)*(DV*ϕ)
        F2 = DU*(DV*η) + (1/η)*(DU*η)*(DV*η) + (1/4)*(1/η)*(a^2)
        F3 = DU*(DV*ϕ) + (1/η)*(DU*η)*(DV*ϕ) + (1/η)*(DV*η)*(DU*ϕ)
        C1 = DU*(DU*η) - (2/a)*(DU*a)*(DU*η) + (4pi*η)*(DU*ϕ)^2
        C2 = DV*(DV*η) - (2/a)*(DV*a)*(DV*η) + (4pi*η)*(DV*ϕ)^2
        C  = sqrt(C1*C1 + C2*C2)
        return λa .* (F1 + C, F2 + C, F3 + C)
    end

    function onaxis(U::NTuple{3, Field{S}})::NTuple{3, Field{S}} where {S}
        (a, η, ϕ) = U
        A1 = (DU*η) * (DV*η) + (1/4) * a^2
        A2 = DU*η + DV*η
        A3 = DU*ϕ - DV*ϕ
        return λb .* (A1, A2, A3)
    end

    function boundary(U::NTuple{3, Field{S}})::NTuple{3, Field{S}} where {S}
        U = (B * U[1], B * U[2], B * U[3])
        return λc .* (U .- combineUVboundary(ubnd, vbnd, :incoming))
    end

    function residual(U::NTuple{3, Field{S}})::NTuple{3, Field{S}} where {S}
        F = offaxis(U)
        A = onaxis(U)
        R = boundary(U)
        return enforcebc!(enforceregularity!(F, A), B, R) 
    end

    function initialguess(PS::ProductSpace{S1, S2})::NTuple{3, Field{ProductSpace{S1, S2}}} where {S1, S2}
        a0 = Field(PS, (u,v)->2) 
        η0 = Field(PS, (u,v)->v-u) 
        ϕ0 = Field(PS, (u,v)->0)
        return (a0, η0, ϕ0)
    end

    function f!(fvec::Array{T,1}, x::Array{T,1}) where {T}
        fvec[:] = reshapeFromTuple(residual(reshapeToTuple(PS, 3, x)))
    end

    U0 = initialguess(PS)

    if debug >= 1
        println("    Before solve            ")
        println("    maximum(Constraints)  = ", norm(constraints(U0)))
    end

    sol = nlsolve(f!, reshapeFromTuple(U0); method=solver, autodiff=:forward, show_trace=true, extended_trace=false, ftol=1e-9, iterations=100)
    ToF = reshapeToTuple(PS, 3, sol.zero)

    if debug >= 1
        println("    Converged?            = ", converged(sol))
        println("    maximum(Constraints)  = ", norm(constraints(ToF)))
    end

    return ToF
end

#----------------------------------------------
# Constraint equations and convergence 
#----------------------------------------------
function constraints(U::NTuple{3, Field})::NTuple{1, Field}
    PS = prolongate(first(U).space)
    DU, DV = derivative(PS)
    a = project(U[1], PS) 
    η = project(U[2], PS) 
    ϕ = project(U[3], PS) 
    C1 = DU*(DU*η) - (2/a)*(DU*a)*(DU*η) + (4pi*η)*(DU*ϕ)^2
    C2 = DV*(DV*η) - (2/a)*(DV*a)*(DV*η) + (4pi*η)*(DV*ϕ)^2
    return (sqrt(C1*C1 + C2*C2), )
end

function logconstraints(U::NTuple{3, Field})::NTuple{1, Field}
    return log.(constraints(U))
end

function pconv(min, max)
    println("Testing p convergence")
    n_ = collect(min:2:max)
    l_ = zeros(size(n_))
    for index in CartesianIndices(n_) 
        np = n_[index]
        nh = 1
        params = Parameters((np, np), (nh,nh), urange, vrange, nfields)
        l_[index] = rmse(extract(map(constraints, distribute(params, computePatch, computeUboundary, computeVboundary)), 1))
        @printf("  n = %i, rmse = %e\n", n_[index], l_[index])
    end
    plotpconv(n_, l_, "$path/collapse-constraints-pconv.pdf")
end

function hconv(min, max)
    println("Testing h convergence")
    n_ = collect(min:max)
    l_ = zeros(size(n_))
    for index in CartesianIndices(n_) 
        nh = n_[index] = 2^n_[index]
        np = 4
        params = Parameters((np, np), (nh, nh), urange, vrange, nfields)
        l_[index] = rmse(extract(map(constraints, distribute(params, computePatch, computeUboundary, computeVboundary)), 1))
        l0 = index[1] > 1 ? l_[index[1] - 1] : 1
        @printf("  n = %3i, rmse[2h] / rmse[h] = %e\n", n_[index], l0 / l_[index])
    end
    plothconv(n_, l_, "$path/collapse-constraints-hconv.pdf")
end

#-----------------------------------------
# Setup simulation grid parameters, 
# call distribute and plot solutions
#-----------------------------------------
npoints  = (18, 18)
npatches = (1, 1) 
urange   = (0.0, 2.0) #(1.0, 3.0)
vrange   = (0.0, 2.0) #(5.0, 7.0)
nfields  = 3
nprocs   = 3
solver   = :trust_region # :trust_region :Newton :Anderson
debug    = 2
idtol    = 1e-5
exclude  = true
params   = Parameters(npoints, npatches, urange, vrange, nfields)
path     = "."
λa       = 1.0
λb       = 1.0 
λc       = 1.0
λ        = 1.0

U = distribute(params, computePatch, computeUboundary, computeVboundary, debug)
@show rmse(extract(map(constraints, U), 1))
# contourf(extract(U, 3), 20, "$path/collase-psi.pdf")
# contourf(extract(map(logconstraints, U), 1), 20, "$path/collapse-constraints.pdf")

#-----------------------------------------
# Check and plot convergence
#-----------------------------------------
# pconv(2, 22)
# hconv(0, 6)

