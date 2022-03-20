#--------------------------------------------------------------------
# DNSS.jl
# Soham 03-2022
# Simulate a self-gravitating scalar field
#--------------------------------------------------------------------

using NLsolve
using PyPlot

# Compute initial data on u0 and v0
function pulse(u::T , v::T)::T where {T<:Number}
    p = 1e-6
    return p*exp(-u^2/0.1) + p*exp(-v^2/0.1)
end

function computeUboundary_(PS::ProductSpace{S1, S2})::NTuple{3, Field{S2}} where {S1, S2}
    a0 = Field(PS, (u,v)->1) 
    η0 = Field(PS, (u,v)->(v-u)/2) 
    ϕ0 = Field(PS, (u,v)->pulse(u,v))
    return computeUboundary((a0, η0, ϕ0))
end

function computeVboundary_(PS::ProductSpace{S1, S2})::NTuple{3, Field{S1}} where {S1, S2}
    a0 = Field(PS, (u,v)->1) 
    η0 = Field(PS, (u,v)->(v-u)/2) 
    ϕ0 = Field(PS, (u,v)->pulse(u,v))
    return computeVboundary((a0, η0, ϕ0))
end

# Recompute initial data at every patch boundary. 
# function computeUboundary(u::NTuple{3, Field{ProductSpace{S1, S2}}})::NTuple{3, Field{S2}} where {S1, S2}
    # return (extractUboundary(u[1], :incoming), solveη(extractUboundary.(u, :incoming)...), extractUboundary(u[3], :incoming))
# end

# function computeVboundary(u::NTuple{3, Field{ProductSpace{S1, S2}}})::NTuple{3, Field{S1}} where {S1, S2}
    # return (extractVboundary(u[1], :incoming), solveη(extractVboundary.(u, :incoming)...), extractVboundary(u[3], :incoming))
# end


function constraints(a::Field{S}, η::Field{S}, ϕ::Field{S})::NTuple{2, Field{S}} where {S}
    DU, DV = derivative(a.space)
    C1 = DU*(DU*η) - (2/a)*(DU*a)*(DU*η) + (4pi*η)*(DU*ϕ)^2
    C2 = DV*(DV*η) - (2/a)*(DV*a)*(DV*η) + (4pi*η)*(DV*ϕ)^2
    return (C1, C2)
end

function residual(a::Field{S}, η::Field{S}, ϕ::Field{S})::NTuple{3, Field{S}} where {S}
    DU, DV = derivative(a.space)
    F1 = DU*(DV*a) - (1/a)*(DU*a)*(DV*a) + (a/η)*(DU*(DV*η)) + (4pi*a)*(DU*ϕ)*(DV*ϕ)
    F2 = DU*(DV*η) + (1/η)*(DU*η)*(DV*η) + (1/4)*(1/η)*(a^2)
    F3 = DU*(DV*ϕ) + (1/η)*(DU*η)*(DV*ϕ) + (1/η)*(DV*η)*(DU*ϕ)
    return (F1, F2, F3)
end

function solveη(a::Field{S}, η::Field{S}, ϕ::Field{S})::Field{S} where {S<:Space{Tag}} where {Tag}
    D = derivative(a.space)
    I = identity(a.space)
    B = incomingboundary(a.space) + outgoingboundary(a.space)
    A = D*D - (2/a)*(D*a)*D + 4pi*(D*ϕ)^2*I
    return solve(A ⊕ B, B*η)
end

function lineconstraint(a::Field{S}, η::Field{S}, ϕ::Field{S})::Number where {S<:Space{Tag}} where {Tag}
    D = derivative(a.space)
    h = D*(D*η) - (2/a)*(D*a)*(D*η) + (4pi*η)*(D*ϕ)^2
    return L2(h)
end



function initialguess(PS::ProductSpace{S1, S2})::NTuple{3, Field{ProductSpace{S1, S2}}} where {S1, S2}
    a0 = Field(PS, (u,v)->1) 
    η0 = Field(PS, (u,v)->(v-u)/2) 
    ϕ0 = Field(PS, (u,v)->pulse(u,v))
    return (a0, η0, ϕ0)
end

# This gets called multiple times, needs to be efficient.
function computePatch_(PS::ProductSpace{S1, S2}, ubnd::NTuple{3, Field{S2}}, 
                 vbnd::NTuple{3, Field{S1}})::NTuple{3, Field{ProductSpace{S1, S2}}} where {S1, S2}

    function initialguess_(PS::ProductSpace{S1, S2})::NTuple{3, Field{ProductSpace{S1, S2}}} where {S1, S2}
        a0 = Field(PS, (u,v)->1) 
        η0 = Field(PS, (u,v)->(v-u)/2) 
        ϕ0 = Field(PS, (u,v)->pulse(u,v))
        return (a0, η0, ϕ0)
    end

    function residualforsolver_(a::Field{S}, η::Field{S}, ϕ::Field{S}, 
                                abnd::Field{S}, ηbnd::Field{S}, ϕbnd::Field{S},
                                DU::Operator{S}, DV::Operator{S}, B::Operator{S}, I::Operator{S})::NTuple{3, Field{S}} where {S}
        # TODO: Can we improve how we compute the residual?
        F1 = DU*(DV*a) - (1/a)*(DU*a)*(DV*a) + (a/η)*(DU*(DV*η)) + (4pi*a)*(DU*ϕ)*(DV*ϕ)
        F2 = DU*(DV*η) + (1/η)*(DU*η)*(DV*η) + (1/4)*(1/η)*(a^2)
        F3 = DU*(DV*ϕ) + (1/η)*(DU*η)*(DV*ϕ) + (1/η)*(DV*η)*(DU*ϕ)

        # Note that this is slower, but more accurate? 
        # Now enforce boundary conditions
        F1 = (I - B) * F1 + B*abnd 
        F2 = (I - B) * F1 + B*ηbnd 
        F3 = (I - B) * F1 + B*ϕbnd 
        return (F1, F2, F3)
    end

    B = incomingboundary(PS)
    I = identity(PS)
    DU, DV = derivative(PS)
    boundarydata = combineUVboundary(ubnd, vbnd, :incoming)

    # Unwrap this further. Can we write it more transparently?  
    function f!(fvec::Array{T,1}, x::Array{T,1}) where {T}
        interiordata = reshapeToTuple(PS, 3, x)
        # fvec[:] = reshapeFromTuple(mix!(residualforsolver_(interiordata..., DU, DV), B, interiordata .- boundarydata))
        fvec[:] = reshapeFromTuple(residualforsolver_(interiordata..., boundarydata..., DU, DV, B, I))
    end

    if true
        println("\t L2(C1, C2) before solve = ",  L2.(constraints_(initialguess(PS)...)))
    end

    F = nlsolve(f!, reshapeFromTuple(initialguess(PS)); 
                                    method=:trust_region, 
                                    autodiff=:forward,
                                    show_trace=false, 
                                    ftol=1e-11, 
                                    iterations=100)

    Tuples = reshapeToTuple(PS, 3, F.zero)

    if true
        println("Solver converged? ", converged(F))
        println("\t L2(C1, C2) after solve = ",  L2.(constraints_(Tuples...)))
    end

    return Tuples
end


# parameters = Parameters((12, 12), (1, 1), (-2.0, 2.0), (-2.0, 2.0))
parameters = stagger(Parameters((12, 12), (2, 2), (-2.0, 2.0), (-2.0, 2.0), 3), 1e-2)
@show parameters
AoT = distribute(parameters,
                 computePatch_, 
                 computeUboundary_, 
                 computeVboundary_)
                 
@show typeof(AoT)


# AoF = extract(AoT, 1)
# contourf(AoF, 20)
# plotaxis(parameters)
# axis("square")
# tight_layout()
# savefig("output/CollapseA.pdf")
# close()

# AoF = extract(AoT, 2)
# contourf(AoF, 20)
# plotaxis(parameters)
# axis("square")
# tight_layout()
# savefig("output/CollapseR.pdf")
# close()

# AoF = extract(AoT, 3)
# contourf(AoF, 20)
# plotaxis(parameters)
# axis("square")
# tight_layout()
# savefig("output/CollapseF.pdf")
# close()

L2C = maximum([maximum(L2.(constraints_(T...))) for T in AoT])
@show L2C 
