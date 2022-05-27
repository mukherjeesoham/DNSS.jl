#--------------------------------------------------------------------
# DNSS.jl
# Soham 03-2022
# Simulate the collapse of a self-gravitating scalar field
# in Vaidya spacetime
#--------------------------------------------------------------------

using ForwardDiff, QuadGK

function psi(x::T)::T where {T<:Real}
    (a, b) = (-0.5, 0.5)
    # (a, b) = (-1, 1)
    if (a <= x <= b) 
        return x * (1-(x/a))^8  * (1 - (x/b))^8 
    else
        return 0.0
    end
end

function psi_spherical(u::T, v::T)::T where {T<:Number}
    t = v + u
    r = v - u
    α = π/2 - k*t
    if  r == 0.0
        return cos(k*t)
    else
        return cos(k*t) * (sin(k*r) / (k*r))
    end
end

function checkinitialdata(a::Field{S}) where {S}
    PS = derivative(a.space)
    dψ = Field(PS, dpsi)
    r  = Field(PS, grr)
    D  = derivative(PS)
    @show norm(-(2 / a) * (D * a) + 4 * π * r * (dψ)^2)
end

function dpsi(x::T)::T where {T<:Real}
    return ForwardDiff.derivative(psi, x)
end

function guv(v::T)::T where {T<:Real}
    integral, error = quadgk(x->(x * dpsi(x)^2), -1, v, rtol=1e-12)
    return exp(2π * integral)
end

function resolve(U::NTuple{3, Field{S}}, var::Symbol)::NTuple{3, Field{S}} where {S}
    PS = first(U).space
    @show PS
    D = derivative(PS)
    I = identity(PS)
    B = incomingboundary(PS) 

    if var == :a
        (a, η, ψ) = U 
        L = (D*η)*D + (4π * η * (D*ψ)^2) * I 
        lna = linsolve(L + B, B*log(a))
        a   = exp(lna)
        return (a, η, ψ)
    end
end
    
# TODO: Test the simplest initial data for spherically symmetric wave
A  = ChebyshevGL{U, 20, Float64}(0.0, 1.0)
B  = ChebyshevGL{V, 20, Float64}(0.0, 1.0)
PS = ProductSpace(A, B)

k  = 0.01
ψ  = extractUboundary(Field(PS, psi_spherical), :incoming) 
r  = extractUboundary(Field(PS, (u,v)->v-u), :incoming)
a  = extractUboundary(Field(PS, (u,v)->2), :incoming)
@show norm(lineconstraints((a, r, ψ)))
