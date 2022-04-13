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

function dpsi(x::T)::T where {T<:Real}
    return ForwardDiff.derivative(psi, x)
end

function guv(v::T)::T where {T<:Real}
    integral, error = quadgk(x->(x * dpsi(x)^2), -1, v, rtol=1e-12)
    return exp(2π * integral)
end
    
A  = ChebyshevGL{U, 200, Float64}(-1.0, 1.0)
r  = Field(A, u->u)
ψ  = Field(A, psi)
dψ = Field(A, dpsi)
a  = Field(A, guv)
@show norm(lineconstraints((a, r, ψ)))

# Check derivatives and their convergence
# TODO: You need about 80 points for the errors to be comparable! 
# Also, setting (a, b) inside the domain wrecks havoc! We loose four orders of accuracy
# Do I really understand spectral methods? 
# If we filter the top half modes, we loose 4 orders of accuracy as well. 
D  = derivative(A)
P1 = D * (D * r)
P2 = -(2 / a) * (D * a) * (D * r)
P3 = 4 * π * r * (D * ψ)^2
P2A = -(2 / a) * (D * a) 
P3A = 4 * π * r * (dψ)^2

@show norm(D * (D * r))
@show norm(-(2 / a) * (D * a) * (D * r) + 4 * π * r * (D * ψ)^2)
@show norm(D * (D * r) -(2 / a) * (D * a) * (D * r) + 4 * π * r * (D * ψ)^2)

P2A = -(2 / a) * (D * a) 
P3A = 4 * π * r * (dψ)^2
@show norm(-(2 / a) * (D * a) + 4 * π * r * (dψ)^2)
plot(P2 + P3A, "./output/computea_error.pdf")

plot(a, "./output/guv.pdf")
plot(r, "./output/grr.pdf")
plot(ψ, "./output/psi.pdf")

# TODO: Check if this works piecewise
B  = ChebyshevGL{U, 200, Float64}(-1.0, 0.0)
C  = ChebyshevGL{U, 200, Float64}( 0.0, 1.0)
ab  = Field(B, guv)
ac  = Field(C, guv)
x   = Field(A, x->x)
xb  = Field(B, x->x)
xc  = Field(C, x->x)


function checkinitialdata(a::Field{S}) where {S}
    PS = derivative(a.space)
    dψ = Field(PS, dpsi)
    r  = Field(PS, grr)
    D  = derivative(PS)
    @show norm(-(2 / a) * (D * a) + 4 * π * r * (dψ)^2)
end

checkinitialdata(a)
checkinitialdata(ab)
checkinitialdata(ac)

# plot(x.value, a.value)
plot(xb.value, ab.value)
plot(xc.value, ac.value)
plt.savefig("./output/piecewise.pdf")
