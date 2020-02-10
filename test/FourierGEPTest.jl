#--------------------------------------------------------------------
# DNSS.jl
# Soham 08-2019
# Find the eigen values and the eigen vectors of the Laplace operator
# on a sphere.
#--------------------------------------------------------------------
# Test Fourier series and it's derivatives
# Code up the metric components and do the sums explicitly.

struct Θ end
struct Φ end 
S = FourierGEP{Θ, 40, Float64}(0, 2π)
D = derivative(S) 
θ = Field(S, u->u)

sinθ = Field(S, θ->sin(θ))
cosθ = D*sinθ

using PyPlot
plot(sinθ)
plot(cosθ)
savefig("FourierGEP.pdf")
close()
