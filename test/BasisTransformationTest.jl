#--------------------------------------------------------------------
# DNSS.jl
# Soham 08-2019
# Finalize the basis transformation routines once and for all
#--------------------------------------------------------------------

S = ChebyshevGL{U, 20, Float64}(-3.0, 1.0)
R = ChebyshevGL{V, 40, Float64}(-1.0, 1.0)

PS = ProductSpace(S, R)
uv = Field(PS, (u,v)->exp(u) + exp(v))
vw = basistransform(uv)
wx = basistransform(vw) 
@test uv ≈ wx
