#--------------------------------------------------------------------
# DNSS.jl
# Soham 08-2019
# Finalize the basis transformation routines once and for all
#--------------------------------------------------------------------

S = ChebyshevGL{U, 20, Float64}(-2.0, 1.0)
R = ChebyshevGL{V, 40, Float64}(-1.0, 3.0)

u = Field(S, u->exp(-u^2))
# v = basistransform(u)
# w = basistransform(v)
# @test u ≈ w

# u = Field(R, u->exp(u))
# v = basistransform(u)
# w = basistransform(v)
# @test u ≈ w

# PS = ProductSpace(S, R)
# uv = Field(PS, (u,v)->exp(u) + exp(v))
# vw = basistransform(uv)
# wx = basistransform(vw) 
# @test uv ≈ wx


function project(PS::Chebyshev{Tag, N, T}, modal::Array{T, 1}, nodal::Array{T,1})::Array{T, 1} where {Tag, N, T}
    PSC   = ChebyshevGL{Tag,N,T}(range(PS)...)
    for index in CartesianIndices(nodal)
        for ordr in 0:order(PSC) 
            nodal[index] += prefactor(ordr+1, order(PSC))*cheb(ordr, naturalcollocation(PSC, index.I[1]))*modal[ordr+1] 
        end
    end
    return nodal
end

plot(u)
plot(project(u, R))
show()

