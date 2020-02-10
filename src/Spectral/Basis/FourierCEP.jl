#--------------------------------------------------------------------
# DNSS.jl
# Soham 02-2020
# Cosine Fourier Series: Endpoint Grid
# See Boyd F.3
#--------------------------------------------------------------------

function collocation(S::FourierCEP{Tag, N, T}, i::Int)::T where {Tag, N, T}
    @assert (i > 0) && (i <= N)
    return (pi/(N-1))*(i-1)
end

function derivative(S::FourierCEP{Tag, N, T}, i::Int, j::Int)::T where {Tag, N, T}
    @assert (i > 0 && i <= N) && (j > 0 && j <=N)
    if j == 1 || j == N 
        return derivative(FourierGEP{Tag, N, T}(range(S)...), i, j) 
    else
        return derivative(FourierGEP{Tag, N, T}(range(S)...), i, j)  + derivative(FourierGEP{Tag, N, T}(range(S)...), i, 2N-j) 
    end
end
