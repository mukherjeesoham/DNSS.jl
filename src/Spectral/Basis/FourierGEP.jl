#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 05-2019
# General Fourier Series: Endpoint Grid
# See Boyd F.2
#--------------------------------------------------------------------

function collocation(S::FourierGEP{Tag, N, T}, i::Int)::T where {Tag, N, T}
    @assert iseven(N) && (i > 0) && (i <= N)
    return (pi/(N/2))*(i-1)
end

function derivative(S::FourierGEP{Tag, N, T}, i::Int, j::Int)::T where {Tag, N, T}
    if i != j
        return ((1/2)*((-1)^(i-j))*cot((1/2)*(collocation(S, i) - collocation(S, j))))
    else
        return 0
    end
end
