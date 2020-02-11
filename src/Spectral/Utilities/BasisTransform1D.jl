#--------------------------------------------------------------------
# DNSS.jl 
# Soham 02-2020
# Basis change functions
#--------------------------------------------------------------------

export linetransform, basistransform, restrict, prolongate

function cheb(m::Int, x::Number)::Number
    if abs(x) <= 1
        return cos(m*acos(x))
    elseif x >= 1
        return cosh(m*acosh(x))
    else
        return ((-1)^m)*cosh(m*acosh(-x))
    end
end

function scaling(i::Int, N::Int)::Number
    (i == 1 || i == N+1) ? (return 1/2) : (return 1)
end

function linetransform(PS::Chebyshev{Tag, N, T}, modal::Array{T, 1}, nodal::Array{T,1})::Array{T, 1} where {Tag, N, T}
    PSC   = ChebyshevGL{Tag,N,T}(range(PS)...)
    for index in CartesianIndices(nodal)
        for ordr in 0:order(PSC) 
            nodal[index] += scaling(ordr+1, order(PSC))*cheb(ordr, naturalcollocation(PSC, index.I[1]))*modal[ordr+1] 
        end
    end
    return nodal
end

function linetransform(PS::ChebyshevGL{Tag, N, T}, nodal::Array{T, 1}, modal::Array{T,1})::Array{T, 1} where {Tag, N, T}
    for ordr in 0:order(PS) 
        for index in CartesianIndices(modal)
            modal[ordr+1] += scaling(index.I[1], order(PS))*cheb(ordr, naturalcollocation(PS, index.I[1]))*nodal[index] 
        end
    end
    return (2/order(PS))*modal
end

function basistransform(u::Field{ChebyshevGL{Tag, N, T}})::Field{Chebyshev{Tag, N, T}} where {Tag, N, T}
    nodal = u.value
    modal = zeros(T, size(nodal))
    return Field(Chebyshev{Tag, N, T}(range(u.space)...), linetransform(u.space, nodal, modal))
end

function basistransform(u::Field{Chebyshev{Tag, N, T}})::Field{ChebyshevGL{Tag, N, T}} where {Tag, N, T}
    modal = u.value
    nodal = zeros(T, size(modal))
    return Field(ChebyshevGL{Tag, N, T}(range(u.space)...), linetransform(u.space, modal, nodal))
end

function restrict(u::Field{Chebyshev{Tag, N1, T}}, N2::Int)::Field{Chebyshev{Tag, N2, T}} where {Tag, N1, T}
    @assert N2 <= N1
    return Field(Chebyshev{Tag, N2, T}(range(u.space)...), u.value[1:N2])
end

function prolongate(u::Field{Chebyshev{Tag, N1, T}}, N2::Int)::Field{Chebyshev{Tag, N2, T}} where {Tag, N1, T}
    @assert N2 >= N1
    v = zeros(T, N2)
    v[1:N1] = u.value
    return Field(Chebyshev{Tag, N2, T}(range(u.space)...), v)
end

function restrict(u::Field{ChebyshevGL{Tag, N1, T}}, N2::Int)::Field{ChebyshevGL{Tag, N2, T}} where {Tag, N1, T}
    return basistransform(restrict(basistransform(u), N2))
end

function prolongate(u::Field{ChebyshevGL{Tag, N1, T}}, N2::Int)::Field{ChebyshevGL{Tag, N2, T}} where {Tag, N1, T}
    return basistransform(prolongate(basistransform(u), N2))
end

