#--------------------------------------------------------------------
# DNSS.jl 
# Soham 02-2020
# Basis change functions
#--------------------------------------------------------------------

export linetransform, basistransform

function cheb(m::Int, x::Number)::Number
    if abs(x) <= 1
        return cos(m*acos(x))
    elseif x >= 1
        return cosh(m*acosh(x))
    else
        return ((-1)^m)*cosh(m*acosh(-x))
    end
end

function prefactor(i::Int, N::Int)::Number
    (i == 1 || i == N+1) ? (return 1/2) : (return 1)
end

function linetransform(PS::Chebyshev{Tag, N, T}, modal::Array{T, 1}, nodal::Array{T,1})::Array{T, 1} where {Tag, N, T}
    PSC   = ChebyshevGL{Tag,N,T}(range(PS)...)
    for index in CartesianIndices(nodal)
        for ordr in 0:order(PSC) 
            nodal[index] += prefactor(ordr+1, order(PSC))*cheb(ordr, naturalcollocation(PSC, index.I[1]))*modal[ordr+1] 
        end
    end
    return nodal
end

function linetransform(PS::ChebyshevGL{Tag, N, T}, nodal::Array{T, 1}, modal::Array{T,1})::Array{T, 1} where {Tag, N, T}
    for ordr in 0:order(PS) 
        for index in CartesianIndices(modal)
            modal[ordr+1] += prefactor(index.I[1], order(PS))*cheb(ordr, naturalcollocation(PS, index.I[1]))*nodal[index] 
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
