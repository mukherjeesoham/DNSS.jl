#--------------------------------------------------------------------
# DNSS.jl 
# Soham 03-2022
# Spectral utitlies for changing the number of basis 
# functions in each element. 
#--------------------------------------------------------------------

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

function embed!(u::Vector{T}, v::Vector{T})::Vector{T} where {T}
    a = length(u)
    b = length(v)
    for i in 1:min(a,b)
        v[i] = u[i]
    end
    return v
end

function embed!(u::Matrix{T}, v::Matrix{T})::Matrix{T} where {T}
    (a, b) = size(u)
    (c, d) = size(v) 
    for i in 1:min(a,c), j in 1:min(b,d)
        v[i,j] = u[i,j]
    end
    return v
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

function project(u::Field{Chebyshev{Tag, N1, T}}, S::Chebyshev{Tag, N2, T})::Field{Chebyshev{Tag, N2, T}} where {Tag, N1, N2, T}
    v = zeros(T, N2)
    return Field(Chebyshev{Tag, N2, T}(range(u.space)...), embed!(u.value, v))
end

function project(u::Field{ChebyshevGL{Tag, N1, T}}, S::ChebyshevGL{Tag, N2, T})::Field{ChebyshevGL{Tag, N2, T}} where {Tag, N1, N2, T}
    return basistransform(project(basistransform(u), basistransform(S)))
end


function Base. filter!(u::Field{ChebyshevGL{Tag, N, T}})::Field{ChebyshevGL{Tag, N, T}} where {Tag, N, T}
    SR = Chebyshev{Tag, Int(ceil(N * 2/3)), T}(range(u.space)...)
    S  = Chebyshev{Tag, N, T}(range(u.space)...)
    return basistransform(project(project(basistransform(u), SR), S))
end

function basistransform(u::Field{S}) where {S<:ProductSpace{ChebyshevGL{Tag1, N1, T}, 
                                                            ChebyshevGL{Tag2, N2, T}}} where {Tag1, Tag2, N1, N2, T}
    nodal   = u.value
    modalS1 = zeros(T, size(nodal))
    modalS2 = zeros(T, size(nodal))

    for index in 1:size(u.space.S2)
        modalS1[:, index] = linetransform(u.space.S1, nodal[:, index], modalS1[:, index])
    end
    
    for index in 1:size(u.space.S1)
        modalS2[index, :] = linetransform(u.space.S2, modalS1[index, :], modalS2[index, :])
    end

    return Field(ProductSpace(Chebyshev{Tag1, N1, T}(u.space.S1.min, u.space.S1.max), 
                              Chebyshev{Tag2, N2, T}(u.space.S2.min, u.space.S2.max)), modalS2)
end

function basistransform(u::Field{S}) where {S<:ProductSpace{Chebyshev{Tag1, N1, T}, 
                                                            Chebyshev{Tag2, N2, T}}} where {Tag1, Tag2, N1, N2, T}
    modal   = u.value
    nodalS1 = zeros(T, size(modal))
    nodalS2 = zeros(T, size(modal))

    for index in 1:size(u.space.S1)
        nodalS1[index, :] = linetransform(u.space.S2, modal[index, :], nodalS1[index, :])
    end

    for index in 1:size(u.space.S2)
        nodalS2[:, index] = linetransform(u.space.S1, nodalS1[:, index], nodalS2[:, index])
    end

    return Field(ProductSpace(ChebyshevGL{Tag1, N1, T}(u.space.S1.min, u.space.S1.max), 
                              ChebyshevGL{Tag2, N2, T}(u.space.S2.min, u.space.S2.max)), nodalS2)
end


function project(u::Field{ProductSpace{Chebyshev{Tag1, N1, T},
                                       Chebyshev{Tag2, N2, T}}}, PS::ProductSpace{Chebyshev{Tag1, N3, T}, 
                                                                                  Chebyshev{Tag2, N4, T}})::Field{ProductSpace{Chebyshev{Tag1, N3, T},
                                                                                                                               Chebyshev{Tag2, N4, T}}} where {Tag1, Tag2, N1, N2, N3, N4, T}
    @assert range(u.space) == range(PS) 
    v = zeros(T, N3, N4)
    return Field(ProductSpace(Chebyshev{Tag1, N3, T}(range(u.space.S1)...),
                              Chebyshev{Tag2, N4, T}(range(u.space.S2)...)), embed!(u.value, v))
end


function project(u::Field{ProductSpace{ChebyshevGL{Tag1, N1, T},
                                       ChebyshevGL{Tag2, N2, T}}}, PS::ProductSpace{ChebyshevGL{Tag1, N3, T},
                                                                                    ChebyshevGL{Tag2, N4, T}})::Field{ProductSpace{ChebyshevGL{Tag1, N3, T},
                                                                                                                                   ChebyshevGL{Tag2, N4, T}}} where {Tag1, Tag2, N1, N2, N3, N4, T}
    return basistransform(project(basistransform(u), basistransform(PS)))
end

function basistransform(PS::ProductSpace{ChebyshevGL{Tag1, N1, T}, ChebyshevGL{Tag2, N2, T}}) where {Tag1, Tag2, N1, N2, T} 
    A = Chebyshev{Tag1, N1, T}(range(PS.S1)...) 
    B = Chebyshev{Tag2, N2, T}(range(PS.S2)...) 
    return ProductSpace(A,B)
end

function basistransform(PS::ProductSpace{Chebyshev{Tag1, N1, T}, Chebyshev{Tag2, N2, T}}) where {Tag1, Tag2, N1, N2, T} 
    A = ChebyshevGL{Tag1, N1, T}(range(PS.S1)...) 
    B = ChebyshevGL{Tag2, N2, T}(range(PS.S2)...) 
    return ProductSpace(A,B)
end

function basistransform(S::ChebyshevGL{Tag, N, T}) where {Tag, N, T} 
    return  Chebyshev{Tag, N, T}(range(S)...) 
end

function basistransform(S::Chebyshev{Tag, N, T}) where {Tag, N, T} 
    return  ChebyshevGL{Tag, N, T}(range(S)...) 
end

function prolongate(S::ChebyshevGL{Tag, N, T}) where {Tag, N, T} 
    return ChebyshevGL{Tag, 2*N, T}(range(S)...)
end

function restrict(S::ChebyshevGL{Tag, N, T}) where {Tag, N, T} 
    N2 = iseven(N) ? div(N, 2) : div(N + 1, 2) 
    return ChebyshevGL{Tag, N2, T}(range(S)...) 
end

function prolongate(PS::ProductSpace{ChebyshevGL{Tag1, N1, T}, ChebyshevGL{Tag2, N2, T}}) where {Tag1, Tag2, N1, N2, T} 
    A = ChebyshevGL{Tag1, 2*N1, T}(range(PS.S1)...) 
    B = ChebyshevGL{Tag2, 2*N2, T}(range(PS.S2)...) 
    return ProductSpace(A,B)
end

function restrict(PS::ProductSpace{ChebyshevGL{Tag1, N1, T}, ChebyshevGL{Tag2, N2, T}}) where {Tag1, Tag2, N1, N2, T} 
    N3 = iseven(N1) ? div(N1, 2) : div(N1 + 1, 2) 
    N4 = iseven(N1) ? div(N2, 2) : div(N2 + 1, 2) 
    A = ChebyshevGL{Tag1, N3, T}(range(PS.S1)...) 
    B = ChebyshevGL{Tag2, N4, T}(range(PS.S2)...) 
    return ProductSpace(A,B)
end

function restrict(u::Field{S}) where {S}
    return project(u, restrict(u.space))
end

function prolongate(u::Field{S}) where {S}
    return project(u, prolongate(u.space))
end

function restrict(U::NTuple{N,Field{S}}) where {N,S}
    return map(restrict, U)
end

function prolongate(U::NTuple{N,Field{S}}) where {N,S}
    return map(prolongate, U)
end

function filtertophalf(u::Field{S}) where {S}
    return (prolongate ∘ restrict)(u)
end

function filtertophalf(U::NTuple{N,Field{S}}) where {N,S}
    return map(filtertophalf, U)
end

