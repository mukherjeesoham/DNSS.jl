#--------------------------------------------------------------------
# DNSS.jl 
# Soham 02-2020
# Basis change functions
#--------------------------------------------------------------------

export basistransform, restrict, prolongate

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

function restrict(u::Field{ProductSpace{Chebyshev{Tag1, N1, T},
                                        Chebyshev{Tag2, N2, T}}}, N3::Int, N4::Int)::Field{ProductSpace{Chebyshev{Tag1, N3, T},
                                                                                                        Chebyshev{Tag2, N4, T}}} where {Tag1, Tag2, N1, N2, T}
    @assert N3 <= N1
    @assert N4 <= N2
    return Field(ProductSpace(Chebyshev{Tag1, N3, T}(range(u.space.S1)...),
                              Chebyshev{Tag2, N4, T}(range(u.space.S2)...)), u.value[1:N3, 1:N4])
end

function prolongate(u::Field{ProductSpace{Chebyshev{Tag1, N1, T},
                                          Chebyshev{Tag2, N2, T}}}, N3::Int, N4::Int)::Field{ProductSpace{Chebyshev{Tag1, N3, T},
                                                                                                          Chebyshev{Tag2, N4, T}}} where {Tag1, Tag2, N1, N2, T}
    @assert N3 >= N1
    @assert N4 >= N2
    v = zeros(T, N3, N4)
    v[1:N1, 1:N2] = u.value
    return Field(ProductSpace(Chebyshev{Tag1, N3, T}(range(u.space.S1)...),
                              Chebyshev{Tag2, N4, T}(range(u.space.S2)...)), v)
end

function restrict(u::Field{ProductSpace{ChebyshevGL{Tag1, N1, T},
                                        ChebyshevGL{Tag2, N2, T}}}, N3::Int, N4::Int)::Field{ProductSpace{ChebyshevGL{Tag1, N3, T},
                                                                                                          ChebyshevGL{Tag2, N4, T}}} where {Tag1, Tag2, N1, N2, T}
    return basistransform(restrict(basistransform(u), N3, N4))
end

function prolongate(u::Field{ProductSpace{ChebyshevGL{Tag1, N1, T},
                                          ChebyshevGL{Tag2, N2, T}}}, N3::Int, N4::Int)::Field{ProductSpace{ChebyshevGL{Tag1, N3, T},
                                                                                                            ChebyshevGL{Tag2, N4, T}}} where {Tag1, Tag2, N1, N2, T}
    return basistransform(prolongate(basistransform(u), N3, N4))
end
