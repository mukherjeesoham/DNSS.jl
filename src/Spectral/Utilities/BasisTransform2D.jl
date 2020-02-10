#--------------------------------------------------------------------
# DNSS.jl 
# Soham 02-2020
# Basis change functions
#--------------------------------------------------------------------

export basistransform

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
