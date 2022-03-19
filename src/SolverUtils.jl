#--------------------------------------------------------------------
# DNSS.jl
# Soham 03-2022
# Utilities used in the solver
#--------------------------------------------------------------------

"""
    Takes a set of tuples, and flattens them into
    a vector.
"""
function reshapeFromTuple(U::NTuple{N, Field}) where {N}
    return vcat(reshape.(U)...)
end

"""
    Takes a vector and reshapes them into tuples. Inverse of
    reshapeFromTuple.
"""
function reshapeToTuple(space::S, N::Int, x::Array{T,1})::NTuple{N, Field}  where {S, T}
    U = reshape(x, :, N)
    return Tuple(reshape(space, U[:, i]) for i in 1:N)
end

"""
    Takes a field, and replaces it's boundary values with the 
    residual.
"""
function enforcebc!(u::Field{S}, A::Operator{S}, bndres::Field{S})::Field{S} where {S}
    for index in CartesianIndices(u.value)
        if A.value[index.I..., index.I...] == 1
            u.value[index] = bndres.value[index]
        end
    end
    return u
end

"""
    Takes a tuple of fields, and replaces it's boundary values with the 
    residual.
"""
function enforcebc!(u::NTuple{N, Field{S}}, A::Operator{S},v:: NTuple{N, Field{S}})::NTuple{3, Field{S}} where {S, N}
    for index in CartesianIndices(u[1].value)
        if A.value[index.I..., index.I...] == 1
            for i in 1:N
                u[i].value[index] = v[i].value[index]
            end
        end
    end
    return u
end
