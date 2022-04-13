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

""" Takes a vector and reshapes them into tuples. Inverse of
    reshapeFromTuple.
"""
function reshapeToTuple(space::S, N::Int, x::Array{T,1})::NTuple{N, Field}  where {S, T}
    U = reshape(x, :, N)
    return Tuple(reshape(space, U[:, i]) for i in 1:N)
end

function enforcebc!(F::NTuple{N, Field{S}}, B::NTuple{N, Field{S}}) where {S,N}
    r = Field(first(F).space, (u,v)->v-u) 
    for index_ in 1:N
        F[index_].value[end, :] = B[index_].value[end, :]
        F[index_].value[:, end] = B[index_].value[:, end]
    end
end

function enforceregularity!(F::NTuple{N, Field{S}}, A::NTuple{N, Field{S}}) where {S,N}
    r = Field(first(F).space, (u,v)->v-u) 
    for index in CartesianIndices(r.value)
        if r.value[index] == 0.0
            for index_ in 1:N
                F[index_].value[index] = A[index_].value[index]
            end
        end
    end
end
