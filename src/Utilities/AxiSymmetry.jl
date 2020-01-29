#--------------------------------------------------------------------
# DNSS.jl
# Soham 04-2019
# Define functions relevant for axis-symmetry
#--------------------------------------------------------------------
export mix! 

function mix!(u::Field{S}, A::Operator{S}, v::Field{S})::Field{S} where {S}
    for index in CartesianIndices(u.value)
        if A.value[index.I..., index.I...] == 1
            u.value[index] = v.value[index]
        end
    end
    return u
end

function mix!(u::NTuple{3, Field{S}}, A::Operator{S},v:: NTuple{3, Field{S}})::NTuple{3, Field{S}} where {S}
    for index in CartesianIndices(u[1].value)
        if A.value[index.I..., index.I...] == 1
            u[1].value[index] = v[1].value[index]
            u[2].value[index] = v[2].value[index]
            u[3].value[index] = v[3].value[index]
        end
    end
    return u
end
