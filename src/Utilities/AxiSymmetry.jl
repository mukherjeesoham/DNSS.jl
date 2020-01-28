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

