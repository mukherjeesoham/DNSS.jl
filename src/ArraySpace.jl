#--------------------------------------------------------------------
# DNSS.jl
# Soham 03-2012
# Define operations for Array spaces
#--------------------------------------------------------------------

"""
    Extract a particular field from the array of Tuples  (AoT)
    Also handles exceptions for when the array index does not contain an array 
    but a product space.
    Input: Array of Tuples (AoT)
    Output: Array of Field (AoF)
"""
function extract(AoT::Array{Union{ProductSpace, NTuple{N, Field}}, 2}, ID::Int)::Array{Field, 2} where {N} 
    AoF = Array{Field, 2}(undef, size(AoT))
    for index in CartesianIndices(AoT)
        try
            AoF[index] = AoT[index][ID]
        catch
            AoF[index] = Field(AoT[index], Array{Missing, 2}(missing, size(AoT[index])))
        end
    end
    return AoF
end

"""
 Find the minimum of a field over all the patches
 Input: Array of Field (AoF)
 Output: minimum (scalar)
"""
function Base. minimum(AoF::Array{Field, 2})::Number 
    VoF = collect(Iterators.flatten(value.(AoF)))
    return minimum(skipmissing(VoF))
end

"""
 Find the maximum of a field over all the patches
 Input: Array of Field (AoF)
 Output: maximum (scalar)
"""
function Base. maximum(AoF::Array{Field, 2})::Number
    VoF = collect(Iterators.flatten(value.(AoF)))
    return maximum(skipmissing(VoF))
end

"""
    Takes in an array of spaces, and returns the value of the function on all of those spaces 
    Input: Array of spaces (AoS)
    Output: Array of fields (AoF) 
"""
function Field(AoS::Array{Union{ProductSpace, NTuple{N,Field}}, 2}, map::Function)  where {N}
    AoF = Array{Field, 2}(undef, size(AoS))
    for index in CartesianIndices(AoF)
        AoF[index] = Field(AoS[index], map)
    end
    return AoF
end

"""
Compute the root-mean-square error over the whole grid. 
For each individual patch use spectral integration, and for the whole
grid use a standard L2 norm.
"""
# function rmse(AoF::Matrix{Field{ProductSpace{S1, S2}}})::Real where {S1, S2}
function rmse(AoF::Matrix{T})::Real  where {T<:Field}
    return norm(norm.(AoF))
end

