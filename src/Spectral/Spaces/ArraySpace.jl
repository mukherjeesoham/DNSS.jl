#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 08-2018
# Define operations for ND spaces
#--------------------------------------------------------------------

export extract, range, Field, stagger

function extract(AoT::Array{Union{ProductSpace, NTuple{3, Field}}, 2}, ID::Int)::Array{Field, 2} 
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

function Base. minimum(AoF::Array{Field, 2})::Number 
    VoF = collect(Iterators.flatten(value.(AoF)))
    return minimum(skipmissing(VoF))
end

function Base. maximum(AoF::Array{Field, 2})::Number
    VoF = collect(Iterators.flatten(value.(AoF)))
    return maximum(skipmissing(VoF))
end

function Field(AoS::T, map::Function) where {T<:Array{Union{NTuple{3, Field}, ProductSpace}, 2}}
    AoF = Array{Field, 2}(undef, size(AoS))
    for index in CartesianIndices(AoF)
        AoF[index] = Field(AoS[index], map)
    end
    return AoF
end

function stagger(params::Parameters, tol::Number; maxiter=100)::Parameters
    ϵ = 0.0
    for iteration in 1:maxiter
        AoS = setup(params)
        AoF = Field(AoS, (u,v)->abs(v-u))
        if minimum(AoF) > tol
            println("Found ϵ that works! ", ϵ)
            println("Distance from the origin ", minimum(AoF))
            return params
        else
            ϵ = ϵ + tol 
            params.vbounds = params.vbounds .+ ϵ
        end
    end
    println("Could not find ϵ that works")
    return params
end
