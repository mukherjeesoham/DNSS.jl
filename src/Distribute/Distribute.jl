#--------------------------------------------------------------------
# DNSS.jl
# Soham 09-2019
# Functions to distribute over multiple patches
# TODO: Implement singularity excision using struct Singularity.
#--------------------------------------------------------------------

export Parameters, distribute, setup

function setup(params::Parameters{T})::Array{Union{ProductSpace, NTuple{3, Field}}} where {T}
    ustops = range(params.ubounds[1], stop=params.ubounds[2], length=params.npatchs[1]+1) 
    vstops = range(params.vbounds[1], stop=params.vbounds[2], length=params.npatchs[2]+1) 
    AoT = Array{Union{ProductSpace, NTuple{3, Field}}}(undef, params.npatchs[1], params.npatchs[2])
    for index in CartesianIndices(AoT)
        AoT[index] = ProductSpace(ChebyshevGL{U, params.npoints[1], T}(ustops[index.I[1]], ustops[index.I[1]+1]),
                                   ChebyshevGL{V, params.npoints[2], T}(vstops[index.I[2]], vstops[index.I[2]+1]))
    end
    return AoT
end

function distribute(params::Parameters, excise::Function, 
                    computeUboundary::Function, computeVboundary::Function)::Array{Union{ProductSpace, NTuple{3, Field}}}
    AoT = setup(params)
    for index in CartesianIndices(AoT)
        if excise(AoT[index]) == true
            println(" Excising patch with bounds: ", range(AoT[index]))
        else
            println("Computing patch with bounds: ", range(fetch(AoT[index])))
            uboundary = index.I[1] == 1 ?  computeUboundary(AoT[index]) : extractUboundary(AoT[index - CartesianIndex((1,0))], :outgoing)
            vboundary = index.I[2] == 1 ?  computeVboundary(AoT[index]) : extractVboundary(AoT[index - CartesianIndex((0,1))], :outgoing)
            AoT[index] = compute(AoT[index], uboundary, vboundary)
        end
    end
    return AoT
end


