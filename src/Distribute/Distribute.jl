#--------------------------------------------------------------------
# DNSS.jl
# Soham 09-2019
# Functions to distribute over multiple patches
#--------------------------------------------------------------------

export Parameters, distribute

struct Parameters{T}
    npoints::NTuple{2, Int}
    npatchs::NTuple{2, Int}
    ubounds::NTuple{2, T}
    vbounds::NTuple{2, T}
end

function setup(params::Parameters{T})::Array{Union{ProductSpace, NTuple{3, Field}}} where {T}
    ustops = range(params.ubounds[1], stop=params.ubounds[2], length=params.npatchs[1]+1) 
    vstops = range(params.vbounds[1], stop=params.vbounds[2], length=params.npatchs[2]+1) 
    tree   = Array{Union{ProductSpace, NTuple{3, Field}}}(undef, params.npatchs[1], params.npatchs[2])
    for index in CartesianIndices(tree)
        tree[index] = ProductSpace(ChebyshevGL{U, params.npoints[1], T}(ustops[index.I[1]], ustops[index.I[1]+1]),
                                   ChebyshevGL{V, params.npoints[2], T}(vstops[index.I[2]], vstops[index.I[2]+1]))
    end
    return tree
end

function distribute(params::Parameters, excise::Function, 
                    computeUboundary::Function, computeVboundary::Function)::Array{Union{ProductSpace, NTuple{3, Field}}}
    tree = setup(params)
    for index in CartesianIndices(tree)
        if excise(tree[index]) == true
            println(" Excising patch with bounds: ", range(tree[index]))
        else
            println("Computing patch with bounds: ", range(fetch(tree[index])))
            uboundary = index.I[1] == 1 ?  computeUboundary(tree[index]) : extractUboundary(tree[index - CartesianIndex((1,0))], :outgoing)
            vboundary = index.I[2] == 1 ?  computeVboundary(tree[index]) : extractVboundary(tree[index - CartesianIndex((0,1))], :outgoing)
            tree[index] = compute(tree[index], uboundary, vboundary)
        end
    end
    return tree
end


