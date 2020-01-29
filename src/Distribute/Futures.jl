#--------------------------------------------------------------------
# DNSS.jl
# Soham 09-2019
# Functions to distribute over multiple patches
#--------------------------------------------------------------------

import Base.Threads.@spawn
export Fdistribute

function Fsetup(params::Parameters{T})::Array{Union{ProductSpace, Task, NTuple{3, Field}}} where {T}
    ustops = range(params.ubounds[1], stop=params.ubounds[2], length=params.npatchs[1]+1) 
    vstops = range(params.vbounds[1], stop=params.vbounds[2], length=params.npatchs[2]+1) 
    tree   = Array{Union{ProductSpace, Task, NTuple{3, Field}}}(undef, params.npatchs[1], params.npatchs[2])
    for index in CartesianIndices(tree)
        tree[index] = @spawn ProductSpace(ChebyshevGL{U, params.npoints[1], T}(ustops[index.I[1]], ustops[index.I[1]+1]),
                                          ChebyshevGL{V, params.npoints[2], T}(vstops[index.I[2]], vstops[index.I[2]+1]))
    end
    return tree
end

function Fdistribute(params::Parameters, excise::Function, 
                    computeUboundary::Function, computeVboundary::Function)::Array{Union{ProductSpace, Task, NTuple{3, Field}}}
    tree = Fsetup(params)
    for index in CartesianIndices(tree)
        if excise(fetch(tree[index])) == true
            println(" Excising patch with bounds: ", range(fetch(tree[index])))
        else
            println("Computing patch with bounds: ", range(fetch(tree[index])))
            uboundary = index.I[1] == 1 ?  FcomputeUboundary(tree[index]) : FextractUboundary(tree[index - CartesianIndex((1,0))], :outgoing)
            vboundary = index.I[2] == 1 ?  FcomputeVboundary(tree[index]) : FextractVboundary(tree[index - CartesianIndex((0,1))], :outgoing)
            tree[index] = Fcompute(tree[index], uboundary, vboundary)
        end
    end
    return fetch.(tree)
end

function FcomputeUboundary(PS::Task)::Task where {S1, S2}
    @spawn computeUboundary(fetch(PS))
end

function FcomputeVboundary(PS::Task)::Task where {S1, S2}
    @spawn computeVboundary(fetch(PS))
end

function FextractUboundary(U::Task, boundarytype::Symbol)::Task
    @spawn extractUboundary(fetch(U), boundarytype)
end

function FextractVboundary(U::Task, boundarytype::Symbol)::Task
    @spawn extractVboundary(fetch(U), boundarytype)
end

function Fcompute(PS::Task, Uboundary::Task, Vboundary::Task)::Task where {S1, S2}
    @spawn compute(fetch(PS), fetch(Uboundary), fetch(Vboundary))
end

