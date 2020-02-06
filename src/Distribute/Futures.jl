#--------------------------------------------------------------------
# DNSS.jl
# Soham 09-2019
# Functions to distribute over multiple patches
#--------------------------------------------------------------------

import Base.Threads.@spawnat :any
export Fdistribute

function Fsetup(params::Parameters{T})::Array{Union{ProductSpace, Task, NTuple{3, Field}}} where {T}
    ustops = range(params.ubounds[1], stop=params.ubounds[2], length=params.npatchs[1]+1) 
    vstops = range(params.vbounds[1], stop=params.vbounds[2], length=params.npatchs[2]+1) 
    AoT = Array{Union{ProductSpace, Task, NTuple{3, Field}}}(undef, params.npatchs[1], params.npatchs[2])
    for index in CartesianIndices(AoT)
        AoT[index] = @spawnat :any ProductSpace(ChebyshevGL{U, params.npoints[1], T}(ustops[index.I[1]], ustops[index.I[1]+1]),
                                          ChebyshevGL{V, params.npoints[2], T}(vstops[index.I[2]], vstops[index.I[2]+1]))
    end
    return AoT
end

function Fdistribute(params::Parameters, excise::Function, 
                    computeUboundary::Function, computeVboundary::Function)::Array{Union{ProductSpace, Task, NTuple{3, Field}}}
    AoT = Fsetup(params)
    for index in CartesianIndices(AoT)
        if excise(fetch(AoT[index])) == true
            println(" Excising patch with bounds: ", range(fetch(AoT[index])))
        else
            println("Computing patch with bounds: ", range(fetch(AoT[index])))
            uboundary = index.I[1] == 1 ?  FcomputeUboundary(AoT[index]) : FextractUboundary(AoT[index - CartesianIndex((1,0))], :outgoing)
            vboundary = index.I[2] == 1 ?  FcomputeVboundary(AoT[index]) : FextractVboundary(AoT[index - CartesianIndex((0,1))], :outgoing)
            AoT[index] = Fcompute(AoT[index], uboundary, vboundary)
        end
    end
    return fetch.(AoT)
end

function FcomputeUboundary(PS::Task)::Task where {S1, S2}
    @spawnat :any computeUboundary(fetch(PS))
end

function FcomputeVboundary(PS::Task)::Task where {S1, S2}
    @spawnat :any computeVboundary(fetch(PS))
end

function FextractUboundary(U::Task, boundarytype::Symbol)::Task
    @spawnat :any extractUboundary(fetch(U), boundarytype)
end

function FextractVboundary(U::Task, boundarytype::Symbol)::Task
    @spawnat :any extractVboundary(fetch(U), boundarytype)
end

function Fcompute(PS::Task, Uboundary::Task, Vboundary::Task)::Task where {S1, S2}
    @spawnat :any compute(fetch(PS), fetch(Uboundary), fetch(Vboundary))
end

