#--------------------------------------------------------------------
# DNSS.jl
# Soham 03-2022
# Distribute the computation over multiple patches 
#--------------------------------------------------------------------

"""
    Set up an abstract computational product spaces given the grid parameters
    Input: Grid parameters 
    Output: Array of computational elements with their individual product spaces.
"""

function setup(params::Parameters{T}) where {T}
    ustops = range(params.ubounds[1], stop=params.ubounds[2], length=params.npatchs[1]+1) 
    vstops = range(params.vbounds[1], stop=params.vbounds[2], length=params.npatchs[2]+1) 
    AoT = Array{Union{ProductSpace, NTuple{params.nfields, Field}}}(undef, params.npatchs[1], params.npatchs[2])
    for index in CartesianIndices(AoT)
        AoT[index] = ProductSpace(ChebyshevGL{U, params.npoints[1], T}(ustops[index.I[1]], ustops[index.I[1]+1]),
                                  ChebyshevGL{V, params.npoints[2], T}(vstops[index.I[2]], vstops[index.I[2]+1]))
    end
    return AoT
end

"""
    Main code to compute the solution over multiple patches
    Input: grid parameters, excision function, uboundary, vboundary
    Output: Array of ProductSpace or NTuple of solutions
"""
function distribute(params::Parameters{T}, computePatch::Function, 
        computeUboundary::Function, computeVboundary::Function, debug::Int=0) where {T}
    AoT = setup(params)
    for index in CartesianIndices(AoT)
        if (debug >= 2) 
            println("Computing patch ", index.I)
            ((umin, vmin), (umax, vmax)) = range(AoT[index])
            @printf("  umin = %.2f umax = %.2f\n", umin, umax)
            @printf("  vmin = %.2f vmax = %.2f\n", vmin, vmax)
        end
        uboundary = index.I[1] == 1 ?  computeUboundary(AoT[index]) : extractUboundary(AoT[index - CartesianIndex((1,0))], :outgoing)
        vboundary = index.I[2] == 1 ?  computeVboundary(AoT[index]) : extractVboundary(AoT[index - CartesianIndex((0,1))], :outgoing)
        AoT[index] = computePatch(AoT[index], uboundary, vboundary)
    end
    return AoT
end

function setup_(params::Parameters{T})::Array{Union{ProductSpace, Future, NTuple{3, Field}}} where {T}
    ustops = range(params.ubounds[1], stop=params.ubounds[2], length=params.npatchs[1]+1) 
    vstops = range(params.vbounds[1], stop=params.vbounds[2], length=params.npatchs[2]+1) 
    AoT = Array{Union{ProductSpace, Future, NTuple{3, Field}}}(undef, params.npatchs[1], params.npatchs[2])
    for index in CartesianIndices(AoT)
        AoT[index] = @spawnat :any ProductSpace(ChebyshevGL{U, params.npoints[1], T}(ustops[index.I[1]], ustops[index.I[1]+1]),
                                         ChebyshevGL{V, params.npoints[2], T}(vstops[index.I[2]], vstops[index.I[2]+1]))
    end
    return AoT
end

function distribute_(params::Parameters{T}, computePatch::Function, 
                     computeUboundary::Function, computeVboundary::Function, debug::Int=0) where {T}
    function computeUboundary_(PS::Future)::Future where {S1, S2}
        @spawnat :any computeUboundary(fetch(PS))
    end
    
    function computeVboundary_(PS::Future)::Future where {S1, S2}
        @spawnat :any computeVboundary(fetch(PS))
    end
    
    function extractUboundary_(U::Future, boundarytype::Symbol)::Future
        @spawnat :any extractUboundary(fetch(U), boundarytype)
    end
    
    function extractVboundary_(U::Future, boundarytype::Symbol)::Future
        @spawnat :any extractVboundary(fetch(U), boundarytype)
    end
    
    function computePatch_(PS::Future, Uboundary::Future, Vboundary::Future)::Future where {S1, S2}
        @spawnat :any computePatch(fetch(PS), fetch(Uboundary), fetch(Vboundary))
    end
    
    AoT = setup_(params)
    for index in CartesianIndices(AoT)
        if (debug >= 2) 
            println("Computing patch ", index.I)
            println("  (umin, vmin) = ", range(fetch(AoT[index]))[1])
            println("  (umax, vmax) = ", range(fetch(AoT[index]))[2])
        end
        uboundary = index.I[1] == 1 ?  computeUboundary_(AoT[index]) : extractUboundary_(AoT[index - CartesianIndex((1,0))], :outgoing)
        vboundary = index.I[2] == 1 ?  computeVboundary_(AoT[index]) : extractVboundary_(AoT[index - CartesianIndex((0,1))], :outgoing)
        AoT[index] = computePatch_(AoT[index], uboundary, vboundary)
    end
    return fetch.(AoT)
end


