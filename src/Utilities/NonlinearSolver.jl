#--------------------------------------------------------------------
# DNSS.jl
# Soham 09-2019
# Functions and function templates to interface with the
# nonlinear solver from NLsolve.jl
#--------------------------------------------------------------------

using NLsolve
export compute

function reshapeFromTuple(U::NTuple{3, Field})
    return vcat(reshape(U[1]), reshape(U[2]), reshape(U[3]))
end

function reshapeToTuple(space::S, x::Array{T,1})::NTuple{3, Field}  where {S, T}
    U = reshape(x, :, 3)
    return (reshape(space, U[:, 1]), reshape(space, U[:, 2]), reshape(space, U[:, 3]))
end

function compute(PS::ProductSpace{S1, S2}, ubnd::NTuple{3, Field{S2}}, 
                 vbnd::NTuple{3, Field{S1}})::NTuple{3, Field{ProductSpace{S1, S2}}} where {S1, S2}

    function f!(fvec::Array{T,1}, x::Array{T,1}) where {T}
        residual = mix!.(residual(reshapeToTuple(PS, x)), 
                         incomingboundary(PS), 
                         combineUVboundary(ubnd, vbnd))
        fvec[:] = reshapeFromTuple(residual)
    end

    us = reshapeToTuple(PS, nlsolve(f!, initialguess(PS); 
                                    method=:trust_region, 
                                    autodiff=:forward,
                                    show_trace=debug, 
                                    ftol=1e-10, 
                                    iterations=100).zero)

    println("\t L2(C1, C2) = ",  L2.(constraints(us...))) : 0
    return us
end

