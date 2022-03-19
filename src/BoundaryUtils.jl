#--------------------------------------------------------------------
# DNSS.jl
# Soham 04-2019
# Boundary utilities
#--------------------------------------------------------------------

export extractUboundary, extractVboundary, combineUVboundary

function extractUboundary(u::Field{ProductSpace{S1, S2}}, boundarytype::Symbol)::Field{S2} where {S1, S2}
    @assert ndims(u.value) == 2
    return (boundarytype == :incoming ? Field(u.space.S2, u.value[end, :]) : Field(u.space.S2, u.value[1, :]))
end

function extractVboundary(u::Field{ProductSpace{S1, S2}}, boundarytype::Symbol)::Field{S1} where {S1, S2}
    @assert ndims(u.value) == 2
    return (boundarytype == :incoming ? Field(u.space.S1, u.value[:, end]) : Field(u.space.S1, (u.value[:, 1])))
end

function combineUVboundary(uboundary::Field{S2}, vboundary::Field{S1}, boundarytype::Symbol)::Field{ProductSpace{S1, S2}} where {S1, S2}
    PS = ProductSpace(vboundary.space, uboundary.space)
    uv  = Field(PS, (u,v)->0)    
    if boundarytype == :incoming
        uv.value[end, :] = uboundary.value
        uv.value[:, end] = vboundary.value
    else 
        uv.value[1, :]   = uboundary.value
        uv.value[:, 1]   = vboundary.value
    end
    return uv 
end

function extractUboundary(u::NTuple{N, Field{ProductSpace{S1, S2}}}, boundarytype::Symbol)::NTuple{N, Field{S2}} where {S1, S2, N}
    return extractUboundary.(u, boundarytype) 
end

function extractVboundary(u::NTuple{N, Field{ProductSpace{S1, S2}}}, boundarytype::Symbol)::NTuple{N, Field{S1}} where {S1, S2, N}
    return extractVboundary.(u, boundarytype) 
end

function combineUVboundary(ubnd::NTuple{N, Field{S2}},
                           vbnd::NTuple{N, Field{S1}}, boundarytype::Symbol)::NTuple{N, Field{ProductSpace{S1, S2}}} where {S1, S2, N}
    return combineUVboundary.(ubnd, vbnd, boundarytype) 
end

