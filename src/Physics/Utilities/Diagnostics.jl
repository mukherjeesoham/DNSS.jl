#--------------------------------------------------------------------
# DNSS.jl
# Soham 04-2019
# Compute diagnostics for spacetime
# TODO: Add Ricci scalar and Kretschmann scalar
#--------------------------------------------------------------------

export outgoingexpansion

function outgoingexpansion(a::Field{S}, η::Field{S})::Field{S} where {S}
    DU, DV = derivative(a.space)
    return 2*(DV*η)/(abs(a^2)*η)
end
