#--------------------------------------------------------------------
# DNSS.jl
# Soham 08-2019
# Test scalar field on Minkowski spacetime
#--------------------------------------------------------------------

using NLsolve, PyPlot

# Set up initial conditions on u0 and v0
function ubnd(PS::ProductSpace{S1, S2})::NTuple{1, Field{S2}} where {S1, S2}
    return (Field(PS.S2, v->sinpi(v) * exp(-v^2)), )
end

function vbnd(PS::ProductSpace{S1, S2})::NTuple{1, Field{S1}} where {S1, S2}
    return (Field(PS.S1, u->sinpi(u) * exp(-u^2)), )
end

# Set up computation on each patch
function compute(PS::ProductSpace{S1, S2}, ubnd::NTuple{1,Field{S2}}, vbnd::NTuple{1,Field{S1}})::NTuple{1, Field{ProductSpace{S1, S2}}} where {S1, S2}

    # Define all the operators we use local to the patch. 
    B = incomingboundary(PS)
    I = identity(PS)
    DU, DV = derivative(PS)
    ϕbnd = combineUVboundary(ubnd, vbnd, :incoming)

    # Compute the residual
    # NOTE: This is slow but works. Sets the residual in the interior 
    # as well as the boundary.
    function residual(ϕ::Field{S})::NTuple{1, Field{S}} where {S}
        F = (I - B) * (DU*(DV*ϕ)) + B * (ϕ - first(ϕbnd))
        return (F,)
    end


    # Specify an initial guess
    function initialguess(PS::ProductSpace{S1, S2})::NTuple{1, Field{ProductSpace{S1, S2}}} where {S1, S2}
        return (Field(PS, (u,v)->1 + rand()), ) 
    end

    # Residual to pass to the non-linear solver in NLsolve
    function f!(fvec::Array{T,1}, x::Array{T,1}) where {T}
        fvec[:]  = reshapeFromTuple(residual(reshapeToTuple(PS, 1, x)...))
    end

    # Compute the initial guess
    U0 = reshapeFromTuple(initialguess(PS))

    # Call the non-linear solver
    U = nlsolve(f!, U0; method=:trust_region, autodiff=:forward, show_trace=true, ftol=1e-9, iterations=100)
    return reshapeToTuple(PS, 1, U.zero)
end

# Setup simulation grid parameters 
npoints = (12, 12)
npatches = (8, 8) 
urange = (-1.0, 1.0)
vrange = (-1.0, 1.0) 
nfields  = 1
Grid  = Parameters(npoints, npatches, urange, vrange, nfields)

# Compute the scalar field over the grid
ϕ_ = distribute(Grid, compute, ubnd, vbnd)
ϕ  = extract(ϕ_, 1)

# Plot the scalar field
contourf(ϕ, 20)
savefig("output/phi_minkowski.pdf")
close()
