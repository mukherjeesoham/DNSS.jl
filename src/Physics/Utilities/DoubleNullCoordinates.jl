#--------------------------------------------------------------------
# DNSS.jl
# Soham 06-2018
# Coordinate transformatins for Schwarzschild background in 
# double-null coordinates. See 'A Relativist's Toolkit, Eric Poisson'
#--------------------------------------------------------------------

using Roots
export find_r_of_UV

function find_t_of_UV(U::T, V::T, M::Number)::T where {T<:Number}
    @assert V > 0   # ensure you're in region I or II
    @assert U*V < T(1) # ensure you don't hit the singularity
    if U*V == 0     # r = 2M 
        t = V       # enforce uniqueness 
    elseif U > 0    # r < 2M
        t = -2M*log(U/V)
    elseif U < 0    # r > 2M
        t = -2M*log(-U/V)
    else
        error("Domain error")
    end
    return t
end

function find_r_of_UV(U::T, V::T, M::Number)::T where {T<:Number}
    @assert V > T(0)       # ensure you're in region I or II
    @assert U*V < T(1)     # ensure you don't hit the singularity
    if U*V == T(0)         # r = 2M 
        r = 2M
    else                # r < 2M or r > 2M 
        f(r) = (r/2M - 1)*exp(r/2M) + U*V 
        try
            # r    = find_zero(f, 2M)
            # NOTE: Introduce bracketing
            if U*V > T(0)
                r = find_zero(f, (eps(T), 2*M))
            else
                # WARNING: Would fail if r > 1000*M
                r = find_zero(f, (2M, 1000*M))
            end
        catch
            @show U, V
            exit()
        end
    end
    @assert r > 0
    # @assert r > 2M    # For testing.
    return r
end

function find_U_of_tr(t::T, r::T, M::Number)::T where {T<:Number}
    @assert r > 0
    if r == 2M
        U = 0
    else
        rstar = r + 2M*log(abs((r/2M)-1))
        u = t - rstar
        r > 2M ?  U = -exp(-u/4M) : U = +exp(-u/4M)
    end
    return U 
end

function find_V_of_tr(t::T, r::T, M::Number)::T where {T<:Number}
    @assert r > 0
    if r == 2M
        V = t           # make the inverse mapping unique. 
    else
        rstar = r + 2M*log(abs((r/2M)-1))
        v = t + rstar
        V = exp(v/4M)
    end
    return V
end
