#-------------------------------------------------------------------- # DNSS.jl 
# Soham 03-2022
# Dump the data using NPZ files to plot using Python
#--------------------------------------------------------------------

using NPZ 

function save(fname::String, w::Field{S}) where {S} 
    u   = Field(w.space, u->u)
    wlm = basistransform(w) 
    npzwrite("output/data/$fname", Dict("u" => u.value, "w" => w.value, "wlm" => wlm.value))
end

function save(fname::String, w::Field{ProductSpace{S1, S2}}) where {S1, S2} 
    u   = Field(w.space.S1, u->u)
    v   = Field(w.space.S2, v->v)
    wlm = basistransform(w) 
    npzwrite("$fname", Dict("u" => u.value, "v" => v.value, "w" => w.value, "wlm" => wlm.value))
end

function save(fname::String, AoF::T) where {T<:Array{Field, 2}}
    for index in CartesianIndices(AoF)
        if eltype(AoF[index].value) <: Number
            (a, b) = index.I 
            save(join(["$fname", "$(10 + a)", "$(10 + b)"], "_"), AoF[index])
        else
            @show index
            "We found a element that does not contain a number?"
        end
    end
end

function save(fname::String, n::Vector, e::Vector)
    npzwrite("$fname", Dict("n" => n, "e" => e))
end
