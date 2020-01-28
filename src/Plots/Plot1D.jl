#--------------------------------------------------------------------
# DNSS.jl 
# Soham 01-2019
# Plotting routines for 1D
#--------------------------------------------------------------------

using PyPlot
export plot

function PyPlot. plot(u::Field{S}; plotstyle="-o", label="") where {S}
    x = Field(u.space, x->x)
    plot(x.value, u.value, plotstyle, label=label) 
end

