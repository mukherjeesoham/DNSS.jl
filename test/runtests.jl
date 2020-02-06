using DNSS, Test

libraries = ["Plots"]

for file in libraries
    @info "Testing $file"
    include("$(file)Test.jl")
end


