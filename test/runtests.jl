using DNSS, Test

libraries = ["Minkowski", "Collapse"]

for file in libraries
    @info "Testing $file"
    include("$(file)Test.jl")
end


