using DNSS, Test

libraries = ["Minkowski"]

for file in libraries
    @info "Testing $file"
    include("$(file)Test.jl")
end


