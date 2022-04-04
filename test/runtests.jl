using DNSS, Test

libraries = ["MinkowskiSphericalSymmetry"]

for file in libraries
    @info "Testing $file"
    include("$(file)Test.jl")
end


