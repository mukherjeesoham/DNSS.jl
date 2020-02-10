using DNSS, Test

libraries = ["Minkowski"]
libraries = ["Collapse"]
libraries = ["FourierGEP"]
libraries = ["BasisTransformation"]

for file in libraries
    @info "Testing $file"
    include("$(file)Test.jl")
end


