using Aqua
using ManifoldFields
using Test

@testset "Code quality (Aqua.jl)" begin
    Aqua.test_all(ManifoldFields)
end

@testset "ManifoldFields.jl" begin
end
