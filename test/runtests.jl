using Aqua
using ManifoldFields
using Test

@testset "Code quality (Aqua.jl)" begin
    Aqua.test_all(ManifoldFields)
end

@testset "ManifoldFields.jl" begin
    include("helpers.jl")
    include("test_locations.jl")
    include("test_field_construction.jl")
    include("test_dimensionaldata_spike.jl")
    include("test_broadcast.jl")
    include("test_interpolate.jl")
    include("test_ugrid.jl")
end
