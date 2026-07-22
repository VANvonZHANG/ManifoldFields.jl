using DimensionalData
using ManifoldFields
using ManifoldMeshes
using Test

include("helpers.jl")

@testset "DimensionalData compatibility spike" begin
    g = small_grid()
    f = DiscreteField(NodeLoc, g, node_time_values(g), node_time_dims(g);
                      name=:temperature, metadata=Dict("units" => "K"))

    @test dims(f) == node_time_dims(g)
    @test size(f) == (num_nodes(g), 3)
    @test f[1, 1] == node_time_values(g)[1, 1]
    @test DimensionalData.data(f) == node_time_values(g)
    @test DimensionalData.name(f) == :temperature
    @test DimensionalData.metadata(f)["units"] == "K"
    @test sum(f) == sum(node_time_values(g))
    @test extrema(f) == extrema(node_time_values(g))
end
