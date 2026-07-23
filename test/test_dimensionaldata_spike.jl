using DimensionalData
using ManifoldFields
using ManifoldMeshes
using Test

@testset "DimensionalData compatibility spike" begin
    g = small_grid()
    f = DiscreteField(NodeLoc, g, node_time_values(g), node_time_dims(g);
        name = :temperature, metadata = Dict("units" => "K"))

    @test dims(f) == node_time_dims(g)
    @test size(f) == (num_nodes(g), 3)
    @test f[1, 1] == node_time_values(g)[1, 1]
    @test DimensionalData.data(f) == node_time_values(g)
    @test DimensionalData.name(f) == :temperature
    @test DimensionalData.metadata(f)["units"] == "K"
    @test sum(f) == sum(node_time_values(g))
    @test extrema(f) == extrema(node_time_values(g))
end

@testset "DiscreteField DimensionalData selector indexing" begin
    g = small_grid()

    f = DiscreteField(NodeLoc, g, node_time_values(g), node_time_dims(g);
        name = :temperature)
    time_slice = f[time = At(2)]
    @test time_slice isa DiscreteField{NodeLoc}
    @test mesh(time_slice) === g
    @test data(time_slice) == node_time_values(g)[:, 2]
    @test dims(time_slice) == node_dims(g)

    dim_selector_slice = f[Dim{:time}(At(2))]
    @test dim_selector_slice isa DiscreteField{NodeLoc}
    @test mesh(dim_selector_slice) === g
    @test data(dim_selector_slice) == node_time_values(g)[:, 2]
    @test dims(dim_selector_slice) == node_dims(g)

    ftn = DiscreteField(NodeLoc, g, time_node_values(g), time_node_dims(g);
        name = :temperature)
    time_first_slice = ftn[time = At(2)]
    @test time_first_slice isa DiscreteField{NodeLoc}
    @test mesh(time_first_slice) === g
    @test data(time_first_slice) == time_node_values(g)[2, :]
    @test dims(time_first_slice) == node_dims(g)
end
