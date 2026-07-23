using DimensionalData
using ManifoldFields
using ManifoldMeshes
using Test

@testset "DiscreteField construction and accessors" begin
    g = small_grid()

    nf = DiscreteField(NodeLoc, g, node_values(g), node_dims(g); name = :temperature,
        metadata = Dict("units" => "K"))
    ef = DiscreteField(EdgeLoc, g, edge_values(g), edge_dims(g); name = :flux)
    cf = DiscreteField(CellLoc, g, cell_values(g), cell_dims(g); name = :area)

    @test mesh(nf) === g
    @test data(nf) == node_values(g)
    @test location(nf) === NodeLoc
    @test DimensionalData.dims(nf) == node_dims(g)
    @test DimensionalData.name(nf) == :temperature
    @test DimensionalData.metadata(nf)["units"] == "K"
    @test parent(nf) == node_values(g)

    @test location(ef) === EdgeLoc
    @test location(cf) === CellLoc

    @test DiscreteField(NodeLoc, g, node_time_values(g), node_time_dims(g)) isa
          DiscreteField{NodeLoc}
    @test DiscreteField(NodeLoc, g, time_node_values(g), time_node_dims(g)) isa
          DiscreteField{NodeLoc}

    @test_throws ArgumentError DiscreteField(NodeLoc, g, node_values(g))
    @test_throws ArgumentError DiscreteField(NodeLoc, g, node_values(g),
        (Dim{:time}(1:num_nodes(g)),))
    @test_throws ArgumentError DiscreteField(NodeLoc, g, node_values(g),
        (Dim{:node}(1:num_nodes(g)),
            Dim{:node}(1:num_nodes(g))))
    @test_throws DimensionMismatch DiscreteField(NodeLoc, g, node_values(g)[1:(end - 1)],
        (Dim{:node}(1:(num_nodes(g) - 1)),))
    @test_throws DimensionMismatch DiscreteField(NodeLoc, g, node_time_values(g),
        (Dim{:node}(1:num_nodes(g)),
            Dim{:time}(1:2)))

    rebuilt = DimensionalData.rebuild(nf; data = 2 .* node_values(g), name = :doubled)
    @test rebuilt isa DiscreteField{NodeLoc}
    @test mesh(rebuilt) === g
    @test data(rebuilt) == 2 .* node_values(g)
    @test DimensionalData.dims(rebuilt) == node_dims(g)
    @test DimensionalData.name(rebuilt) == :doubled

    sim = similar(nf)
    @test sim isa DiscreteField{NodeLoc}
    @test mesh(sim) === g
    @test DimensionalData.dims(sim) == node_dims(g)
end
