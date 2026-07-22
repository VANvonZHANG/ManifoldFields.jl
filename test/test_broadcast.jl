using DimensionalData
using ManifoldFields
using ManifoldMeshes
using Test

@testset "DiscreteField broadcast arithmetic" begin
    g = small_grid()
    f = DiscreteField(NodeLoc, g, node_values(g), node_dims(g); name=:a,
                      metadata=Dict("units" => "K"))
    h = DiscreteField(NodeLoc, g, 2 .* node_values(g), node_dims(g); name=:b,
                      metadata=Dict("units" => "K"))

    r = f .+ h
    @test r isa DiscreteField{NodeLoc}
    @test mesh(r) === g
    @test location(r) === NodeLoc
    @test data(r) == 3 .* node_values(g)
    @test dims(r) == node_dims(g)
    @test DimensionalData.name(r) == :a
    @test DimensionalData.metadata(r)["units"] == "K"

    s = sin.(f)
    @test s isa DiscreteField{NodeLoc}
    @test mesh(s) === g
    @test data(s) == sin.(node_values(g))

    t = 2 .* f
    @test t isa DiscreteField{NodeLoc}
    @test data(t) == 2 .* node_values(g)
    @test DimensionalData.name(t) == :a

    arr_result = f .+ ones(num_nodes(g))
    @test arr_result isa DiscreteField{NodeLoc}
    @test data(arr_result) == node_values(g) .+ 1

    dim_array = DimArray(ones(num_nodes(g)), node_dims(g))
    dim_result = f .+ dim_array
    @test dim_result isa DiscreteField{NodeLoc}
    @test mesh(dim_result) === g
    @test data(dim_result) == node_values(g) .+ 1

    g2 = small_grid()
    same_shape_other_mesh = DiscreteField(NodeLoc, g2, node_values(g2), node_dims(g2))
    @test_throws DimensionMismatch f .+ same_shape_other_mesh
    @test_throws DimensionMismatch f .+ same_shape_other_mesh .+ dim_array

    cell_f = DiscreteField(CellLoc, g, cell_values(g), cell_dims(g))
    @test_throws DimensionMismatch f .+ cell_f

    rebound = withmesh(same_shape_other_mesh, mesh(f))
    @test data(f .+ rebound) == node_values(g) .+ node_values(g2)
end
