using DimensionalData
using ManifoldFields
using ManifoldMeshes
using Test

include("helpers.jl")

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

    s = sin.(f)
    @test s isa DiscreteField{NodeLoc}
    @test mesh(s) === g
    @test data(s) == sin.(node_values(g))

    t = 2 .* f
    @test t isa DiscreteField{NodeLoc}
    @test data(t) == 2 .* node_values(g)

    g2 = small_grid()
    same_shape_other_mesh = DiscreteField(NodeLoc, g2, node_values(g2), node_dims(g2))
    @test_throws DimensionMismatch f .+ same_shape_other_mesh

    cell_f = DiscreteField(CellLoc, g, cell_values(g), cell_dims(g))
    @test_throws DimensionMismatch f .+ cell_f

    rebound = withmesh(same_shape_other_mesh, mesh(f))
    @test data(f .+ rebound) == node_values(g) .+ node_values(g2)
end
