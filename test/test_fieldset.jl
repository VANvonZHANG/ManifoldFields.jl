using DimensionalData
using ManifoldFields
using ManifoldMeshes
using Test

@testset "location derivation" begin
    @test ManifoldFields._loc_from_layerdims((Dim{:node}(1:5),)) === NodeLoc
    @test ManifoldFields._loc_from_layerdims((Dim{:time}(1:3), Dim{:cell}(1:4))) === CellLoc
    @test ManifoldFields._loc_from_layerdims((Dim{:edge}(1:6),)) === EdgeLoc
    @test_throws ArgumentError ManifoldFields._loc_from_layerdims((Dim{:time}(1:3),))
    @test_throws ArgumentError ManifoldFields._loc_from_layerdims(
        (Dim{:node}(1:5), Dim{:cell}(1:4)))
end
