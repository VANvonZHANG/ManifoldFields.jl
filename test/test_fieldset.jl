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

using ManifoldFields: FieldSet, field_names

@testset "FieldSet construction" begin
    g = small_grid()
    u = DiscreteField(NodeLoc, g, node_time_values(g), node_time_dims(g); name = :u)
    v = DiscreteField(CellLoc, g, cell_values(g), cell_dims(g); name = :v)

    fs = FieldSet(g, (u = u, v = v))
    @test fs isa FieldSet
    @test mesh(fs) === g
    @test field_names(fs) == (:u, :v)

    fs_pairs = FieldSet(g, :u => u, :v => v)
    @test field_names(fs_pairs) == (:u, :v)

    fs_dict = FieldSet(g, Dict(:u => u, :v => v))
    @test Set(field_names(fs_dict)) == Set([:u, :v])
    @test issorted(field_names(fs_dict))

    fs_conv = FieldSet(:u => u, :v => v)
    @test mesh(fs_conv) === g

    w = DiscreteField(EdgeLoc, g, edge_values(g), edge_dims(g); name = :w)
    fs3 = FieldSet(g, :u => u, :v => v, :w => w)
    @test field_names(fs3) == (:u, :v, :w)
end

@testset "FieldSet construction validation" begin
    g = small_grid()
    other = small_grid(; nlat = 3)
    u = DiscreteField(NodeLoc, g, node_time_values(g), node_time_dims(g); name = :u)

    @test_throws ArgumentError FieldSet(g, NamedTuple())
    @test_throws ArgumentError FieldSet()

    wrong_mesh = DiscreteField(
        NodeLoc, other, collect(Float64, 1:num_nodes(other)), node_dims(other); name = :u)
    @test_throws DimensionMismatch FieldSet(g, :u => wrong_mesh)

    # a layer whose dims contain no location dimension is rejected by the
    # inner-constructor check (exercised via the direct fieldwise constructor)
    @test_throws ArgumentError ManifoldFields.FieldSet(
        (plain = collect(Float64, 1:3),),
        (Dim{:time}(1:3),),
        (),
        (plain = (Dim{:time}(1:3),),),
        DimensionalData.NoMetadata(),
        (plain = DimensionalData.NoMetadata(),),
        g
    )

    # location dim length mismatch vs mesh (inner-constructor path)
    n = num_nodes(g)
    @test_throws DimensionMismatch ManifoldFields.FieldSet(
        (u = zeros(n - 1, 3),),
        (Dim{:node}(1:(n - 1)), Dim{:time}(1:3)),
        (),
        (u = (Dim{:node}(1:(n - 1)), Dim{:time}(1:3)),),
        DimensionalData.NoMetadata(),
        (u = DimensionalData.NoMetadata(),),
        g
    )
end
