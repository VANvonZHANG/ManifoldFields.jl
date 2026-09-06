using DimensionalData
using ManifoldFields
using ManifoldMeshes
using Statistics
using Test

@testset "location derivation" begin
    @test ManifoldFields._loc_from_layerdims((Dim{:node}(1:5),)) === NodeLoc
    @test ManifoldFields._loc_from_layerdims((Dim{:time}(1:3), Dim{:cell}(1:4))) === CellLoc
    @test ManifoldFields._loc_from_layerdims((Dim{:edge}(1:6),)) === EdgeLoc
    @test_throws ArgumentError ManifoldFields._loc_from_layerdims((Dim{:time}(1:3),))
    @test_throws ArgumentError ManifoldFields._loc_from_layerdims(
        (Dim{:node}(1:5), Dim{:cell}(1:4)))
end

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

@testset "FieldSet access" begin
    g = small_grid()
    fs = sample_fieldset(g)

    u = fs[:u]
    @test u isa DiscreteField{NodeLoc}
    @test mesh(u) === g
    @test data(u) == node_time_values(g)
    @test DimensionalData.dims(u) == node_time_dims(g)
    @test DimensionalData.name(u) == :u

    v = fs.v
    @test v isa DiscreteField{CellLoc}
    @test mesh(v) === g

    all = fields(fs)
    @test all isa NamedTuple{(:u, :v)}
    @test all[:v] isa DiscreteField{CellLoc}

    @test_throws KeyError fs[:missing]

    io = IOBuffer()
    show(io, fs)
    str = String(take!(io))
    @test occursin("FieldSet", str)
    @test occursin(":u", str)

    # linear indexing pins DimensionalData semantics (per-field scalars, not an error)
    @test fs[1] isa NamedTuple
    @test keys(fs[1]) == (:u, :v)
end

@testset "FieldSet dimensional slicing" begin
    g = small_grid()
    u = DiscreteField(NodeLoc, g, node_time_values(g), node_time_dims(g); name = :u)
    w = DiscreteField(NodeLoc, g, time_node_values(g), time_node_dims(g); name = :w)
    fs = FieldSet(g, :u => u, :w => w)   # both layers have :time

    sliced = fs[Dim{:time}(1:2)]
    @test sliced isa FieldSet
    @test mesh(sliced) === g
    @test size(DimensionalData.data(sliced)[:u]) == (num_nodes(g), 2)
    @test size(DimensionalData.data(sliced)[:w]) == (2, num_nodes(g))

    @test_throws ArgumentError fs[Dim{:node}(1:2)]
    @test_throws ArgumentError fs[Dim{:cell}(1:2)]   # not present anywhere

    # nowhere-present, non-location dim: warn and return an unchanged copy
    fs_typo = @test_warn "not found in any field" fs[Dim{:typo}(1:2)]
    @test fs_typo isa FieldSet
    @test fs_typo == fs
end

@testset "FieldSet rebuild overrides" begin
    g = small_grid()
    fs = sample_fieldset(g)

    r1 = DimensionalData.rebuild(fs)
    @test r1 isa FieldSet
    @test mesh(r1) === g
    @test data(r1[:v]) == data(fs[:v])

    r2 = DimensionalData.rebuild(fs, DimensionalData.data(fs))
    @test r2 isa FieldSet
    @test mesh(r2) === g

    das = fields(fs)
    r3 = DimensionalData.rebuild_from_arrays(fs, das)
    @test r3 isa FieldSet
    @test mesh(r3) === g
    @test field_names(r3) == (:u, :v)

    r4 = DimensionalData.rebuild_from_arrays(fs, Tuple(das))
    @test r4 isa FieldSet
    @test mesh(r4) === g
end

@testset "FieldSet broadcast" begin
    g = small_grid()
    fs = sample_fieldset(g)

    fs2 = 2 .* fs
    @test fs2 isa FieldSet
    @test mesh(fs2) === g
    @test data(fs2[:u]) == 2 .* node_time_values(g)
    @test data(fs2[:v]) == 2 .* cell_values(g)

    fs3 = sin.(fs)
    @test fs3 isa FieldSet
    @test data(fs3[:v]) == sin.(cell_values(g))

    fs4 = fs .+ fs
    @test fs4 isa FieldSet
    @test data(fs4[:u]) == 2 .* node_time_values(g)

    # mixed locations coexist through broadcast
    @test location(fs4[:u]) === NodeLoc
    @test location(fs4[:v]) === CellLoc

    # same-name-different-location pairing fails (u is NodeLoc here, NodeLoc ok;
    # construct a CellLoc field named :u to force the mismatch)
    u_cell = DiscreteField(CellLoc, g, cell_values(g), cell_dims(g); name = :u)
    v_cell = DiscreteField(CellLoc, g, cell_values(g), cell_dims(g); name = :v)
    fs_cell = FieldSet(g, :u => u_cell, :v => v_cell)
    @test_throws DimensionMismatch fs .+ fs_cell

    # name-set mismatch
    other = FieldSet(g, :u => DiscreteField(
        NodeLoc, g, node_time_values(g), node_time_dims(g); name = :u))
    @test_throws DimensionMismatch fs .+ other

    # mesh mismatch
    g2 = small_grid()
    fs_g2 = FieldSet(g2,
        :u => DiscreteField(
            NodeLoc, g2, node_time_values(g2), node_time_dims(g2); name = :u),
        :v => DiscreteField(CellLoc, g2, cell_values(g2), cell_dims(g2); name = :v))
    @test_throws DimensionMismatch fs .+ fs_g2

    # dim-order mismatch within same field name (consistent with DiscreteField)
    u_tn = DiscreteField(NodeLoc, g, time_node_values(g), time_node_dims(g); name = :u)
    fs_tn = FieldSet(g, :u => u_tn, :v => DiscreteField(
        CellLoc, g, cell_values(g), cell_dims(g); name = :v))
    @test_throws DimensionMismatch fs .+ fs_tn
end

@testset "FieldSet merge and equality" begin
    g = small_grid()
    fs = sample_fieldset(g)
    w = DiscreteField(EdgeLoc, g, edge_values(g), edge_dims(g); name = :w)

    fs3 = merge(fs, :w => w)
    @test field_names(fs3) == (:u, :v, :w)
    @test fs3[:w] isa DiscreteField{EdgeLoc}
    @test mesh(fs3) === g

    fs2 = merge(fs, (w = w,))
    @test field_names(fs2) == (:u, :v, :w)

    other = small_grid(; nlat = 3)
    w_wrong = DiscreteField(
        EdgeLoc, other, collect(Float64, 1:num_edges(other)), edge_dims(other); name = :w)
    @test_throws DimensionMismatch merge(fs, :w => w_wrong)
    @test_throws DimensionMismatch merge(fs, (w = w_wrong,))

    fs_copy = merge(fs)
    @test field_names(fs_copy) == (:u, :v)

    @test fs == sample_fieldset(g)
    @test fs != 2 .* fs
    @test fs != FieldSet(g, :x => DiscreteField(
        NodeLoc, g, node_time_values(g), node_time_dims(g); name = :x))
    g3 = small_grid(; nlat = 3)
    fs_g3 = FieldSet(g3,
        :u => DiscreteField(
            NodeLoc, g3, node_time_values(g3), node_time_dims(g3); name = :u),
        :v => DiscreteField(CellLoc, g3, cell_values(g3), cell_dims(g3); name = :v))
    @test fs != fs_g3

    # stack-stack merge validates mesh identity (same counts, different objects)
    g_other = small_grid()   # same nlat/nlon as g, distinct object
    fs_other = FieldSet(g_other,
        :u => DiscreteField(
            NodeLoc, g_other, node_time_values(g_other), node_time_dims(g_other); name = :u),
        :v => DiscreteField(
            CellLoc, g_other, cell_values(g_other), cell_dims(g_other); name = :v))
    @test_throws DimensionMismatch merge(fs, fs_other)

    # setindex! style update validates mesh identity
    w_other = DiscreteField(
        EdgeLoc, g_other, edge_values(g_other), edge_dims(g_other); name = :w)
    @test_throws DimensionMismatch Base.setindex(fs, w_other, :w)

    # valid setindex adds the field
    fs_set = Base.setindex(fs, w, :w)
    @test field_names(fs_set) == (:u, :v, :w)
    @test mesh(fs_set) === g
end

@testset "FieldSet reductions" begin
    g = small_grid()
    fs = sample_fieldset(g)   # :u has dims (node, time); :v has dims (cell,)

    s = sum(fs; dims = Dim{:time})
    @test s isa FieldSet
    @test mesh(s) === g
    @test data(s[:u]) == dropdims(
        sum(node_time_values(g); dims = 2); dims = 2)

    total = sum(fs)
    @test total isa NamedTuple
    @test total[:v] == sum(cell_values(g))

    @test_throws DimensionMismatch sum(fs; dims = Dim{:node})
end

@testset "FieldSet keepdims reductions" begin
    g = small_grid()
    fs = sample_fieldset(g)

    s = sum(fs; dims = Dim{:time})                       # default: drop
    @test DimensionalData.dims(s[:u]) == (Dim{:node}(1:num_nodes(g)),)
    @test data(s[:u]) == vec(sum(node_time_values(g); dims = 2))

    sk = sum(fs; dims = Dim{:time}, keepdims = true)     # keep length-1
    dsk = DimensionalData.dims(sk[:u])
    @test map(DimensionalData.name, dsk) == (:node, :time)   # time retained, not dropped
    @test length(dsk[2]) == 1   # stock DimensionalData reducelookup keeps a value, not 1:1
    @test size(data(sk[:u])) == (num_nodes(g), 1)

    fk = Statistics.mean(fs; dims = Dim{:time}, keepdims = true)
    @test size(data(fk[:u])) == (num_nodes(g), 1)
end

@testset "FieldSet cat mesh validation" begin
    g = small_grid()
    fs1 = sample_fieldset(g)
    fs2 = 2 .* sample_fieldset(g)

    c = Base.cat(fs1, fs2; dims = Dim{:time})
    @test c isa FieldSet
    @test mesh(c) === g
    @test size(data(c[:u])) == (num_nodes(g), 6)

    g_other = small_grid()   # same counts, different object
    fs3 = FieldSet(g_other,
        :u => DiscreteField(
            NodeLoc, g_other, node_time_values(g), node_time_dims(g); name = :u),
        :v => DiscreteField(CellLoc, g_other, cell_values(g), cell_dims(g); name = :v))
    @test_throws DimensionMismatch Base.cat(fs1, fs3; dims = Dim{:time})
end
