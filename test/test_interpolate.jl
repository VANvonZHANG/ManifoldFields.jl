using ManifoldFields
using ManifoldMeshes
using StaticArrays
using Test

@testset "NodeLoc interpolate" begin
    g = small_grid(; nlat = 4, nlon = 8)
    values = node_values(g)
    f = DiscreteField(NodeLoc, g, values, node_dims(g))

    cid = 6
    lat, lon = cell_centroid(g, cid) |> ManifoldMeshes._cartesian_to_latlon
    nodes, weights = interpolation_weights(g, cid, lat, lon)
    expected = sum(values[nodes[i]] * weights[i] for i in eachindex(nodes))
    @test interpolate(f, lat, lon) ≈ expected

    p = SVector{3}(cell_centroid(g, cid))
    @test interpolate(f, p) ≈ expected

    ft = DiscreteField(NodeLoc, g, node_time_values(g), node_time_dims(g))
    result = interpolate(ft, lat, lon)
    @test result isa Vector
    @test length(result) == 3
    @test result ≈
          [sum(node_time_values(g)[nodes[i], t] * weights[i] for i in eachindex(nodes))
           for t in 1:3]

    ftn = DiscreteField(NodeLoc, g, time_node_values(g), time_node_dims(g))
    result2 = interpolate(ftn, lat, lon)
    @test result2 isa Vector
    @test length(result2) == 3
    @test result2 ≈
          [sum(time_node_values(g)[t, nodes[i]] * weights[i] for i in eachindex(nodes))
           for t in 1:3]

    @test_throws ArgumentError interpolate(f, 95.0, 0.0)
    @test_throws MethodError interpolate(
        DiscreteField(CellLoc, g, cell_values(g), cell_dims(g)),
        lat,
        lon
    )
end

@testset "NodeLoc batch interpolate" begin
    g = small_grid(; nlat = 4, nlon = 8)
    values = node_values(g)
    f = DiscreteField(NodeLoc, g, values, node_dims(g))

    lats = [-45.0, 0.0, 30.0, 60.0]
    lons = [10.0, 90.0, 200.0, 350.0]
    result = interpolate(f, lats, lons)

    @test result isa DimensionalData.AbstractDimArray
    @test size(result) == (4,)
    @test DimensionalData.name(result) == DimensionalData.name(f)
    @test result ≈ [interpolate(f, lats[i], lons[i]) for i in 1:4]

    # linear in latitude is reproduced by bilinear weight — independent check
    lin = DiscreteField(NodeLoc, g,
        [ManifoldMeshes._cartesian_to_latlon(node_coordinates(g, n))[1]
         for n in 1:num_nodes(g)],
        node_dims(g))
    @test interpolate(lin, lats, lons) ≈ lats rtol = 1e-12

    # trailing dimension is preserved, point axis first
    ft = DiscreteField(NodeLoc, g, node_time_values(g), node_time_dims(g))
    rt = interpolate(ft, lats, lons)
    @test size(rt) == (4, 3)
    for i in 1:4
        @test collect(selectdim(DimensionalData.data(rt), 1, i)) ≈
              interpolate(ft, lats[i], lons[i])
    end

    # time-first storage is handled as well
    ftn = DiscreteField(NodeLoc, g, time_node_values(g), time_node_dims(g))
    rtn = interpolate(ftn, lats, lons)
    @test size(rtn) == (4, 3)
    for i in 1:4
        @test collect(selectdim(DimensionalData.data(rtn), 1, i)) ≈
              interpolate(ftn, lats[i], lons[i])
    end

    # a mesh other than the field's own can be sampled pointwise
    other = small_grid(; nlat = 3, nlon = 6)
    lats2, lons2 = ManifoldMeshes._cartesian_to_latlon(cell_centroid(other, 3))
    @test interpolate(f, [lats2], [lons2])[1] ≈ interpolate(f, lats2, lons2)

    # empty query is empty, not an error
    @test length(interpolate(f, Float64[], Float64[])) == 0

    @test_throws DimensionMismatch interpolate(f, [0.0, 1.0], [0.0])
    @test_throws ArgumentError interpolate(f, [95.0], [0.0])
end
