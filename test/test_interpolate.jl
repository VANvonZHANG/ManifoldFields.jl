using ManifoldFields
using ManifoldMeshes
using StaticArrays
using Test

@testset "NodeLoc interpolate" begin
    g = small_grid(; nlat=4, nlon=8)
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
    @test result ≈ [
        sum(node_time_values(g)[nodes[i], t] * weights[i] for i in eachindex(nodes))
        for t in 1:3
    ]

    ftn = DiscreteField(NodeLoc, g, time_node_values(g), time_node_dims(g))
    result2 = interpolate(ftn, lat, lon)
    @test result2 isa Vector
    @test length(result2) == 3
    @test result2 ≈ [
        sum(time_node_values(g)[t, nodes[i]] * weights[i] for i in eachindex(nodes))
        for t in 1:3
    ]

    @test_throws ArgumentError interpolate(f, 95.0, 0.0)
    @test_throws MethodError interpolate(
        DiscreteField(CellLoc, g, cell_values(g), cell_dims(g)),
        lat,
        lon,
    )
end
