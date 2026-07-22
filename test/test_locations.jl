using DimensionalData
using ManifoldMeshes
using ManifoldFields
using Test

@testset "location trait methods" begin
    g = LatLonGrid(
        lat_edges = collect(range(-90.0, 90.0; length = 3)),
        lon_edges = collect(range(0.0, 360.0; length = 5))
    )

    @test location_dimname(NodeLoc) == Dim{:node}
    @test location_dimname(EdgeLoc) == Dim{:edge}
    @test location_dimname(CellLoc) == Dim{:cell}

    @test ManifoldFields.expected_location_length(NodeLoc, g) == num_nodes(g)
    @test ManifoldFields.expected_location_length(EdgeLoc, g) == num_edges(g)
    @test ManifoldFields.expected_location_length(CellLoc, g) == num_cells(g)

    @test ManifoldFields.ugrid_location(NodeLoc) == "node"
    @test ManifoldFields.ugrid_location(EdgeLoc) == "edge"
    @test ManifoldFields.ugrid_location(CellLoc) == "face"

    @test ManifoldFields.ugrid_dimname(NodeLoc) == "n_node"
    @test ManifoldFields.ugrid_dimname(EdgeLoc) == "n_edge"
    @test ManifoldFields.ugrid_dimname(CellLoc) == "n_face"
end
