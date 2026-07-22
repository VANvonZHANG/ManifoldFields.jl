using ManifoldFields
using ManifoldMeshes
using Test

@testset "UGRID topology writer" begin
    g = small_grid()
    ds = to_ugrid(g)

    @test ds.attributes["Conventions"] == "CF-1.11 UGRID-1.0"

    @test haskey(ds.variables, "Mesh2")
    meshvar = ds.variables["Mesh2"]
    @test meshvar.attrs["cf_role"] == "mesh_topology"
    @test meshvar.attrs["topology_dimension"] == 2
    @test meshvar.attrs["node_coordinates"] == "Mesh2_node_lon Mesh2_node_lat"
    @test meshvar.attrs["face_node_connectivity"] == "Mesh2_face_nodes"
    @test meshvar.attrs["face_dimension"] == "n_face"

    node_lon = ds.variables["Mesh2_node_lon"]
    @test node_lon.dims == ("n_node",)
    @test size(node_lon.data) == (num_nodes(g),)
    @test node_lon.attrs["standard_name"] == "longitude"
    @test node_lon.attrs["units"] == "degrees_east"

    node_lat = ds.variables["Mesh2_node_lat"]
    @test node_lat.dims == ("n_node",)
    @test size(node_lat.data) == (num_nodes(g),)
    @test node_lat.attrs["standard_name"] == "latitude"
    @test node_lat.attrs["units"] == "degrees_north"

    face_nodes = ds.variables["Mesh2_face_nodes"]
    @test face_nodes.dims == ("n_face", "n_max_face_nodes")
    @test size(face_nodes.data) == (num_cells(g), 4)
    @test face_nodes.attrs["cf_role"] == "face_node_connectivity"
    @test face_nodes.attrs["start_index"] == 1
    for c in 1:num_cells(g)
        @test Tuple(face_nodes.data[c, :]) == cell_nodes(g, c)
    end
end
