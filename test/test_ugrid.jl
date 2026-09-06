using DimensionalData
using LinearAlgebra
using ManifoldFields
using ManifoldMeshes
using StaticArrays
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

@testset "UGRID field writer" begin
    g = small_grid()

    nf = DiscreteField(NodeLoc, g, node_values(g), node_dims(g); name = :node_temp)
    node_ds = to_ugrid(nf)
    @test haskey(node_ds.variables, "node_temp")
    node_temp = node_ds.variables["node_temp"]
    @test node_temp.attrs["mesh"] == "Mesh2"
    @test node_temp.attrs["location"] == "node"
    @test node_temp.attrs["coordinates"] == "Mesh2_node_lon Mesh2_node_lat"
    @test node_temp.dims == ("n_node",)
    @test node_temp.data == data(nf)

    cf = DiscreteField(CellLoc, g, cell_values(g), cell_dims(g); name = :cell_area)
    cell_ds = to_ugrid(cf)
    @test haskey(cell_ds.variables, "cell_area")
    cell_area = cell_ds.variables["cell_area"]
    @test cell_area.attrs["mesh"] == "Mesh2"
    @test cell_area.attrs["location"] == "face"
    @test cell_area.attrs["coordinates"] == "Mesh2_face_lon Mesh2_face_lat"
    @test cell_area.dims == ("n_face",)
    @test cell_area.data == data(cf)

    cell_topology = cell_ds.variables["Mesh2"]
    @test cell_topology.attrs["face_coordinates"] == "Mesh2_face_lon Mesh2_face_lat"
    @test haskey(cell_ds.variables, "Mesh2_face_lon")
    @test haskey(cell_ds.variables, "Mesh2_face_lat")
    face_lon = cell_ds.variables["Mesh2_face_lon"]
    face_lat = cell_ds.variables["Mesh2_face_lat"]
    @test face_lon.dims == ("n_face",)
    @test face_lat.dims == ("n_face",)
    @test size(face_lon.data) == (num_cells(g),)
    @test size(face_lat.data) == (num_cells(g),)

    ef = DiscreteField(EdgeLoc, g, edge_values(g), edge_dims(g); name = :edge_flux)
    edge_ds = to_ugrid(ef)
    @test haskey(edge_ds.variables, "edge_flux")
    edge_flux = edge_ds.variables["edge_flux"]
    @test edge_flux.attrs["mesh"] == "Mesh2"
    @test edge_flux.attrs["location"] == "edge"
    @test edge_flux.attrs["coordinates"] == "Mesh2_edge_lon Mesh2_edge_lat"
    @test edge_flux.dims == ("n_edge",)
    @test edge_flux.data == data(ef)

    edge_topology = edge_ds.variables["Mesh2"]
    @test edge_topology.attrs["edge_node_connectivity"] == "Mesh2_edge_nodes"
    @test edge_topology.attrs["edge_coordinates"] == "Mesh2_edge_lon Mesh2_edge_lat"
    @test edge_topology.attrs["edge_dimension"] == "n_edge"
    @test haskey(edge_ds.variables, "Mesh2_edge_nodes")
    edge_nodes_var = edge_ds.variables["Mesh2_edge_nodes"]
    @test edge_nodes_var.dims == ("n_edge", "Two")
    @test edge_nodes_var.attrs["cf_role"] == "edge_node_connectivity"
    @test edge_nodes_var.attrs["start_index"] == 1
    @test size(edge_nodes_var.data) == (num_edges(g), 2)
    for e in 1:num_edges(g)
        @test Tuple(edge_nodes_var.data[e, :]) == edge_nodes(g, e)
    end
    @test haskey(edge_ds.variables, "Mesh2_edge_lon")
    @test haskey(edge_ds.variables, "Mesh2_edge_lat")
    edge_lon = edge_ds.variables["Mesh2_edge_lon"]
    edge_lat = edge_ds.variables["Mesh2_edge_lat"]
    @test edge_lon.dims == ("n_edge",)
    @test edge_lat.dims == ("n_edge",)
    @test size(edge_lon.data) == (num_edges(g),)
    @test size(edge_lat.data) == (num_edges(g),)

    ntf = DiscreteField(NodeLoc, g, node_time_values(g), node_time_dims(g);
        name = :node_temp_by_time)
    nt_ds = to_ugrid(ntf)
    @test nt_ds.variables["node_temp_by_time"].dims == ("n_node", "time")
end

@testset "UGRID own-file in-memory round-trip" begin
    g = small_grid()
    f = DiscreteField(CellLoc, g, cell_values(g), cell_dims(g); name = :cell_area,
        metadata = Dict("units" => "m2"))
    ds = to_ugrid(f)
    mesh_attrs = ds.variables["Mesh2"].attrs

    @test mesh_attrs["manifoldfields_grid_type"] == "LatLonGrid"
    @test mesh_attrs["manifoldfields_lat_edges"] == g.lat_edges
    @test mesh_attrs["manifoldfields_lon_edges"] == g.lon_edges
    @test mesh_attrs["manifoldfields_radius"] == g.R

    loaded = from_ugrid(ds)
    @test loaded isa FieldSet
    @test loaded[:cell_area] isa DiscreteField{CellLoc}
    @test location(loaded[:cell_area]) === CellLoc
    @test data(loaded[:cell_area]) == data(f)
    @test DimensionalData.dims(loaded[:cell_area]) == DimensionalData.dims(f)
    @test DimensionalData.metadata(loaded[:cell_area])["units"] == "m2"
    @test DimensionalData.metadata(loaded[:cell_area])["mesh"] == "Mesh2"
    @test num_cells(mesh(loaded)) == num_cells(g)
    @test num_nodes(mesh(loaded)) == num_nodes(g)

    missing_radius_ds = to_ugrid(g)
    delete!(missing_radius_ds.variables["Mesh2"].attrs, "manifoldfields_radius")
    @test_throws ArgumentError from_ugrid_mesh(missing_radius_ds)

    nf = DiscreteField(NodeLoc, g, node_time_values(g), node_time_dims(g);
        name = :node_temp_by_time)
    loaded_node = from_ugrid(to_ugrid(nf))
    @test loaded_node isa FieldSet
    @test loaded_node[:node_temp_by_time] isa DiscreteField{NodeLoc}
    @test location(loaded_node[:node_temp_by_time]) === NodeLoc
    @test data(loaded_node[:node_temp_by_time]) == data(nf)
    @test DimensionalData.dims(loaded_node[:node_temp_by_time]) ==
          DimensionalData.dims(nf)
end

@testset "UGRID NetCDF save/load round-trip" begin
    g = small_grid()
    f = DiscreteField(CellLoc, g, cell_values(g), cell_dims(g); name = :cell_area,
        metadata = Dict("units" => "m2"))
    path = tempname() * ".nc"

    save_ugrid(f, path)
    loaded = load_ugrid(path)
    @test loaded isa FieldSet
    @test loaded[:cell_area] isa DiscreteField{CellLoc}
    @test location(loaded[:cell_area]) === CellLoc
    @test data(loaded[:cell_area]) == data(f)
    @test DimensionalData.dims(loaded[:cell_area]) == DimensionalData.dims(f)
    @test DimensionalData.metadata(loaded[:cell_area])["units"] == "m2"
    @test DimensionalData.metadata(loaded[:cell_area])["mesh"] == "Mesh2"

    loaded_mesh = load_ugrid_mesh(path)
    @test num_cells(loaded_mesh) == num_cells(g)
    @test num_nodes(loaded_mesh) == num_nodes(g)

    int_values = fill(Int64(3_000_000_000), num_nodes(g))
    int_field = DiscreteField(NodeLoc, g, int_values, node_dims(g); name = :node_count)
    int_path = tempname() * ".nc"
    save_ugrid(int_field, int_path)
    loaded_int = load_ugrid(int_path)
    @test data(loaded_int[:node_count]) == data(int_field)
    @test eltype(data(loaded_int[:node_count])) == Int64

    float32_values = Float32.(node_values(g))
    float32_field = DiscreteField(NodeLoc, g, float32_values, node_dims(g);
        name = :node_temp_float32)
    float32_path = tempname() * ".nc"
    save_ugrid(float32_field, float32_path)
    loaded_float32 = load_ugrid(float32_path)
    @test data(loaded_float32[:node_temp_float32]) == data(float32_field)
    @test eltype(data(loaded_float32[:node_temp_float32])) == Float32

    unnamed_dimarray = DimArray(node_values(g), node_dims(g))
    unnamed_field = DiscreteField(NodeLoc, g, unnamed_dimarray)
    @test DimensionalData.name(unnamed_field) == :field
    unnamed_path = tempname() * ".nc"
    save_ugrid(unnamed_field, unnamed_path)
    unnamed_ds = to_ugrid(unnamed_field)
    @test haskey(unnamed_ds.variables, "field")
    loaded_unnamed = load_ugrid(unnamed_path)
    @test loaded_unnamed isa FieldSet
    @test DimensionalData.name(loaded_unnamed[:field]) == :field
    @test data(loaded_unnamed[:field]) == node_values(g)

    @test_throws ArgumentError save_ugrid(f, tempname() * ".nc"; format = :unknown)
end

@testset "UGRID NetCDF file signature" begin
    g = small_grid()
    f = DiscreteField(NodeLoc, g, node_values(g), node_dims(g); name = :node_temp)
    path = tempname() * ".nc"

    save_ugrid(f, path)
    open(path, "r") do io
        magic = read(io, 4)
        @test magic == UInt8[0x43, 0x44, 0x46, 0x01] ||
              magic == UInt8[0x43, 0x44, 0x46, 0x02] ||
              magic == UInt8[0x89, 0x48, 0x44, 0x46]
    end
end

@testset "UGRID reader rejects invalid datasets" begin
    empty = UGridDataset(
        Dict{String, ManifoldFields.UGridVariable}(),
        Dict{String, Any}("Conventions" => "CF-1.11 UGRID-1.0")
    )
    @test_throws ArgumentError from_ugrid_mesh(empty)

    g = small_grid()
    f = DiscreteField(NodeLoc, g, node_values(g), node_dims(g); name = :node_temp)

    multi_data_ds = to_ugrid(f)
    multi_data_ds.variables["other_node_temp"] = ManifoldFields.UGridVariable(
        node_values(g),
        ("n_node",),
        Dict{String, Any}(
            "mesh" => "Mesh2",
            "location" => "node",
            "coordinates" => "Mesh2_node_lon Mesh2_node_lat"
        )
    )
    loaded_multi = from_ugrid(multi_data_ds)
    @test loaded_multi isa FieldSet
    @test Set(field_names(loaded_multi)) == Set([:node_temp, :other_node_temp])

    missing_location_ds = to_ugrid(f)
    delete!(missing_location_ds.variables["node_temp"].attrs, "location")
    @test_throws ArgumentError from_ugrid(missing_location_ds)

    unsupported_location_ds = to_ugrid(f)
    unsupported_location_ds.variables["node_temp"].attrs["location"] = "volume"
    @test_throws ArgumentError from_ugrid(unsupported_location_ds)

    missing_cf_role_ds = to_ugrid(g)
    delete!(missing_cf_role_ds.variables["Mesh2"].attrs, "cf_role")
    @test_throws ArgumentError from_ugrid_mesh(missing_cf_role_ds)

    wrong_cf_role_ds = to_ugrid(g)
    wrong_cf_role_ds.variables["Mesh2"].attrs["cf_role"] = "not_mesh_topology"
    @test_throws ArgumentError from_ugrid_mesh(wrong_cf_role_ds)

    missing_face_nodes_ds = to_ugrid(g)
    delete!(missing_face_nodes_ds.variables, "Mesh2_face_nodes")
    @test_throws ArgumentError from_ugrid_mesh(missing_face_nodes_ds)

    for attr in (
        "manifoldfields_grid_type",
        "manifoldfields_lat_edges",
        "manifoldfields_lon_edges",
        "manifoldfields_radius"
    )
        missing_metadata_ds = to_ugrid(g)
        delete!(missing_metadata_ds.variables["Mesh2"].attrs, attr)
        @test_throws ArgumentError from_ugrid_mesh(missing_metadata_ds)
    end

    mismatched_rank_ds = to_ugrid(f)
    mismatched_rank_ds.variables["node_temp"] = ManifoldFields.UGridVariable(
        node_values(g),
        ("n_node", "time"),
        copy(mismatched_rank_ds.variables["node_temp"].attrs)
    )
    @test_throws DimensionMismatch from_ugrid(mismatched_rank_ds)
end

@testset "UGRID FieldSet writer" begin
    g = small_grid()
    fs = FieldSet(
        :u => DiscreteField(NodeLoc, g, node_values(g), node_dims(g); name = :u),
        :v => DiscreteField(CellLoc, g, cell_values(g), cell_dims(g);
            name = :v, metadata = Dict("units" => "m2")),
        :w => DiscreteField(EdgeLoc, g, edge_values(g), edge_dims(g); name = :w)
    )

    ds = to_ugrid(fs)
    @test haskey(ds.variables, "Mesh2")
    @test ds.variables["Mesh2"].attrs["face_coordinates"] == "Mesh2_face_lon Mesh2_face_lat"
    @test ds.variables["Mesh2"].attrs["edge_node_connectivity"] == "Mesh2_edge_nodes"

    for (key, locstr) in ((:u, "node"), (:v, "face"), (:w, "edge"))
        @test haskey(ds.variables, String(key))
        var = ds.variables[String(key)]
        @test var.attrs["mesh"] == "Mesh2"
        @test var.attrs["location"] == locstr
    end
    @test ds.variables["v"].attrs["units"] == "m2"
end

@testset "UGRID topology validation" begin
    g = small_grid()
    f = DiscreteField(NodeLoc, g, node_values(g), node_dims(g); name = :node_temp)
    ds = to_ugrid(f)

    # own-file start_index=1 round-trips through validation
    @test ManifoldFields._validate_topology!(g, ds) === g

    # 0-based external variant also validates
    ds0 = to_ugrid(f)
    fn0 = ds0.variables["Mesh2_face_nodes"]
    fn0.data .= fn0.data .- 1
    fn0.attrs["start_index"] = 0
    @test ManifoldFields._validate_topology!(g, ds0) === g

    # missing start_index defaults to 0 (UGRID default)
    ds_noattr = to_ugrid(f)
    fn = ds_noattr.variables["Mesh2_face_nodes"]
    fn.data .= fn.data .- 1
    delete!(fn.attrs, "start_index")
    @test ManifoldFields._validate_topology!(g, ds_noattr) === g

    # perturbed connectivity fails
    ds_bad = to_ugrid(f)
    ds_bad.variables["Mesh2_face_nodes"].data[1, 1] += 1
    @test_throws ArgumentError ManifoldFields._validate_topology!(g, ds_bad)

    # wrong node count fails
    ds_short = to_ugrid(f)
    ds_short.variables["Mesh2_node_lon"].data = ds_short.variables["Mesh2_node_lon"].data[1:(end - 1)]
    @test_throws ArgumentError ManifoldFields._validate_topology!(g, ds_short)
end

@testset "UGRID multi-field round-trip" begin
    g = small_grid()
    fs = sample_fieldset(g)
    ds = to_ugrid(fs)

    loaded = from_ugrid(ds)
    @test loaded isa FieldSet
    @test mesh(loaded) === mesh(loaded[:u])
    @test field_names(loaded) == (:u, :v)
    @test loaded[:u] isa DiscreteField{NodeLoc}
    @test loaded[:v] isa DiscreteField{CellLoc}
    @test data(loaded[:u]) == node_time_values(g)
    @test data(loaded[:v]) == cell_values(g)
    @test DimensionalData.dims(loaded[:u]) == node_time_dims(g)
    @test num_nodes(mesh(loaded)) == num_nodes(g)

    empty_data_ds = to_ugrid(g)
    @test_throws ArgumentError from_ugrid(empty_data_ds)
end

@testset "UGRID external file via injected mesh" begin
    g = small_grid()
    f = DiscreteField(NodeLoc, g, node_values(g), node_dims(g); name = :node_temp)

    ext = to_ugrid(f)
    ma = ext.variables["Mesh2"].attrs
    for k in ("manifoldfields_grid_type", "manifoldfields_lat_edges",
        "manifoldfields_lon_edges", "manifoldfields_radius")
        delete!(ma, k)
    end
    fn = ext.variables["Mesh2_face_nodes"]
    fn.data .= fn.data .- 1
    fn.attrs["start_index"] = 0

    # no mesh source -> error
    @test_throws ArgumentError from_ugrid_mesh(ext)
    @test_throws ArgumentError from_ugrid(ext)

    # injected mesh validates (0-based external connectivity)
    m = from_ugrid_mesh(ext; mesh = g)
    @test m === g
    fs = from_ugrid(ext; mesh = g)
    @test fs isa FieldSet
    @test data(fs[:node_temp]) == node_values(g)

    # mismatched mesh fails validation
    g_perturbed = small_grid(; nlat = 3)
    @test_throws ArgumentError from_ugrid_mesh(ext; mesh = g_perturbed)

    # file-based load_ugrid with mesh injection (kwarg forwarding coverage)
    ext_path = tempname() * ".nc"
    save_ugrid(f, ext_path)
    fs_file = load_ugrid(ext_path; mesh = g)
    @test fs_file isa FieldSet
    @test mesh(fs_file) === g
    @test data(fs_file[:node_temp]) == node_values(g)
    @test_throws ArgumentError load_ugrid(ext_path; mesh = small_grid(; nlat = 3))
end

@testset "UGRID NetCDF multi-field round-trip" begin
    g = small_grid()
    fs = sample_fieldset(g)
    path = tempname() * ".nc"

    save_ugrid(fs, path)
    loaded = load_ugrid(path)
    @test loaded isa FieldSet
    @test mesh(loaded) === mesh(loaded[:u])
    @test data(loaded[:u]) == node_time_values(g)
    @test data(loaded[:v]) == cell_values(g)
    @test DimensionalData.dims(loaded[:u]) == node_time_dims(g)

    # dimension-order variant: (time, node) instead of (node, time)
    fs_tn = FieldSet(
        :u => DiscreteField(NodeLoc, g, time_node_values(g), time_node_dims(g); name = :u),
        :v => DiscreteField(CellLoc, g, cell_values(g), cell_dims(g); name = :v)
    )
    tn_path = tempname() * ".nc"
    save_ugrid(fs_tn, tn_path)
    loaded_tn = load_ugrid(tn_path)
    @test DimensionalData.dims(loaded_tn[:u]) == time_node_dims(g)
    @test data(loaded_tn[:u]) == time_node_values(g)
end

@testset "UGRID own-file mesh metadata writers" begin
    cs = CubedSphereGrid(n = 2)
    ds = to_ugrid(cs)
    a = ds.variables["Mesh2"].attrs
    @test a["manifoldfields_grid_type"] == "CubedSphereGrid"
    @test a["manifoldfields_n"] == 2
    @test a["manifoldfields_projection"] == "gnomonic"
    @test a["manifoldfields_radius"] == 1.0
    @test a["manifoldfields_rotation"] == collect(vec(SMatrix{3, 3, Float64, 9}(I)))

    rg = ReducedGaussianGrid(nlat = 4)
    a = to_ugrid(rg).variables["Mesh2"].attrs
    @test a["manifoldfields_grid_type"] == "ReducedGaussianGrid"
    @test a["manifoldfields_nlat"] == 4
    @test a["manifoldfields_radius"] == 1.0

    hp = HEALPixGrid(nside = 1)
    a = to_ugrid(hp).variables["Mesh2"].attrs
    @test a["manifoldfields_grid_type"] == "HEALPixGrid"
    @test a["manifoldfields_nside"] == 1
    @test a["manifoldfields_ordering"] == "ring"
    @test a["manifoldfields_rotation"] == collect(vec(SMatrix{3, 3, Float64, 9}(I)))
    @test a["manifoldfields_radius"] == 1.0
end

@testset "UGRID own-file round-trips for all grid types" begin
    for g in (CubedSphereGrid(n = 2), ReducedGaussianGrid(nlat = 4), HEALPixGrid(nside = 1))
        f = DiscreteField(NodeLoc, g, collect(Float64, 1:num_nodes(g)),
            (Dim{:node}(1:num_nodes(g)),); name = :psi)
        ds = to_ugrid(f)
        loaded = from_ugrid(ds)
        @test loaded isa FieldSet
        @test typeof(mesh(loaded)) == typeof(g)
        @test num_nodes(mesh(loaded)) == num_nodes(g)
        @test num_cells(mesh(loaded)) == num_cells(g)
        @test data(loaded[:psi]) == data(f)
    end

    # HEALPix nested ordering round-trips too
    hp = HEALPixGrid(nside = 1, ordering = :nested)
    f = DiscreteField(NodeLoc, hp, collect(Float64, 1:num_nodes(hp)),
        (Dim{:node}(1:num_nodes(hp)),); name = :psi)
    loaded = from_ugrid(to_ugrid(f))
    @test typeof(mesh(loaded)) == typeof(hp)
    @test mesh(loaded).ordering == :nested

    # rotated CubedSphereGrid round-trips exactly (rotation is serialized)
    rot = SMatrix{3, 3, Float64, 9}(0.0, 1.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0, 1.0)
    csr = CubedSphereGrid(n = 2, rotation = rot)
    f = DiscreteField(NodeLoc, csr, collect(Float64, 1:num_nodes(csr)),
        (Dim{:node}(1:num_nodes(csr)),); name = :psi)
    loaded = from_ugrid(to_ugrid(f))
    @test typeof(mesh(loaded)) == typeof(csr)
    @test mesh(loaded).rotation == rot
    @test data(loaded[:psi]) == data(f)
    for n in 1:num_nodes(csr)
        @test node_coordinates(mesh(loaded), n) ≈ node_coordinates(csr, n)
    end

    # coordinate mismatch fails loudly
    g0 = ReducedGaussianGrid(nlat = 4)
    f0 = DiscreteField(NodeLoc, g0, collect(Float64, 1:num_nodes(g0)),
        (Dim{:node}(1:num_nodes(g0)),); name = :psi)
    ds0 = to_ugrid(f0)
    ds0.variables["Mesh2_node_lon"].data .+= 1.0
    @test_throws ArgumentError from_ugrid_mesh(ds0)
end

@testset "UGRID FieldSet global-attribute read-back" begin
    g = small_grid()
    u = DiscreteField(NodeLoc, g, node_values(g), node_dims(g); name = :u)
    v = DiscreteField(CellLoc, g, cell_values(g), cell_dims(g); name = :v)
    fs = FieldSet(g, :u => u, :v => v)
    fs = DimensionalData.rebuild(fs; metadata = Dict("title" => "demo", "institution" => "X"))

    loaded = from_ugrid(to_ugrid(fs))
    @test DimensionalData.metadata(loaded)["title"] == "demo"
    @test DimensionalData.metadata(loaded)["institution"] == "X"
    @test !haskey(DimensionalData.metadata(loaded), "Conventions")

    # no metadata written -> reconstructed FieldSet keeps NoMetadata (not an empty Dict)
    plain = from_ugrid(to_ugrid(FieldSet(g, :u => u, :v => v)))
    @test DimensionalData.metadata(plain) isa DimensionalData.NoMetadata
end
