import NCDatasets
import DimensionalData
import ManifoldMeshes:
                       CellLoc,
                       CubedSphereGrid,
                       EdgeLoc,
                       Equiangular,
                       Gnomomic,
                       HEALPixGrid,
                       LatLonGrid,
                       NodeLoc,
                       ProjectionStyle,
                       ReducedGaussianGrid,
                       UnstructuredMesh,
                       cell_centroid,
                       cell_nodes,
                       edge_nodes,
                       node_coordinates,
                       num_cells,
                       num_edges,
                       num_nodes
import StaticArrays: SMatrix

mutable struct UGridVariable{T}
    data::T
    dims::Tuple{Vararg{String}}
    attrs::Dict{String, Any}
end

struct UGridDataset
    variables::Dict{String, UGridVariable}
    attributes::Dict{String, Any}
end

function _lonlat(p)
    r = sqrt(sum(abs2, p))
    lat = asind(p[3] / r)
    lon = mod(rad2deg(atan(p[2], p[1])), 360)
    return lon, lat
end

function to_ugrid(mesh)
    variables = Dict{String, UGridVariable}()

    mesh_attrs = Dict{String, Any}(
        "cf_role" => "mesh_topology",
        "topology_dimension" => 2,
        "node_coordinates" => "Mesh2_node_lon Mesh2_node_lat",
        "face_node_connectivity" => "Mesh2_face_nodes",
        "face_dimension" => "n_face"
    )
    _add_own_file_mesh_metadata!(mesh_attrs, mesh)

    variables["Mesh2"] = UGridVariable(0, (), mesh_attrs)

    node_lon = Vector{Float64}(undef, num_nodes(mesh))
    node_lat = Vector{Float64}(undef, num_nodes(mesh))
    for n in 1:num_nodes(mesh)
        node_lon[n], node_lat[n] = _lonlat(node_coordinates(mesh, n))
    end

    variables["Mesh2_node_lon"] = UGridVariable(
        node_lon,
        ("n_node",),
        Dict{String, Any}("standard_name" => "longitude", "units" => "degrees_east")
    )
    variables["Mesh2_node_lat"] = UGridVariable(
        node_lat,
        ("n_node",),
        Dict{String, Any}("standard_name" => "latitude", "units" => "degrees_north")
    )

    ks = [length(cell_nodes(mesh, c)) for c in 1:num_cells(mesh)]
    max_k = maximum(ks)
    fill_value = -1
    face_nodes = fill(fill_value, num_cells(mesh), max_k)
    for c in 1:num_cells(mesh)
        face_nodes[c, 1:ks[c]] .= collect(cell_nodes(mesh, c))
    end
    fn_attrs = Dict{String, Any}(
        "cf_role" => "face_node_connectivity",
        "start_index" => 1
    )
    if !all(ks .== max_k)
        fn_attrs["_FillValue"] = fill_value
    end
    variables["Mesh2_face_nodes"] = UGridVariable(
        face_nodes,
        ("n_face", "n_max_face_nodes"),
        fn_attrs
    )

    return UGridDataset(
        variables,
        Dict{String, Any}("Conventions" => "CF-1.11 UGRID-1.0")
    )
end

_add_own_file_mesh_metadata!(attrs, mesh) = attrs

function _add_own_file_mesh_metadata!(attrs, mesh::LatLonGrid)
    attrs["manifoldfields_grid_type"] = "LatLonGrid"
    attrs["manifoldfields_lat_edges"] = copy(mesh.lat_edges)
    attrs["manifoldfields_lon_edges"] = copy(mesh.lon_edges)
    attrs["manifoldfields_radius"] = mesh.R
    return attrs
end

_projection_name(::Gnomomic) = "gnomonic"
_projection_name(::Equiangular) = "equiangular"

function _add_own_file_mesh_metadata!(attrs, mesh::CubedSphereGrid)
    attrs["manifoldfields_grid_type"] = "CubedSphereGrid"
    attrs["manifoldfields_n"] = mesh.n
    attrs["manifoldfields_projection"] = _projection_name(ProjectionStyle(typeof(mesh)))
    attrs["manifoldfields_rotation"] = collect(vec(mesh.rotation))
    attrs["manifoldfields_radius"] = mesh.R
    return attrs
end

function _add_own_file_mesh_metadata!(attrs, mesh::ReducedGaussianGrid)
    attrs["manifoldfields_grid_type"] = "ReducedGaussianGrid"
    attrs["manifoldfields_nlat"] = mesh.nlat
    attrs["manifoldfields_radius"] = mesh.R
    return attrs
end

function _add_own_file_mesh_metadata!(attrs, mesh::HEALPixGrid)
    attrs["manifoldfields_grid_type"] = "HEALPixGrid"
    attrs["manifoldfields_nside"] = mesh.nside
    attrs["manifoldfields_ordering"] = String(mesh.ordering)
    attrs["manifoldfields_rotation"] = collect(vec(mesh.rotation))
    attrs["manifoldfields_radius"] = mesh.R
    return attrs
end

function _add_own_file_mesh_metadata!(attrs, mesh::UnstructuredMesh)
    attrs["manifoldfields_grid_type"] = "UnstructuredMesh"
    attrs["manifoldfields_radius"] = mesh.R
    return attrs
end

function _add_field_variable!(ds::UGridDataset, f::DiscreteField)
    Loc = location(f)
    varname = String(Symbol(DimensionalData.name(f)))
    ds.variables[varname] = UGridVariable(
        data(f),
        _ugrid_field_dims(f, Loc),
        _data_var_attrs(f, Loc, _ugrid_coordinates(Loc))
    )
    return ds
end

function to_ugrid(f::DiscreteField)
    ds = to_ugrid(mesh(f))
    Loc = location(f)
    _add_location_coordinates!(ds, mesh(f), Loc)
    return _add_field_variable!(ds, f)
end

function _metadata_attrs(metadata)
    attrs = Dict{String, Any}()
    metadata === nothing && return attrs
    metadata isa AbstractDict || return attrs
    for (k, v) in metadata
        attrs[String(k)] = v
    end
    return attrs
end

function _data_var_attrs(f::DiscreteField, Loc, coordinates)
    attrs = _metadata_attrs(DimensionalData.metadata(f))
    attrs["mesh"] = "Mesh2"
    attrs["location"] = ugrid_location(Loc)
    attrs["coordinates"] = coordinates
    return attrs
end

_ugrid_coordinates(::Type{NodeLoc}) = "Mesh2_node_lon Mesh2_node_lat"
_ugrid_coordinates(::Type{CellLoc}) = "Mesh2_face_lon Mesh2_face_lat"
_ugrid_coordinates(::Type{EdgeLoc}) = "Mesh2_edge_lon Mesh2_edge_lat"

function _ugrid_field_dims(f::DiscreteField, Loc)
    locname = DimensionalData.name(location_dimname(Loc))
    return Tuple(
        DimensionalData.name(d) == locname ?
        ugrid_dimname(Loc) :
        String(Symbol(DimensionalData.name(d))) for d in DimensionalData.dims(f)
    )
end

_add_location_coordinates!(::UGridDataset, mesh, ::Type{NodeLoc}) = nothing

function _add_location_coordinates!(ds::UGridDataset, mesh, ::Type{CellLoc})
    lon = Vector{Float64}(undef, num_cells(mesh))
    lat = Vector{Float64}(undef, num_cells(mesh))
    for c in 1:num_cells(mesh)
        lon[c], lat[c] = _lonlat(cell_centroid(mesh, c))
    end

    ds.variables["Mesh2_face_lon"] = UGridVariable(
        lon,
        ("n_face",),
        Dict{String, Any}("standard_name" => "longitude", "units" => "degrees_east")
    )
    ds.variables["Mesh2_face_lat"] = UGridVariable(
        lat,
        ("n_face",),
        Dict{String, Any}("standard_name" => "latitude", "units" => "degrees_north")
    )
    ds.variables["Mesh2"].attrs["face_coordinates"] = "Mesh2_face_lon Mesh2_face_lat"
    return nothing
end

function _add_location_coordinates!(ds::UGridDataset, mesh, ::Type{EdgeLoc})
    edge_node_data = Matrix{Int}(undef, num_edges(mesh), 2)
    lon = Vector{Float64}(undef, num_edges(mesh))
    lat = Vector{Float64}(undef, num_edges(mesh))

    for e in 1:num_edges(mesh)
        n1, n2 = edge_nodes(mesh, e)
        edge_node_data[e, :] .= (n1, n2)
        p1 = node_coordinates(mesh, n1)
        p2 = node_coordinates(mesh, n2)
        midpoint = p1 + p2
        radius = sqrt(sum(abs2, midpoint))
        point = iszero(radius) ? p1 : midpoint / radius
        lon[e], lat[e] = _lonlat(point)
    end

    ds.variables["Mesh2_edge_nodes"] = UGridVariable(
        edge_node_data,
        ("n_edge", "Two"),
        Dict{String, Any}("cf_role" => "edge_node_connectivity", "start_index" => 1)
    )
    ds.variables["Mesh2_edge_lon"] = UGridVariable(
        lon,
        ("n_edge",),
        Dict{String, Any}("standard_name" => "longitude", "units" => "degrees_east")
    )
    ds.variables["Mesh2_edge_lat"] = UGridVariable(
        lat,
        ("n_edge",),
        Dict{String, Any}("standard_name" => "latitude", "units" => "degrees_north")
    )

    topology_attrs = ds.variables["Mesh2"].attrs
    topology_attrs["edge_node_connectivity"] = "Mesh2_edge_nodes"
    topology_attrs["edge_dimension"] = "n_edge"
    topology_attrs["edge_coordinates"] = "Mesh2_edge_lon Mesh2_edge_lat"
    return nothing
end

function _topology_candidates(ds::UGridDataset)
    return [name
            for (name, var) in ds.variables
            if get(var.attrs, "cf_role", "") == "mesh_topology"]
end

function _discover_topology(ds::UGridDataset)
    candidates = _topology_candidates(ds)
    isempty(candidates) && throw(ArgumentError(
        "UGRID dataset has no mesh topology container (cf_role=mesh_topology)"))
    length(candidates) == 1 && return candidates[1]
    "Mesh2" in candidates && return "Mesh2"
    sort!(candidates)
    throw(ArgumentError(
        "UGRID dataset has multiple mesh topology containers ($(join(candidates, ", "))); " *
        "ManifoldFields reads datasets with a single topology or one named Mesh2"))
end

function _node_coordinate_vars(ds::UGridDataset, meshvar::UGridVariable)
    names = String.(split(String(_require_attr(meshvar, "node_coordinates"))))
    length(names) == 2 || throw(ArgumentError(
        "UGRID node_coordinates must name exactly two variables, got '$(join(names, " "))'"))
    v1, v2 = _require_var(ds, names[1]), _require_var(ds, names[2])
    function is_lat(v, name)
        sn = get(v.attrs, "standard_name", "")
        return sn == "latitude" ||
               (sn != "longitude" && occursin("lat", lowercase(name)))
    end
    lat1, lat2 = is_lat(v1, names[1]), is_lat(v2, names[2])
    lat1 == lat2 && throw(ArgumentError(
        "cannot order node_coordinates '$(names[1]) $(names[2])' as longitude/latitude"))
    return lat2 ? (v1, v2) : (v2, v1)
end

function _topology_variables(ds::UGridDataset)
    name = _discover_topology(ds)
    meshvar = ds.variables[name]
    node_lon, node_lat = _node_coordinate_vars(ds, meshvar)
    face_nodes = _require_var(ds, String(_require_attr(meshvar, "face_node_connectivity")))
    face_nodes = _orient_face_nodes(face_nodes, meshvar)
    return (name = name, meshvar = meshvar, node_lon = node_lon, node_lat = node_lat,
        face_nodes = face_nodes)
end

# UGRID makes the face dimension the slowest-varying (first) dimension of
# face_node_connectivity, but some writers store it node-major — e.g. UXarray's
# `to_xarray` emits `(n_max_face_nodes, n_face)` — so orient by dimension name
# and every consumer sees one row per face.
function _orient_face_nodes(face_nodes::UGridVariable, meshvar::UGridVariable)
    face_dim = String(get(meshvar.attrs, "face_dimension", "n_face"))
    if length(face_nodes.dims) == 2 && face_nodes.dims[2] == face_dim
        return UGridVariable(permutedims(face_nodes.data),
            (face_nodes.dims[2], face_nodes.dims[1]), copy(face_nodes.attrs))
    end
    return face_nodes
end

function _topology_aux_names(meshvar::UGridVariable, topology_name::String)
    aux = Set{String}([topology_name])
    for (k, v) in meshvar.attrs
        (endswith(k, "_coordinates") || endswith(k, "_connectivity")) || continue
        for tok in split(String(v))
            push!(aux, String(tok))
        end
    end
    return aux
end

function _find_data_vars(ds::UGridDataset, topology_name::String, aux)
    return sort!([name
                  for (name, var) in ds.variables
                  if get(var.attrs, "mesh", "") == topology_name && !(name in aux) &&
                     !haskey(var.attrs, "cf_role")])
end

function _infer_location(var::UGridVariable, meshvar::UGridVariable)
    loc_dim(attr, default) = String(get(meshvar.attrs, attr, default))
    loc_dim("node_dimension", "n_node") in var.dims && return NodeLoc
    loc_dim("edge_dimension", "n_edge") in var.dims && return EdgeLoc
    loc_dim("face_dimension", "n_face") in var.dims && return CellLoc
    throw(ArgumentError(
        "UGRID data variable $(var.dims) has no location attribute and its dimensions do not identify node/edge/face"))
end

function _require_var(ds::UGridDataset, name::String)
    haskey(ds.variables, name) ||
        throw(ArgumentError("UGRID dataset is missing variable $name"))
    return ds.variables[name]
end

function _require_attr(var::UGridVariable, name::String)
    haskey(var.attrs, name) ||
        throw(ArgumentError("UGRID variable is missing attribute $name"))
    return var.attrs[name]
end

function _loc_from_ugrid(location)
    location == "node" && return NodeLoc
    location == "edge" && return EdgeLoc
    location == "face" && return CellLoc
    throw(ArgumentError("unsupported UGRID location $location"))
end

function _dims_from_ugrid(var::UGridVariable, Loc)
    length(var.dims) == ndims(var.data) ||
        throw(
            DimensionMismatch(
            "UGRID variable dims rank $(length(var.dims)) does not match data ndims $(ndims(var.data))",
        ),
        )
    return Tuple(
        DimensionalData.Dim{Symbol(d == ugrid_dimname(Loc) ?
                                   DimensionalData.name(location_dimname(Loc)) : d)}(
            1:size(var.data, i),
        ) for (i, d) in enumerate(var.dims)
    )
end

function _metadata_vector(attrs, name::String)
    haskey(attrs, name) || throw(ArgumentError("UGRID Mesh2 is missing attribute $name"))
    value = attrs[name]
    value isa AbstractVector && return Vector{Float64}(value)
    throw(ArgumentError("UGRID Mesh2 attribute $name must be a vector"))
end

function _metadata_float(attrs, name::String)
    haskey(attrs, name) || throw(ArgumentError("UGRID Mesh2 is missing attribute $name"))
    value = attrs[name]
    value isa Number && return Float64(value)
    throw(ArgumentError("UGRID Mesh2 attribute $name must be numeric"))
end

function _metadata_int(attrs, name::String)
    haskey(attrs, name) || throw(ArgumentError("UGRID Mesh2 is missing attribute $name"))
    value = attrs[name]
    value isa Integer && return Int(value)
    throw(ArgumentError("UGRID Mesh2 attribute $name must be an integer"))
end

function _metadata_string(attrs, name::String)
    haskey(attrs, name) || throw(ArgumentError("UGRID Mesh2 is missing attribute $name"))
    value = attrs[name]
    value isa AbstractString && return String(value)
    throw(ArgumentError("UGRID Mesh2 attribute $name must be a string"))
end

function _metadata_rotation(attrs)
    v = _metadata_vector(attrs, "manifoldfields_rotation")
    length(v) == 9 || throw(ArgumentError(
        "UGRID Mesh2 attribute manifoldfields_rotation must have 9 entries"))
    return SMatrix{3, 3, Float64, 9}(reshape(v, 3, 3))
end

function _validate_topology!(m, ds::UGridDataset)
    topo = _topology_variables(ds)
    n_node = length(topo.node_lon.data)
    n_face = size(topo.face_nodes.data, 1)
    num_nodes(m) == n_node || throw(ArgumentError(
        "mesh node count $(num_nodes(m)) != UGRID node count $n_node"))
    num_cells(m) == n_face || throw(ArgumentError(
        "mesh cell count $(num_cells(m)) != UGRID face count $n_face"))
    start_index = Int(get(topo.face_nodes.attrs, "start_index", 0))
    fillv = haskey(topo.face_nodes.attrs, "_FillValue") ?
            Int(topo.face_nodes.attrs["_FillValue"]) : nothing
    for c in 1:num_cells(m)
        expected = collect(Int, cell_nodes(m, c))
        row = Int.(collect(topo.face_nodes.data[c, :]))
        active = fillv === nothing ? row : row[row .!= fillv]
        actual = active .+ (1 - start_index)
        actual == expected || throw(ArgumentError(
            "UGRID face_node_connectivity row $c does not match mesh connectivity (start_index-normalized)"))
    end
    return m
end

function _validate_node_coordinates!(m, ds::UGridDataset; atol = 1e-8)
    topo = _topology_variables(ds)
    for n in 1:num_nodes(m)
        lon, lat = _lonlat(node_coordinates(m, n))
        abs(mod(lon - topo.node_lon.data[n] + 180.0, 360.0) - 180.0) <= atol ||
            throw(ArgumentError(
                "reconstructed mesh node $n longitude does not match the UGRID file; the file describes a different geometry — pass mesh= explicitly"))
        abs(lat - topo.node_lat.data[n]) <= atol ||
            throw(ArgumentError(
                "reconstructed mesh node $n latitude does not match the UGRID file; the file describes a different geometry — pass mesh= explicitly"))
    end
    return m
end

function _unstructured_from_topology(topo, attrs)
    fn = topo.face_nodes
    fillv = haskey(fn.attrs, "_FillValue") ? Int(fn.attrs["_FillValue"]) : nothing
    ks = [fillv === nothing ? size(fn.data, 2) :
          count(!=(fillv), Int.(fn.data[c, :])) for c in 1:size(fn.data, 1)]
    all(3 .<= ks) || throw(ArgumentError(
        "face_node_connectivity of topology '$(topo.name)' has cells with fewer than 3 nodes"))
    start_index = Int(get(fn.attrs, "start_index", 0))
    conn = Matrix{Int}(fn.data)
    for c in 1:size(conn, 1)
        for k in 1:ks[c]
            conn[c, k] += (1 - start_index)
        end
    end
    R = haskey(attrs, "manifoldfields_radius") ?
        Float64(attrs["manifoldfields_radius"]) : 1.0
    return UnstructuredMesh(collect(Float64, topo.node_lon.data),
        collect(Float64, topo.node_lat.data), conn; R = R, start_index = 1,
        fill_value = fillv === nothing ? -1 : fillv)
end

function from_ugrid_mesh(ds::UGridDataset; grid_type = nothing, mesh = nothing)
    topo = _topology_variables(ds)
    meshvar = topo.meshvar
    _require_attr(meshvar, "topology_dimension") == 2 ||
        throw(ArgumentError("mesh topology $(topo.name) must have topology_dimension=2"))
    _require_attr(meshvar, "face_dimension")

    if mesh !== nothing
        return _validate_topology!(mesh, ds)
    end

    attrs = meshvar.attrs
    if !haskey(attrs, "manifoldfields_grid_type")
        # Foreign topology: reconstruct generically as an UnstructuredMesh when
        # the connectivity forms valid cells; otherwise stay actionable.
        m = try
            _unstructured_from_topology(topo, attrs)
        catch err
            err isa ArgumentError || rethrow()
            throw(ArgumentError(
                "cannot reconstruct a ManifoldMeshes grid from topology '$(topo.name)': it carries no manifoldfields_* grid metadata and its connectivity cannot form an UnstructuredMesh ($(err.msg)); pass mesh= explicitly if the geometry matches a supported grid type"))
        end
        return _validate_node_coordinates!(m, ds)
    end
    requested_grid_type = grid_type === nothing ? attrs["manifoldfields_grid_type"] :
                          grid_type
    if requested_grid_type == "LatLonGrid"
        lat_edges = _metadata_vector(attrs, "manifoldfields_lat_edges")
        lon_edges = _metadata_vector(attrs, "manifoldfields_lon_edges")
        radius = _metadata_float(attrs, "manifoldfields_radius")
        mesh = LatLonGrid(lat_edges = lat_edges, lon_edges = lon_edges; R = radius)
        return _validate_node_coordinates!(_validate_topology!(mesh, ds), ds)
    elseif requested_grid_type == "CubedSphereGrid"
        n = _metadata_int(attrs, "manifoldfields_n")
        projection = Symbol(_metadata_string(attrs, "manifoldfields_projection"))
        rotation = _metadata_rotation(attrs)
        radius = _metadata_float(attrs, "manifoldfields_radius")
        mesh = CubedSphereGrid(n = n, projection = projection, rotation = rotation; R = radius)
        return _validate_node_coordinates!(_validate_topology!(mesh, ds), ds)
    elseif requested_grid_type == "ReducedGaussianGrid"
        nlat = _metadata_int(attrs, "manifoldfields_nlat")
        radius = _metadata_float(attrs, "manifoldfields_radius")
        mesh = ReducedGaussianGrid(nlat = nlat; R = radius)
        return _validate_node_coordinates!(_validate_topology!(mesh, ds), ds)
    elseif requested_grid_type == "HEALPixGrid"
        nside = _metadata_int(attrs, "manifoldfields_nside")
        ordering = Symbol(_metadata_string(attrs, "manifoldfields_ordering"))
        rotation = _metadata_rotation(attrs)
        radius = _metadata_float(attrs, "manifoldfields_radius")
        mesh = HEALPixGrid(nside = nside, ordering = ordering, rotation = rotation; R = radius)
        return _validate_node_coordinates!(_validate_topology!(mesh, ds), ds)
    elseif requested_grid_type == "UnstructuredMesh"
        m = _unstructured_from_topology(topo, attrs)
        return _validate_node_coordinates!(m, ds)
    end

    throw(ArgumentError("cannot reconstruct mesh without supported ManifoldFields mesh metadata"))
end

function from_ugrid(ds::UGridDataset; grid_type = nothing, mesh = nothing)
    topology_name = _discover_topology(ds)
    meshvar = ds.variables[topology_name]
    data_vars = _find_data_vars(ds, topology_name, _topology_aux_names(meshvar, topology_name))
    isempty(data_vars) && throw(ArgumentError(
        "UGRID dataset has no data variables (variables with a mesh attribute referencing topology '$topology_name')"))
    m = mesh === nothing ? from_ugrid_mesh(ds; grid_type = grid_type) :
        _validate_topology!(mesh, ds)
    # A reconstructed UnstructuredMesh derives its edge numbering from the
    # face-node table, not the file's edge ordering (edge_node_connectivity is
    # never read), so edge data would silently bind to wrong edges.
    foreign = mesh === nothing && !haskey(meshvar.attrs, "manifoldfields_grid_type")
    fields_nt = NamedTuple{Tuple(Symbol.(data_vars))}(map(data_vars) do varname
        var = ds.variables[varname]
        Loc = haskey(var.attrs, "location") ?
              _loc_from_ugrid(_require_attr(var, "location")) :
              _infer_location(var, meshvar)
        Loc === EdgeLoc && foreign && throw(ArgumentError(
            "foreign topology '$(topology_name)': edge data variables cannot be attached to a reconstructed UnstructuredMesh (its edge numbering is derived from the face-node table, not the file's edge ordering); pass mesh= explicitly or read edge variables separately"))
        return DiscreteField(Loc, m, var.data, _dims_from_ugrid(var, Loc);
            name = Symbol(varname), metadata = var.attrs)
    end)
    fs = FieldSet(m, fields_nt)
    attrs = copy(ds.attributes)
    delete!(attrs, "Conventions")
    isempty(attrs) && return fs
    return DimensionalData.rebuild(fs; metadata = attrs)
end

function _define_dimensions!(nc, ds::UGridDataset)
    lengths = Dict{String, Int}()
    for var in values(ds.variables)
        for (i, d) in enumerate(var.dims)
            len = size(var.data, i)
            if haskey(lengths, d)
                lengths[d] == len ||
                    throw(DimensionMismatch("dimension $d has inconsistent lengths $(lengths[d]) and $len"))
            else
                lengths[d] = len
            end
        end
    end
    for (name, len) in lengths
        NCDatasets.defDim(nc, name, len)
    end
    return nc
end

function _field_nc_type(data)
    data isa Number && return typeof(data)
    return eltype(data)
end

function _nc_type(var::UGridVariable)
    haskey(var.attrs, "mesh") && return _field_nc_type(var.data)
    data = var.data
    data isa Integer && return Int32
    data isa AbstractFloat && return Float64
    eltype(data) <: Integer && return Int32
    eltype(data) <: AbstractFloat && return Float64
    return eltype(data)
end

function _write_attrs!(ncvar, attrs)
    for (k, v) in attrs
        ncvar.attrib[k] = v
    end
    return ncvar
end

function _write_one_variable!(nc, name::String, var::UGridVariable)
    T = _nc_type(var)
    ncvar = NCDatasets.defVar(nc, name, T, var.dims)
    # NetCDF requires _FillValue to be defined before any data is written
    # (error -122) and with the variable's on-disk type (error -45). All
    # other attributes are written after the data, preserving the on-disk
    # layout of fill-free files byte-for-byte.
    attrs = var.attrs
    if haskey(attrs, "_FillValue")
        ncvar.attrib["_FillValue"] = T(attrs["_FillValue"])
        attrs = copy(attrs)
        delete!(attrs, "_FillValue")
    end
    if var.dims == ()
        ncvar[] = var.data
    else
        ncvar[:] = var.data
    end
    _write_attrs!(ncvar, attrs)
    return nc
end

function _write_variables!(nc, ds::UGridDataset)
    for (name, var) in ds.variables
        _write_one_variable!(nc, name, var)
    end
    return nc
end

function save_ugrid(x, path::AbstractString; format = :netcdf)
    format == :netcdf || throw(ArgumentError("v0 only supports format=:netcdf"))
    ds = to_ugrid(x)
    NCDatasets.NCDataset(path, "c") do nc
        for (k, v) in ds.attributes
            nc.attrib[k] = v
        end
        _define_dimensions!(nc, ds)
        _write_variables!(nc, ds)
    end
    return nothing
end

function _read_var_data(nc, name)
    # Connectivity variables are read raw: `_FillValue` entries stay as their
    # on-disk value instead of being masked to `missing` by the CF layer, so
    # the fill-aware topology readers see the fill value itself. Everything
    # else is CF-decoded (`_FillValue` per CF, scale/offset unpacked).
    raw = occursin("face_node_connectivity",
        String(get(nc[name].attrib, "cf_role", "")))
    v = raw ? NCDatasets.variable(nc, name) : nc[name]
    names = NCDatasets.dimnames(v)
    indices = ntuple(_ -> Colon(), length(names))
    return isempty(indices) ? v[] : v[indices...]
end

function _read_ugrid_dataset(path::AbstractString)
    vars = Dict{String, UGridVariable}()
    attrs = Dict{String, Any}()
    NCDatasets.NCDataset(path, "r") do nc
        for (k, v) in nc.attrib
            attrs[String(k)] = v
        end
        for name in keys(nc)
            raw = NCDatasets.variable(nc, String(name))
            vars[String(name)] = UGridVariable(
                _read_var_data(nc, String(name)),
                Tuple(String.(NCDatasets.dimnames(raw))),
                Dict{String, Any}(String(k) => val for (k, val) in raw.attrib)
            )
        end
    end
    return UGridDataset(vars, attrs)
end

"""
    load_ugrid(path; grid_type = nothing, mesh = nothing) -> FieldSet

Read a UGRID file as a `FieldSet`. `_FillValue` on data variables is applied
per CF (fill entries surface as `missing`); the face-node connectivity variable
is read raw so fill values stay integers. Topologies carrying `manifoldfields_*`
metadata reconstruct their original grid; foreign topologies reconstruct as
an `UnstructuredMesh` (at unit radius when `manifoldfields_radius` is absent).
Passing `mesh=` bypasses reconstruction.
"""
function load_ugrid(path::AbstractString; grid_type = nothing, mesh = nothing)
    return from_ugrid(_read_ugrid_dataset(path); grid_type = grid_type, mesh = mesh)
end

function load_ugrid_mesh(path::AbstractString; grid_type = nothing, mesh = nothing)
    from_ugrid_mesh(_read_ugrid_dataset(path); grid_type = grid_type, mesh = mesh)
end

function to_ugrid(fs::FieldSet)
    ds = to_ugrid(mesh(fs))
    fs_fields = fields(fs)
    locs = unique(location(f) for f in Tuple(fs_fields))
    for Loc in locs
        _add_location_coordinates!(ds, mesh(fs), Loc)
    end
    for f in Tuple(fs_fields)
        _add_field_variable!(ds, f)
    end
    for (k, v) in _metadata_attrs(DimensionalData.metadata(fs))
        ds.attributes[String(k)] = v
    end
    return ds
end
