import NCDatasets
import DimensionalData
import ManifoldMeshes:
    CellLoc,
    EdgeLoc,
    NodeLoc,
    cell_centroid,
    cell_nodes,
    edge_nodes,
    node_coordinates,
    num_cells,
    num_edges,
    num_nodes

struct UGridVariable{T}
    data::T
    dims::Tuple{Vararg{String}}
    attrs::Dict{String,Any}
end

struct UGridDataset
    variables::Dict{String,UGridVariable}
    attributes::Dict{String,Any}
end

function _lonlat(p)
    r = sqrt(sum(abs2, p))
    lat = asind(p[3] / r)
    lon = mod(rad2deg(atan(p[2], p[1])), 360)
    return lon, lat
end

function to_ugrid(mesh)
    variables = Dict{String,UGridVariable}()

    variables["Mesh2"] = UGridVariable(
        0,
        (),
        Dict{String,Any}(
            "cf_role" => "mesh_topology",
            "topology_dimension" => 2,
            "node_coordinates" => "Mesh2_node_lon Mesh2_node_lat",
            "face_node_connectivity" => "Mesh2_face_nodes",
            "face_dimension" => "n_face",
        ),
    )

    node_lon = Vector{Float64}(undef, num_nodes(mesh))
    node_lat = Vector{Float64}(undef, num_nodes(mesh))
    for n in 1:num_nodes(mesh)
        node_lon[n], node_lat[n] = _lonlat(node_coordinates(mesh, n))
    end

    variables["Mesh2_node_lon"] = UGridVariable(
        node_lon,
        ("n_node",),
        Dict{String,Any}("standard_name" => "longitude", "units" => "degrees_east"),
    )
    variables["Mesh2_node_lat"] = UGridVariable(
        node_lat,
        ("n_node",),
        Dict{String,Any}("standard_name" => "latitude", "units" => "degrees_north"),
    )

    face_nodes = Matrix{Int}(undef, num_cells(mesh), 4)
    for c in 1:num_cells(mesh)
        face_nodes[c, :] .= cell_nodes(mesh, c)
    end
    variables["Mesh2_face_nodes"] = UGridVariable(
        face_nodes,
        ("n_face", "n_max_face_nodes"),
        Dict{String,Any}("cf_role" => "face_node_connectivity", "start_index" => 1),
    )

    return UGridDataset(
        variables,
        Dict{String,Any}("Conventions" => "CF-1.11 UGRID-1.0"),
    )
end

function to_ugrid(f::DiscreteField)
    ds = to_ugrid(mesh(f))
    Loc = location(f)
    coordinates = _ugrid_coordinates(Loc)
    _add_location_coordinates!(ds, mesh(f), Loc)

    varname = String(Symbol(DimensionalData.name(f)))
    var_dims = _ugrid_field_dims(f, Loc)
    ds.variables[varname] = UGridVariable(
        data(f),
        var_dims,
        Dict{String,Any}(
            "mesh" => "Mesh2",
            "location" => ugrid_location(Loc),
            "coordinates" => coordinates,
        ),
    )
    return ds
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
        Dict{String,Any}("standard_name" => "longitude", "units" => "degrees_east"),
    )
    ds.variables["Mesh2_face_lat"] = UGridVariable(
        lat,
        ("n_face",),
        Dict{String,Any}("standard_name" => "latitude", "units" => "degrees_north"),
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
        Dict{String,Any}("cf_role" => "edge_node_connectivity", "start_index" => 1),
    )
    ds.variables["Mesh2_edge_lon"] = UGridVariable(
        lon,
        ("n_edge",),
        Dict{String,Any}("standard_name" => "longitude", "units" => "degrees_east"),
    )
    ds.variables["Mesh2_edge_lat"] = UGridVariable(
        lat,
        ("n_edge",),
        Dict{String,Any}("standard_name" => "latitude", "units" => "degrees_north"),
    )

    topology_attrs = ds.variables["Mesh2"].attrs
    topology_attrs["edge_node_connectivity"] = "Mesh2_edge_nodes"
    topology_attrs["edge_dimension"] = "n_edge"
    topology_attrs["edge_coordinates"] = "Mesh2_edge_lon Mesh2_edge_lat"
    return nothing
end

function from_ugrid(args...; kwargs...)
    throw(ErrorException("from_ugrid is implemented in a later task"))
end

function from_ugrid_mesh(args...; kwargs...)
    throw(ErrorException("from_ugrid_mesh is implemented in a later task"))
end

function save_ugrid(args...; kwargs...)
    throw(ErrorException("save_ugrid is implemented in a later task"))
end

function load_ugrid(args...; kwargs...)
    throw(ErrorException("load_ugrid is implemented in a later task"))
end

function load_ugrid_mesh(args...; kwargs...)
    throw(ErrorException("load_ugrid_mesh is implemented in a later task"))
end
