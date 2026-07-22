import NCDatasets
import ManifoldMeshes: cell_nodes, node_coordinates, num_cells, num_nodes

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
