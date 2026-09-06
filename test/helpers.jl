using DimensionalData
using ManifoldMeshes

function small_grid(; nlat = 2, nlon = 4)
    return LatLonGrid(
        lat_edges = collect(range(-90.0, 90.0; length = nlat + 1)),
        lon_edges = collect(range(0.0, 360.0; length = nlon + 1))
    )
end

node_dims(g) = (Dim{:node}(1:num_nodes(g)),)
edge_dims(g) = (Dim{:edge}(1:num_edges(g)),)
cell_dims(g) = (Dim{:cell}(1:num_cells(g)),)
node_time_dims(g) = (Dim{:node}(1:num_nodes(g)), Dim{:time}(1:3))
time_node_dims(g) = (Dim{:time}(1:3), Dim{:node}(1:num_nodes(g)))

node_values(g) = collect(Float64, 1:num_nodes(g))
edge_values(g) = collect(Float64, 1:num_edges(g))
cell_values(g) = collect(Float64, 1:num_cells(g))
node_time_values(g) = reshape(collect(Float64, 1:(num_nodes(g) * 3)), num_nodes(g), 3)
time_node_values(g) = permutedims(node_time_values(g), (2, 1))

function sample_fieldset(g)
    return FieldSet(
        :u => DiscreteField(NodeLoc, g, node_time_values(g), node_time_dims(g); name = :u),
        :v => DiscreteField(CellLoc, g, cell_values(g), cell_dims(g); name = :v)
    )
end
