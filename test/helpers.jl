using DimensionalData
using ManifoldMeshes

small_grid() = LatLonGrid(2, 4)

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
