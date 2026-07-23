import DimensionalData: Dim
import ManifoldMeshes: CellLoc, EdgeLoc, NodeLoc, num_cells, num_edges, num_nodes

location_dimname(::Type{NodeLoc}) = Dim{:node}
location_dimname(::Type{EdgeLoc}) = Dim{:edge}
location_dimname(::Type{CellLoc}) = Dim{:cell}

expected_location_length(::Type{NodeLoc}, mesh) = num_nodes(mesh)
expected_location_length(::Type{EdgeLoc}, mesh) = num_edges(mesh)
expected_location_length(::Type{CellLoc}, mesh) = num_cells(mesh)

ugrid_location(::Type{NodeLoc}) = "node"
ugrid_location(::Type{EdgeLoc}) = "edge"
ugrid_location(::Type{CellLoc}) = "face"

ugrid_dimname(::Type{NodeLoc}) = "n_node"
ugrid_dimname(::Type{EdgeLoc}) = "n_edge"
ugrid_dimname(::Type{CellLoc}) = "n_face"
