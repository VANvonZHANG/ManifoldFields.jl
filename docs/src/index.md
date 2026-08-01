# ManifoldFields.jl

`ManifoldFields.jl` provides discrete field containers for
`ManifoldMeshes.jl` meshes, with location-aware variables and UGRID NetCDF IO.

## Quick Start

```@example quickstart
using DimensionalData
using ManifoldFields
using ManifoldMeshes

mesh = LatLonGrid(
    lat_edges = collect(range(-90.0, 90.0; length = 3)),
    lon_edges = collect(range(0.0, 360.0; length = 5)),
)

values = collect(Float64, 1:num_nodes(mesh))
dims = (Dim{:node}(1:num_nodes(mesh)),)
field = DiscreteField(NodeLoc, mesh, values, dims; name = :temperature)

interpolate(field, -45.0, 45.0)
```

## Scope

The initial package surface is deliberately narrow:

- `DiscreteField{NodeLoc}`, `DiscreteField{EdgeLoc}`, and `DiscreteField{CellLoc}`
- Node-location interpolation
- Own-file UGRID NetCDF round trips
- Broadcast preservation of mesh and location metadata

See the [roadmap](roadmap.md) for the planned FieldSet, UGRID, interpolation,
remapping, and lazy-data work.

