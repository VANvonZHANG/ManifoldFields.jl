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
- Multi-variable `FieldSet` containers on a shared mesh
- Node-location interpolation
- Own-file UGRID NetCDF round trips
- Broadcast preservation of mesh and location metadata

See the [roadmap](roadmap.md) for the planned interpolation,
remapping, and lazy-data work.

## FieldSet

Multiple `DiscreteField`s sharing one mesh (mixed locations allowed):

```julia
using DimensionalData
using ManifoldFields
using ManifoldMeshes

g = LatLonGrid(lat_edges = collect(-90.0:90.0:90.0), lon_edges = collect(0.0:90.0:360.0))
u = DiscreteField(NodeLoc, g, randn(num_nodes(g)), (Dim{:node}(1:num_nodes(g)),); name = :u)
v = DiscreteField(CellLoc, g, randn(num_cells(g)), (Dim{:cell}(1:num_cells(g)),); name = :v)

fs = FieldSet(g, :u => u, :v => v)   # or FieldSet(:u => u, :v => v)
fs[:u]                               # DiscreteField on the shared mesh
2 .* fs                              # FieldSet (per-field broadcast)
save_ugrid(fs, "fields.nc")
fs2 = load_ugrid("fields.nc")        # FieldSet round-trip
```

See `docs/superpowers/specs/2026-09-04-fieldset-design.md` for the design
rationale.

