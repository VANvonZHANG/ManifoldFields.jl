# ManifoldFields.jl

Discrete fields on `ManifoldMeshes.jl` meshes.

[![CI](https://github.com/VANvonZHANG/ManifoldFields.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/VANvonZHANG/ManifoldFields.jl/actions/workflows/CI.yml)
[![Documentation](https://github.com/VANvonZHANG/ManifoldFields.jl/actions/workflows/Documentation.yml/badge.svg)](https://github.com/VANvonZHANG/ManifoldFields.jl/actions/workflows/Documentation.yml)

```julia
using DimensionalData
using ManifoldFields
using ManifoldMeshes

g = LatLonGrid(
    lat_edges=collect(range(-90.0, 90.0; length=3)),
    lon_edges=collect(range(0.0, 360.0; length=5)),
)

values = collect(Float64, 1:num_nodes(g))
dims = (Dim{:node}(1:num_nodes(g)),)

f = DiscreteField(NodeLoc, g, values, dims; name=:temperature)
temperature_at_point = interpolate(f, -45.0, 45.0)

path = tempname() * ".nc"
save_ugrid(f, path)
loaded = load_ugrid(path)
```

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

`DiscreteField{NodeLoc}`, `DiscreteField{EdgeLoc}`, and `DiscreteField{CellLoc}`
objects remain the single-variable building blocks for node, edge, and cell
variables; `FieldSet` groups them on a shared mesh.

UGRID IO in v0 is focused on own-file NetCDF round trips. UGRID vocabulary is
kept at the IO boundary; field-level location semantics use `ManifoldMeshes`
location tags such as `NodeLoc`, `EdgeLoc`, and `CellLoc`.

## Package Status

This repository is currently treated as a private/development Julia package.
General Registry registration is deferred until the delivered `FieldSet` API and
UGRID multi-field IO contracts have stabilized. See `ROADMAP.md`.
