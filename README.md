# ManifoldFields.jl

Discrete fields on `ManifoldMeshes.jl` meshes.

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

v0 supports separate `DiscreteField{NodeLoc}`, `DiscreteField{EdgeLoc}`, and
`DiscreteField{CellLoc}` objects for node, edge, and cell variables. A
multi-variable `FieldSet` container is deferred.

UGRID IO in v0 is focused on own-file NetCDF round trips. UGRID vocabulary is
kept at the IO boundary; field-level location semantics use `ManifoldMeshes`
location tags such as `NodeLoc`, `EdgeLoc`, and `CellLoc`.
