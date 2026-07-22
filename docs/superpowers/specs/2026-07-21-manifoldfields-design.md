# ManifoldFields.jl — Discrete Fields on Manifold Meshes — Design

- **Date:** 2026-07-21
- **Status:** Revised after design review
- **Scope:** `ManifoldFields.jl` v0 — the field container that pairs a `ManifoldMeshes`
  mesh with a `YAXArrays` data cube, with element-wise arithmetic, point-gather
  interpolation, and UGRID NetCDF own-file round-trip IO.
- **Supersedes:** none. Named consumer of the locate/interpolation primitives
  (locate spec §4b / §11); gives `AbstractLocation` its first real use.

## 1. Motivation

`ManifoldMeshes.jl` answers "what does the mesh look like?" and, since v0.5/v0.6,
"which cell contains this point, and with what bilinear weights?". The missing
piece is the **field**: the data that lives on the mesh. `ManifoldFields.jl` is
the thin coupling that pairs a mesh with a data array.

It is the Julia analogue of Python's **UXarray** `UxDataArray`, built on
**YAXArrays.jl** (the Julia xarray) instead of xarray. It is the named consumer
of the locate/interpolation primitives, the bridge to NetCDF/Zarr cube storage,
and the prerequisite for a future `ManifoldRegrid.jl`. It also gives
`AbstractLocation` (`NodeLoc` / `CellLoc` / `EdgeLoc`, defined in
`ManifoldMeshes` `traits.jl` but currently *unused*) its first real consumer.

## 2. Scope

### In scope (v0)

- `DiscreteField{Loc} <: AbstractDimArray` — a mesh + a `YAXArrays` data cube
- Constructor with location-axis dimension check (mismatch fails at construction)
- Accessors (`mesh`, `data`, `location`) plus the operations inherited from
  `DimensionalData` (dims, indexing, reductions, non-spatial selectors)
- Element-wise arithmetic via inherited DD broadcast, with mesh-identity
  enforcement
- `interpolate(f, lat, lon)` point-gather (+ `SVector{3}` overload), `NodeLoc`
  only, multi-dimensional (reduces over the location axis, preserves all
  non-location axes)
- NetCDF UGRID IO: `to_ugrid` / `from_ugrid` / `save_ugrid` / `load_ugrid`
  (write topology + one data variable; read our own single-field files)

### Explicitly out of scope (v0)

- **`FieldSet`** (the `UxDataset` analogue — one mesh + named fields) — v0.1
- **`CellLoc` / `EdgeLoc` gather** — `interpolation_weights` is `NodeLoc`-only
  today; other locations carry more design choices and no consumer
- **Conservative remap, analysis operators** (`gradient`, `integrate`) — these
  live in `ManifoldRegrid.jl` / a future operators layer
- **Generic external-UGRID read of arbitrary unstructured meshes** — blocked on
  the `IsMesh` concrete type (locate spec §11-4)
- **Zarr UGRID conformance** — deferred until the NetCDF contract is solid
- **Multi-variable UGRID datasets** — require `FieldSet` (v0.1)
- **In-place ops (`setindex!`)** — a field is not an object written through
  (footgun for disk-backed cubes)

## 3. Architectural Context

```
ManifoldRegrid.jl      ── bilinear + conservative remap orchestration (future)
       │  consumes
ManifoldFields.jl      ── DiscreteField{Loc} = (mesh, YAXArray) + UGRID IO   ← this package
       ├── ManifoldMeshes.jl   locate_cell, interpolation_weights, num_*,
       │                        cell_nodes, node_coordinates (pure geometry)
       └── YAXArrays.jl / DimensionalData   Cube (an AbstractDimArray), axes, IO
```

### Reference designs (studied, not copied)

- **UXarray** ([data structures](https://uxarray.readthedocs.io/en/stable/user-guide/data-structures.html)):
  `Grid` (pure topology) / `UxDataArray = DataArray + .uxgrid` / `UxDataset`
  (one shared grid, many variables). We mirror this split as
  `AbstractManifoldMesh` / `DiscreteField` / `FieldSet`.
- **UGRID conventions** ([spec](https://ugrid-conventions.github.io/ugrid-conventions/)):
  the interchange format. We borrow its **structure** — topology/data
  separation, staggered `location`, `face_node_connectivity` conventions,
  shared-topology multi-variable model — **not its vocabulary**.

### YAXArrays reality (verified from source)

`YAXArray` (a.k.a. `Cube`) is declared `struct YAXArray{...} <:
AbstractDimArray{...}` — YAXArrays is built on `DimensionalData`. It is a
concrete (final) struct, so it **cannot be subtyped**. The Julia-honest
equivalent of UXarray's `UxDataArray <: xr.DataArray` is therefore to subtype
the **abstract** base: `DiscreteField <: AbstractDimArray`, a sibling of
`YAXArray`, inheriting DimensionalData's broadcast / dims / selectors /
reductions.

### Layering rules (non-negotiable)

- **`ManifoldMeshes` stays pure geometry** — no field data, no
  NetCDF/UGRID/YAXArrays dependency. It never touches `f.data`.
- **`ManifoldFields` owns the field + UGRID IO** — it calls mesh primitives and
  YAXArrays; it never re-implements geometry.
- **One vocabulary throughout the Julia stack: `cell` / `node` / `edge`**
  (matches `ManifoldMeshes`). UGRID's `face` appears **only** inside
  `to_ugrid` / `from_ugrid`. This reaffirms locate spec §3 ("internal names
  stay as `cell`; the `cell ↔ face` mapping lives only at the serialization
  boundary").

### Dependencies

`ManifoldMeshes` (the mesh + locate primitives), `YAXArrays` (the Cube substrate
+ IO, which brings `DimensionalData` transitively — the layer we subtype), and
`StaticArrays` (only for the `interpolate(f, p::SVector{3})` overload). No
`Makie`, no direct `Manifolds` / `ManifoldsBase`. Compat versions are pinned
during implementation.

## 4. Public Interface (summary)

```julia
# Construction
DiscreteField(::Type{Loc}, mesh, data::AbstractDimArray; name, metadata) where {Loc<:AbstractLocation}
DiscreteField(::Type{Loc}, mesh, values::AbstractArray, dims; name, metadata) where {Loc<:AbstractLocation}

# Example location tags
DiscreteField(NodeLoc, mesh, node_data, node_dims)
DiscreteField(EdgeLoc, mesh, edge_data, edge_dims)
DiscreteField(CellLoc, mesh, cell_data, cell_dims)

# Accessors
mesh(f::DiscreteField)          -> AbstractManifoldMesh
data(f::DiscreteField)          -> AbstractArray     # the underlying cube/array
location(f::DiscreteField)      -> Type{<:AbstractLocation}
withmesh(f::DiscreteField, mesh) -> DiscreteField     # explicit rebinding

# Inherited from AbstractDimArray: dims, size, f[...], non-spatial selectors,
# reductions, and element-wise broadcast. Spatial queries use interpolate(...).

# Point gather (NodeLoc only in v0)
interpolate(f::DiscreteField, lat::Real, lon::Real)
interpolate(f::DiscreteField, p::StaticArrays.SVector{3})

# UGRID IO
to_ugrid(mesh::AbstractManifoldMesh)                    -> UGridDataset
to_ugrid(field::DiscreteField)                          -> UGridDataset
from_ugrid(ds::UGridDataset; grid_type=nothing)         -> DiscreteField
from_ugrid_mesh(ds::UGridDataset; grid_type=nothing)    -> AbstractManifoldMesh
save_ugrid(x, path::AbstractString; format=:netcdf)     -> nothing
load_ugrid(path::AbstractString; grid_type=nothing)     -> DiscreteField
load_ugrid_mesh(path::AbstractString; grid_type=nothing) -> AbstractManifoldMesh
```

## 5. Core Type

```julia
using DimensionalData: AbstractDimArray, Dimension, Dim
using ManifoldMeshes: AbstractLocation, AbstractManifoldMesh

struct DiscreteField{
    Loc<:AbstractLocation,
    T,N,D,A<:AbstractArray{T,N},
    M<:AbstractManifoldMesh,
} <: AbstractDimArray{T,N,D,A}
    data::A            # the cube / array (AbstractDimArray `data` field)
    dims::D            # DD Dimensions, one of which is the location axis
    refdims
    name::Symbol
    metadata
    mesh::M            # the manifold mesh this field lives on
end
```

- **`Loc`** is a type parameter (`NodeLoc` / `CellLoc` / `EdgeLoc`, imported
  from `ManifoldMeshes`) — the in-Julia equivalent of UGRID's `location`
  attribute, promoted to a type. The public constructor takes the location as a
  **type tag** (`DiscreteField(NodeLoc, ...)`), not an instance tag
  (`DiscreteField(NodeLoc(), ...)`). `Loc` selects the location axis and drives
  the dimension check.
- **Location axis** is a `DimensionalData` `Dimension` named `:node` / `:cell` /
  `:edge` — chosen by `Loc`, so the type and the dim agree (`CellLoc ↔ :cell`).
  It may appear in any array position, but must appear exactly once.
- **Non-location axes** are free (`:time`, `:layer`, …) and follow CF
  conventions; multi-dimensional fields `(loc, time, layer)` and
  `(time, loc, layer)` are first-class, not deferred.
- **`data`** is any `AbstractArray`; the blessed substrate is a `YAXArray`
  (`Cube`), so NetCDF/Zarr-backed out-of-core cubes plug in with no adapter.

### Location trait table (defined in `ManifoldFields`)

This is the single, centralized dispatch table — `ManifoldMeshes` stays
pure-geometry-vocab; `ManifoldFields` owns field/IO location semantics. This is
plain Julia tag dispatch, not a heavyweight trait system:

```julia
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
```

An internal `LocationKind` hierarchy (`PointLocation` / `LineLocation` /
`AreaLocation`) is intentionally deferred. v0 does not need dispatch over those
coarser geometric classes; adding it now would be structure without a consumer.

### Constructor — dimension check at construction time

```julia
function DiscreteField(::Type{Loc}, mesh::M, data::AbstractDimArray; name=:field, metadata=nothing) where {Loc<:AbstractLocation,M}
    locdim = location_dimname(Loc)
    expected = expected_location_length(Loc, mesh)
    # reuse dims(data); require exactly one matching location axis of expected length
    …
    DiscreteField{Loc,…}(data, dims, …, mesh)
end

function DiscreteField(::Type{Loc}, mesh::M, values::AbstractArray, dims; name=:field, metadata=nothing) where {Loc<:AbstractLocation,M}
    locdim = location_dimname(Loc)
    expected = expected_location_length(Loc, mesh)
    # ordinary arrays do not carry names; caller must provide dims explicitly
    …
    DiscreteField{Loc,…}(values, dims, …, mesh)
end
```

The whole point of carrying `Loc`: a field/grid mismatch fails at
**construction**, not at first use.

Constructor rules:

- The location argument is a subtype object such as `NodeLoc`, not an instance
  such as `NodeLoc()`.
- Plain `AbstractArray` input without explicit `dims` is an `ArgumentError`.
- The location dimension must be present exactly once.
- The location dimension length must equal `num_nodes`, `num_cells`, or
  `num_edges`, selected by `Loc`.
- The location dimension is looked up by dimension name/type, not by array
  position.

### Required `AbstractDimArray` compatibility spike

Before committing the v0 type shape, implementation must validate the exact
`DimensionalData` contract against the pinned `DimensionalData` / `YAXArrays`
versions. The spike must prove:

- `dims(f)`, `size(f)`, and basic indexing work.
- `DimensionalData.data(f)`, `refdims(f)`, `name(f)`, and `metadata(f)` return
  the intended values.
- `sum(f)`, `mean(f)`, and `extrema(f)` work through DD.
- `sin.(f)` and scalar broadcast return `DiscreteField` while preserving mesh,
  location, dims, name, and metadata according to §7.
- `f .+ g` triggers the mesh/`Loc` checks in §7.

If this spike shows that the custom `AbstractDimArray` subtype requires a
broader or unstable hook surface, the fallback design is the wrapper approach
(F2 in §12), revisited before implementation proceeds.

## 6. Accessors & Inherited Behaviour

- `mesh(f)`, `data(f)`, `location(f)` (returns the `Loc` **type**).
- **Inherited from `AbstractDimArray`** (this is the payoff of F1): `dims(f)`,
  `size(f)`, indexing `f[...]`, reductions `sum` / `mean` / `extrema`,
  selectors over non-spatial axes such as `time`, and element-wise broadcast.
  Spatial lookup is not a DD selector; it goes through `interpolate`.
- **No `setindex!`.** A `DiscreteField` is not written through (disk-cube
  safety). Mutation, when needed, is via `data(f)` at the caller's risk.

## 7. Arithmetic (inherited from DimensionalData)

Because `DiscreteField <: AbstractDimArray`, DD's broadcast machinery preserves
type + dims through element-wise ops: `f .+ g`, `2 .* f`, `sin.(f)`, `f .* g`
all return a `DiscreteField` with the same axes.

- **Single custom hook** — override the result-build (`similar` / broadcast
  `copy`) so that, across all `DiscreteField` args: (a) meshes are identical
  (`===`), else `DimensionMismatch`; (b) `Loc` matches, else
  `DimensionMismatch`; (c) the result carries the shared mesh + `Loc`; (d)
  dims, name, metadata, and refdims preservation is explicit and covered by
  tests. Scalars and plain numbers combine naturally (no mesh check against
  them).
- **Mesh-identity rule:** `f.mesh === g.mesh` (object identity), not structural
  `==`. O(1), and it catches the real bug (fields on different grids). Users
  who reconstruct "the same" grid share the object.
- **Explicit mesh rebinding:** `withmesh(f::DiscreteField, mesh)` returns a new
  `DiscreteField` over the same data/dims/name/metadata with a different mesh,
  after rerunning the location-length constructor check. It performs no
  structural-equivalence proof; callers use it when they know an IO round-trip
  or reconstruction produced the same grid but not the same object.
- **Operator forms:** Julia idiom is dotted broadcast (`f .+ g`). Scalar `*` /
  `/` work natively (`*(::Number, ::AbstractArray)`). Non-dotted binary `+` /
  `-` may be added as thin aliases to broadcast for ergonomics (optional).

## 8. Point Gather

```julia
interpolate(f::DiscreteField, lat::Real, lon::Real)
interpolate(f::DiscreteField, p::StaticArrays.SVector{3})   # converts to (lat,lon), calls primary
```

Composes the locate primitives (it adds **no** geometry):

```julia
cid    = locate_cell(mesh(f), lat, lon)
nodes, w = interpolation_weights(mesh(f), cid, lat, lon)   # NTuple{4,Int}, NTuple{4,Float64}
# index the location axis at `nodes`, weight by `w`, sum over the location axis
```

- **Multi-dimensional semantics:** the gather reduces over the **location axis**
  wherever it appears and preserves all non-location axes. A `NodeLoc` field on
  `(node,)` returns a scalar; fields on `(node, time)` and `(time, node)` both
  return a length-`T` vector (the interpolated time series at the query point).
- **`NodeLoc` only in v0** — `interpolation_weights` is `NodeLoc`-only today.
  Calling `interpolate` on a `CellLoc` / `EdgeLoc` field errors with an
  informative message (deferred, see §12).
- **Error propagation:** out-of-domain coordinates reuse `locate_cell`'s
  `ArgumentError`; `missing`/`NaN` values in the gathered data follow ordinary
  Julia arithmetic propagation for the underlying array element type.
- This is the direct analogue of UXarray's `.interp(lon_lat=(…), method=)`.

## 9. UGRID NetCDF IO (v0)

We borrow UGRID's **structure** and confine the `cell ↔ face` vocabulary
translation to these functions. `ManifoldMeshes` stays pure (no NetCDF/UGRID
dependency); the IO lives in `ManifoldFields`, which calls mesh primitives and
YAXArrays for the bytes. UGRID location strings and dimension names are obtained
through the location trait table (`ugrid_location`, `ugrid_dimname`), so IO does
not need `if loc == :node` branches.

The v0 writer targets **UGRID v1.0-conformant NetCDF** for 2D meshes, not merely
"UGRID-inspired" output. The compliance target is: satisfy UGRID required
attributes/dimensions for mesh variables, mesh coordinate variables, mesh
connectivity variables, and mesh data variables; additionally emit the selected
CF/UGRID recommendations below where they are cheap and deterministic.

```julia
# in-memory conversion (DiscreteField ⇄ UGRID dataset)
to_ugrid(mesh::AbstractManifoldMesh)   -> UGridDataset   # topology only
to_ugrid(field::DiscreteField)         -> UGridDataset   # topology + one data var
from_ugrid(ds::UGridDataset; grid_type=nothing) -> DiscreteField
from_ugrid_mesh(ds::UGridDataset; grid_type=nothing) -> AbstractManifoldMesh

# file IO (bytes delegated to YAXArrays; ManifoldFields adds UGRID scaffolding)
save_ugrid(x, path; format=:netcdf)
load_ugrid(path; grid_type=nothing) = from_ugrid(YAXArrays.read(path))
load_ugrid_mesh(path; grid_type=nothing) = from_ugrid_mesh(YAXArrays.read(path))
```

### `UGridDataset` and return shapes

`UGridDataset` is the in-memory UGRID representation: a `YAXArrays`-backed
dataset holding the scalar mesh topology variable (`Mesh2`), the coordinate and
connectivity cubes (`Mesh2_node_lon` / `Mesh2_node_lat`,
`Mesh2_face_nodes`, and conditionally `Mesh2_edge_nodes`), and the data cube
carrying `mesh` / `location` / `coordinates` attributes. Its exact Julia type
is a plan-phase detail (a `YAXArrays.Dataset` or a thin wrapper).

`from_ugrid` returns a single `DiscreteField`, because the field already carries
its mesh. `from_ugrid_mesh` is the topology-only escape hatch. In v0 (no
`FieldSet`) multi-data-variable datasets raise `ArgumentError("multi-variable
UGRID datasets require FieldSet")`.

### Vocabulary translation (lives *only* in `to_ugrid` / `from_ugrid`)

| Internal (`ManifoldMeshes` vocab) | UGRID NetCDF | Direction |
|---|---|---|
| `NodeLoc` / `:node` axis | `location="node"` / `n_node` dim | both |
| `EdgeLoc` / `:edge` axis | `location="edge"` / `n_edge` dim | both |
| `CellLoc` / `:cell` axis | `location="face"` / `n_face` dim | both |
| `cell_nodes(g,c)` `(SW,SE,NE,NW)` | `Mesh2_face_nodes`, `cf_role="face_node_connectivity"`, `start_index=1` | write |
| `node_coordinates(g,i)` `SVector{3}` (Cartesian) | `Mesh2_node_lon` / `Mesh2_node_lat` (via `_cartesian_to_latlon`) | write |
| `cell_centroid(g,c)` | `Mesh2_face_lon` / `Mesh2_face_lat` characteristic coordinates | write for `CellLoc` data |
| `edge_nodes(g,e)` midpoint | `Mesh2_edge_lon` / `Mesh2_edge_lat` characteristic coordinates | write for `EdgeLoc` data |
| `edge_nodes(g,e)` | `Mesh2_edge_nodes`, `cf_role="edge_node_connectivity"`, `start_index=1` | write when needed |
| `field.mesh` | data-var `mesh="Mesh2"` attribute | write |
| `num_cells` / `num_nodes` / `num_edges` | `n_face` / `n_node` / `n_edge` dim lengths | both |

**Verified correctness facts:** (1) `cell_nodes` order `(SW,SE,NE,NW)` is
already anticlockwise "as viewed from above" (outward) — UGRID-compatible with
no reordering. (2) All current grids are `IsUniform{4}`, so
`n_max_face_nodes = 4` and `_FillValue` is unused; the writer still handles
`IsMixed` (ragged) correctly for future grids.

### UGRID v1.0 conformance contract

The NetCDF writer must produce these UGRID structures:

- **Global attribute:** `Conventions = "CF-1.11 UGRID-1.0"` unless the caller
  supplies a stricter compatible value.
- **Mesh topology variable:** scalar integer variable `Mesh2` with
  `cf_role="mesh_topology"`, `topology_dimension=2`,
  `node_coordinates="Mesh2_node_lon Mesh2_node_lat"`, and
  `face_node_connectivity="Mesh2_face_nodes"`.
- **Element dimensions:** `n_node`, `n_face`, and, when edge topology or
  `EdgeLoc` data are written, `n_edge`. The writer may also emit
  `face_dimension="n_face"` and `edge_dimension="n_edge"` on `Mesh2` to make
  the element dimension mapping explicit even when using the default
  connectivity dimension order.
- **Node coordinate variables:** `Mesh2_node_lon(n_node)` and
  `Mesh2_node_lat(n_node)` with CF-compatible `standard_name` and `units`
  (`longitude` / `degrees_east`, `latitude` / `degrees_north`).
- **Face connectivity variable:** `Mesh2_face_nodes(n_face, n_max_face_nodes)`
  with integer type, `cf_role="face_node_connectivity"`, and `start_index=1`.
  For fixed quadrilateral grids `n_max_face_nodes = 4` and there are no missing
  indices. For future mixed grids, faces with fewer nodes use an explicit
  negative integer `_FillValue`; every face must still have at least three
  non-missing node indices.
- **Edge connectivity variable:** `Mesh2_edge_nodes(n_edge, Two)` with integer
  type, `cf_role="edge_node_connectivity"`, and `start_index=1` whenever edge
  data are written. It must not contain missing indices and should not carry
  `_FillValue`.
- **Characteristic coordinates:** when writing `CellLoc` data, emit
  `Mesh2_face_lon(n_face)` / `Mesh2_face_lat(n_face)` and set
  `Mesh2:face_coordinates = "Mesh2_face_lon Mesh2_face_lat"`. When writing
  `EdgeLoc` data, emit `Mesh2_edge_lon(n_edge)` / `Mesh2_edge_lat(n_edge)` and
  set `Mesh2:edge_coordinates = "Mesh2_edge_lon Mesh2_edge_lat"`.
- **Mesh data variable:** every field data variable has
  `mesh="Mesh2"`, `location=ugrid_location(Loc)`, and a CF `coordinates`
  attribute pointing to the characteristic coordinates for its UGRID location:
  node data → `Mesh2_node_lon Mesh2_node_lat`; face data →
  `Mesh2_face_lon Mesh2_face_lat`; edge data →
  `Mesh2_edge_lon Mesh2_edge_lat`.
- **Data variable dimensions:** the data variable must contain exactly one mesh
  element dimension, and that dimension must be the one corresponding to
  `location`: `n_node` for `NodeLoc`, `n_edge` for `EdgeLoc`, `n_face` for
  `CellLoc`. Non-location dimensions (`time`, `layer`, …) are preserved as
  ordinary CF/YAX dimensions.

### NetCDF round-trip mechanism

`to_ugrid` embeds a `grid_type` attribute (e.g. `"LatLonGrid"`) plus the
construction parameters on the topology variable, so `from_ugrid`
**losslessly reconstructs the exact `ManifoldMeshes` grid** for our own files.
For **external** UGRID files (no such attribute), v0 may do best-effort
structural detection if the case is trivial (e.g. a regular lon×lat structure →
`LatLonGrid`), but this is not part of the v0 contract. Otherwise it errors —
fully general read of arbitrary unstructured meshes is blocked on the `IsMesh`
concrete type (locate spec §11-4).

**v0 milestone:** NetCDF write path (`to_ugrid` / `save_ugrid`) solid; NetCDF
read path (`from_ugrid` / `load_ugrid`) handles our own single-field files via
`grid_type`. Zarr UGRID conformance is future work.

## 10. Error Handling

- **`DimensionMismatch`:** location-axis length ≠ `num_*`; mesh identity
  (`===`) violation in arithmetic; `Loc` mismatch in arithmetic.
- **`ArgumentError`:** data has no location axis; data is not 1-axis-compatible
  with the chosen `Loc`; plain `AbstractArray` was constructed without explicit
  dims; UGRID read has no `mesh_topology` variable; grid type undetectable on
  external read; v0 sees multiple UGRID data variables; UGRID read is missing
  required `cf_role`, `topology_dimension`, `node_coordinates`,
  `face_node_connectivity`, `mesh`, or `location` attributes.
- **`MethodError` / informative error:** `interpolate` on a non-`NodeLoc` field
  (CellLoc/EdgeLoc gather deferred).

## 11. Testing Plan

New package test suite. Follows the project's spot-check pattern plus
round-trip and conformance checks.

- **Construction + dim check:** one field per grid × per `Loc` type tag
  (`NodeLoc`, `EdgeLoc`, `CellLoc`) (4 × 3); mismatched-length →
  `DimensionMismatch`; missing/repeated location dim → `ArgumentError`; plain
  `Array` without dims → `ArgumentError`; multi-dim `(loc, time)` and
  `(time, loc)` construct.
- **Accessors:** `mesh`/`data`/`location` return the right objects; inherited
  `dims`/indexing/`sum` behave; non-location selectors work where axes are
  labeled.
- **Arithmetic:** `f .+ g` returns a `DiscreteField` with preserved axes; scalar
  `2 .* f`; chained `(f .+ g) .* 2` keeps the mesh; **different-mesh** →
  `DimensionMismatch`; **same structure but different mesh object** →
  `DimensionMismatch`; `withmesh(g, mesh(f))` then succeeds; **mixed `Loc`** →
  error; `sin.(f)` re-wraps; name/metadata/refdims preservation is fixed.
- **Gather:** round-trip at cell centroids (reuse the locate spec's parity
  guards); scalar result for `(node,)`; time-vector result for `(node, time)`
  and `(time, node)`; out-of-domain lat/lon errors match `locate_cell`;
  `CellLoc`/`EdgeLoc` gather errors informatively.
- **UGRID IO round-trip:** `to_ugrid(f) |> from_ugrid` yields equal mesh + data
  for each grid; `save_ugrid` → `load_ugrid` NetCDF file round-trip;
  dims/name/metadata/location/grid construction parameters survive;
  multi-variable UGRID raises the v0 `ArgumentError`.
- **UGRID conformance:** assert global `Conventions` contains `UGRID-1.0`;
  `Mesh2` is scalar with `cf_role="mesh_topology"`,
  `topology_dimension=2`, `node_coordinates`, and
  `face_node_connectivity`; `Mesh2_face_nodes` has
  `cf_role="face_node_connectivity"` and `start_index=1`; connectivity is
  anticlockwise and has no missing values for current quadrilateral grids;
  `CellLoc` writes `location="face"` and exactly one `n_face` data dimension;
  `NodeLoc` writes `location="node"` and exactly one `n_node` data dimension;
  `EdgeLoc` writes `edge_node_connectivity`, `edge_coordinates`,
  `location="edge"`, and exactly one `n_edge` data dimension; data variables
  include `mesh`, `location`, and `coordinates`.
- **Aqua** code-quality checks.

## 12. Decisions Log

| Decision | Chosen | Rejected alternative | Why |
|---|---|---|---|
| Field realization | **F1: `DiscreteField <: AbstractDimArray`** | F2 wrap `YAXArray` + custom BroadcastStyle; F3 hybrid | `YAXArray` is final; subtyping `AbstractDimArray` inherits DD broadcast/dims/selectors/reductions — least custom code, truest UXarray analogue |
| v0 scope | Container + arithmetic + gather + **NetCDF UGRID own-file round-trip** | Pure container only (literal §4b); full NetCDF+Zarr+external UGRID | User: build on YAXArrays from day 1; IO proves the UGRID borrowing and is the storage contract, but v0 stays narrow enough to finish solidly |
| Internal vocabulary | `cell` / `node` / `edge` everywhere | UGRID-native `face` in the field | Single vocab matches `ManifoldMeshes`; UGRID is an interchange format handled by explicit IO functions (reaffirms locate spec §3) |
| Location public API | type tag `DiscreteField(NodeLoc, ...)` + location trait methods | instance tag `DiscreteField(NodeLoc(), ...)`; symbol enums; heavyweight trait hierarchy | `DiscreteField{NodeLoc}` carries staggered location at the type level; constructors/IO dispatch cleanly; future locations add methods without changing `ManifoldMeshes` |
| Location representation | Flat named axis `:cell`/`:node`/`:edge` + non-location axes | Grid-native multi-axis (`:lon,:lat`) | Uniform across grids; composes with flat `interpolation_weights` node IDs; matches UXarray's `n_face` model |
| Data substrate | `YAXArrays` hard dep from v0 | weak dep / deferred | User: cubes are the substrate (NetCDF/Zarr, multi-dim, labeled axes) |
| Arithmetic inheritance | DD broadcast (free) | hand-rolled operators | Scales to any element-wise op; result stays wrapped; mesh-identity check in one place |
| Mesh-identity check | `===` (object identity) | structural `==` | O(1); catches the real bug; users share the mesh object |
| Mesh rebinding | explicit `withmesh(f, mesh)` | automatic equivalent-grid matching in broadcast | Keeps arithmetic strict and cheap; gives IO/reconstruction users an intentional escape hatch |
| `cell ↔ face` mapping location | Inside `to_ugrid`/`from_ugrid` only | field-layer rename | Keeps the field internally consistent; UGRID vocab only at the IO boundary |
| UGRID writer contract | UGRID v1.0-conformant NetCDF required attrs + selected CF/UGRID recommendations | merely UGRID-inspired structure | Lets v0 output be checked by UGRID conformance rules instead of relying on informal structural similarity |
| UGRID read generality (v0) | Our-files-first via `grid_type` attr | full arbitrary-mesh read now | Arbitrary read blocked on `IsMesh` concrete type (§11-4) |

## 13. Future Work (out of scope, recorded)

1. **`FieldSet`** — one mesh + named fields; the `UxDataset` analogue. `DiscreteField`
   is shaped to compose into it without redesign.
2. **`CellLoc` / `EdgeLoc` gather** — when `interpolation_weights` is extended
   to those locations.
3. **Conservative remap primitives** (`cell_intersection_area`, spherical
   polygon clipping) in the mesh layer + orchestration in `ManifoldRegrid.jl`.
4. **Analysis operators** (`gradient`, `integrate`) — on the field, using mesh
   primitives (`cell_volume`, `edge_outward_normal`). `cell_volume` already
   exists (≡ UXarray `face_areas`); `integrate(f) = sum(f .* cell_volume)` is a
   one-liner candidate.
5. **Generic external-UGRID read** of arbitrary unstructured meshes — blocked on
   the `IsMesh` concrete type (§11-4).
6. **Zarr UGRID conformance** — once the NetCDF UGRID model is stable.
7. **Out-of-core `mapCube` compute** operating directly on fields.
8. **`ManifoldRegrid.jl`** — the consumer: `remap(src_field, dest_mesh; method)`
   → new `DiscreteField` on the destination grid, using `locate_cell` +
   `interpolation_weights` + gather.
