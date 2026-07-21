# ManifoldFields.jl — Discrete Fields on Manifold Meshes — Design

- **Date:** 2026-07-21
- **Status:** Ready for review
- **Scope:** `ManifoldFields.jl` v0 — the field container that pairs a `ManifoldMeshes`
  mesh with a `YAXArrays` data cube, with element-wise arithmetic, point-gather
  interpolation, and bidirectional UGRID NetCDF IO.
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
  `DimensionalData` (dims, indexing, reductions, selectors)
- Element-wise arithmetic via inherited DD broadcast, with mesh-identity
  enforcement
- `interpolate(f, lat, lon)` point-gather (+ `SVector{3}` overload), `NodeLoc`
  only, multi-dimensional (reduces over the location axis, preserves trailing)
- Bidirectional UGRID IO: `to_ugrid` / `from_ugrid` / `save_ugrid` / `load_ugrid`
  (write path solid in v0; read path handles our own files first)

### Explicitly out of scope (v0)

- **`FieldSet`** (the `UxDataset` analogue — one mesh + named fields) — v0.1
- **`CellLoc` / `EdgeLoc` gather** — `interpolation_weights` is `NodeLoc`-only
  today; other locations carry more design choices and no consumer
- **Conservative remap, analysis operators** (`gradient`, `integrate`) — these
  live in `ManifoldRegrid.jl` / a future operators layer
- **Generic external-UGRID read of arbitrary unstructured meshes** — blocked on
  the `IsMesh` concrete type (locate spec §11-4)
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
DiscreteField(::Loc, mesh, data; name, metadata) where {Loc<:AbstractLocation}

# Accessors
mesh(f::DiscreteField)          -> AbstractManifoldMesh
data(f::DiscreteField)          -> AbstractArray     # the underlying cube/array
location(f::DiscreteField)      -> Type{<:AbstractLocation}

# Inherited from AbstractDimArray: dims, f[i], f[Near(lat)], sum/mean/extrema, f .+ g, ...

# Point gather (NodeLoc only in v0)
interpolate(f::DiscreteField, lat::Real, lon::Real)
interpolate(f::DiscreteField, p::StaticArrays.SVector{3})

# UGRID IO
to_ugrid(mesh::AbstractManifoldMesh)                -> UGridDataset
to_ugrid(field::DiscreteField)                      -> UGridDataset
from_ugrid(ds::UGridDataset; grid_type=nothing)     -> Tuple{mesh, fields...}
save_ugrid(x, path::AbstractString; format)         -> nothing
load_ugrid(path::AbstractString; grid_type=nothing) -> Tuple{mesh, fields...}
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
  attribute, promoted to a type. It selects the location axis and drives the
  dimension check.
- **Location axis** is a `DimensionalData` `Dimension` named `:node` / `:cell` /
  `:edge` — chosen by `Loc`, so the type and the dim agree (`CellLoc ↔ :cell`).
- **Trailing axes** are free (`:time`, `:layer`, …) and follow CF conventions;
  multi-dimensional fields `(loc, time, layer)` are first-class, not deferred.
- **`data`** is any `AbstractArray`; the blessed substrate is a `YAXArray`
  (`Cube`), so NetCDF/Zarr-backed out-of-core cubes plug in with no adapter.

### Loc ↔ dimension correspondence (defined in `ManifoldFields`)

This is the single, centralized dispatch table — `ManifoldMeshes` stays
pure-geometry-vocab; `ManifoldFields` owns the (trivial) `Loc` semantics:

```julia
location_dimname(::Type{NodeLoc}) = Dim{:node}
location_dimname(::Type{CellLoc}) = Dim{:cell}
location_dimname(::Type{EdgeLoc}) = Dim{:edge}

_expected_len(::NodeLoc, g) = num_nodes(g)
_expected_len(::CellLoc, g) = num_cells(g)
_expected_len(::EdgeLoc, g) = num_edges(g)
```

### Constructor — dimension check at construction time

```julia
function DiscreteField(::Loc, mesh::M, data; name=:field, metadata=nothing) where {Loc<:AbstractLocation,M}
    locdim = location_dimname(Loc)
    expected = _expected_len(Loc(), mesh)
    # require a location axis of matching length; construct dims accordingly
    …
    DiscreteField{Loc,…}(data, dims, …, mesh)
end
```

The whole point of carrying `Loc`: a field/grid mismatch fails at
**construction**, not at first use.

> **Implementation note.** The exact `AbstractDimArray` interface contract
> (required accessors: `DimensionalData.data/dims/refdims/name/metadata`, plus
> the broadcast `similar`/style hooks) must be validated against the installed
> `DimensionalData` version during implementation. The field set above is the
> expected shape; details are a plan-phase item.

## 6. Accessors & Inherited Behaviour

- `mesh(f)`, `data(f)`, `location(f)` (returns the `Loc` **type**).
- **Inherited from `AbstractDimArray`** (this is the payoff of F1): `dims(f)`,
  indexing `f[i]` / `f[Di=…]`, reductions `sum` / `mean` / `extrema`, DD
  selectors `f[Near(lat)]`, and element-wise broadcast.
- **No `setindex!`.** A `DiscreteField` is not written through (disk-cube
  safety). Mutation, when needed, is via `data(f)` at the caller's risk.

## 7. Arithmetic (inherited from DimensionalData)

Because `DiscreteField <: AbstractDimArray`, DD's broadcast machinery preserves
type + dims through element-wise ops: `f .+ g`, `2 .* f`, `sin.(f)`, `f .* g`
all return a `DiscreteField` with the same axes.

- **Single custom hook** — override the result-build (`similar` / broadcast
  `copy`) so that, across all `DiscreteField` args: (a) meshes are identical
  (`===`), else `DimensionMismatch`; (b) `Loc` matches, else
  `DimensionMismatch`; (c) the result carries the shared mesh + `Loc`. Scalars
  and plain numbers combine naturally (no mesh check against them).
- **Mesh-identity rule:** `f.mesh === g.mesh` (object identity), not structural
  `==`. O(1), and it catches the real bug (fields on different grids). Users
  who reconstruct "the same" grid share the object.
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
  and **preserves trailing axes**. A `NodeLoc` field on `(node,)` returns a
  scalar; a field on `(node, time)` returns a length-`T` vector (the
  interpolated time series at the query point); etc.
- **`NodeLoc` only in v0** — `interpolation_weights` is `NodeLoc`-only today.
  Calling `interpolate` on a `CellLoc` / `EdgeLoc` field errors with an
  informative message (deferred, see §12).
- This is the direct analogue of UXarray's `.interp(lon_lat=(…), method=)`.

## 9. UGRID IO (v0)

We borrow UGRID's **structure** and confine the `cell ↔ face` vocabulary
translation to these functions. `ManifoldMeshes` stays pure (no NetCDF/UGRID
dependency); the IO lives in `ManifoldFields`, which calls mesh primitives and
YAXArrays for the bytes.

```julia
# in-memory conversion (DiscreteField ⇄ UGRID dataset)
to_ugrid(mesh::AbstractManifoldMesh)   -> UGridDataset   # topology only
to_ugrid(field::DiscreteField)         -> UGridDataset   # topology + one data var
from_ugrid(ds::UGridDataset; grid_type=nothing) -> Tuple{mesh, fields...}

# file IO (bytes delegated to YAXArrays; ManifoldFields adds UGRID scaffolding)
save_ugrid(x, path; format ∈ {":netcdf", ":zarr"})
load_ugrid(path; grid_type=nothing) = from_ugrid(YAXArrays.read(path))
```

### `UGridDataset` and return shapes

`UGridDataset` is the in-memory UGRID representation: a `YAXArrays`-backed
dataset holding the scalar `mesh_topology` variable, the coordinate and
connectivity cubes (`node_lon` / `node_lat`, `face_node_connectivity`), and the
data cube(s) carrying `mesh` / `location` attributes. Its exact Julia type is a
plan-phase detail (a `YAXArrays.Dataset` or a thin wrapper).

`from_ugrid` returns `(mesh, fields)`. In v0 (no `FieldSet`) it handles
single-data-variable datasets — returning `(mesh, field)`; multi-variable
datasets require `FieldSet` (v0.1).

### Vocabulary translation (lives *only* in `to_ugrid` / `from_ugrid`)

| Internal (`ManifoldMeshes` vocab) | UGRID NetCDF | Direction |
|---|---|---|
| `CellLoc` / `:cell` axis | `location="face"` / `n_face` dim | both |
| `cell_nodes(g,c)` `(SW,SE,NE,NW)` | `face_node_connectivity`, `start_index=1`, `_FillValue` for ragged | write |
| `node_coordinates(g,i)` `SVector{3}` (Cartesian) | `node_lon` / `node_lat` (via `_cartesian_to_latlon`) | write |
| `edge_nodes(g,e)` | `edge_node_connectivity` | write (optional) |
| `field.mesh` | data-var `mesh="Mesh2"` attribute | write |
| `num_cells` / `num_nodes` / `num_edges` | `n_face` / `n_node` / `n_edge` dim lengths | both |

**Verified correctness facts:** (1) `cell_nodes` order `(SW,SE,NE,NW)` is
already anticlockwise "as viewed from above" (outward) — UGRID-compatible with
no reordering. (2) All current grids are `IsUniform{4}`, so
`n_max_face_nodes = 4` and `_FillValue` is unused; the writer still handles
`IsMixed` (ragged) correctly for future grids.

### Round-trip mechanism

`to_ugrid` embeds a `grid_type` attribute (e.g. `"LatLonGrid"`) plus the
construction parameters on the topology variable, so `from_ugrid`
**losslessly reconstructs the exact `ManifoldMeshes` grid** for our own files.
For **external** UGRID files (no such attribute), `from_ugrid` does best-effort
structural detection (e.g. a regular lon×lat structure → `LatLonGrid`);
otherwise it errors — fully general read of arbitrary unstructured meshes is
blocked on the `IsMesh` concrete type (locate spec §11-4).

**v0 milestone:** write path (`to_ugrid` / `save_ugrid`) solid; read path
(`from_ugrid` / `load_ugrid`) handles our own files (via `grid_type`) first,
with best-effort external detection.

## 10. Error Handling

- **`DimensionMismatch`:** location-axis length ≠ `num_*`; mesh identity
  (`===`) violation in arithmetic; `Loc` mismatch in arithmetic.
- **`ArgumentError`:** data has no location axis; data is not 1-axis-compatible
  with the chosen `Loc`; UGRID read has no `mesh_topology` variable; grid type
  undetectable on external read.
- **`MethodError` / informative error:** `interpolate` on a non-`NodeLoc` field
  (CellLoc/EdgeLoc gather deferred).

## 11. Testing Plan

New package test suite. Follows the project's spot-check pattern plus
round-trip and conformance checks.

- **Construction + dim check:** one field per grid × per `Loc` (4 × 3);
  mismatched-length → `DimensionMismatch`; multi-dim `(loc, time)` constructs.
- **Accessors:** `mesh`/`data`/`location` return the right objects; inherited
  `dims`/indexing/`sum` behave.
- **Arithmetic:** `f .+ g` returns a `DiscreteField` with preserved axes; scalar
  `2 .* f`; chained `(f .+ g) .* 2` keeps the mesh; **different-mesh** →
  `DimensionMismatch`; **mixed `Loc`** → error; `sin.(f)` re-wraps.
- **Gather:** round-trip at cell centroids (reuse the locate spec's parity
  guards); scalar result for `(node,)`; time-vector result for `(node, time)`.
- **UGRID IO round-trip:** `to_ugrid(f) |> from_ugrid` yields equal mesh + data
  for each grid; conformance — `cf_role="mesh_topology"`,
  `topology_dimension=2`, `start_index=1`, `location="face"` for `CellLoc`,
  anticlockwise connectivity; `save_ugrid` → `load_ugrid` file round-trip.
- **Aqua** code-quality checks.

## 12. Decisions Log

| Decision | Chosen | Rejected alternative | Why |
|---|---|---|---|
| Field realization | **F1: `DiscreteField <: AbstractDimArray`** | F2 wrap `YAXArray` + custom BroadcastStyle; F3 hybrid | `YAXArray` is final; subtyping `AbstractDimArray` inherits DD broadcast/dims/selectors/reductions — least custom code, truest UXarray analogue |
| v0 scope | Container + arithmetic + gather + **UGRID IO** | Pure container only (literal §4b) | User: build on YAXArrays from day 1; IO proves the UGRID borrowing and is the storage contract |
| Internal vocabulary | `cell` / `node` / `edge` everywhere | UGRID-native `face` in the field | Single vocab matches `ManifoldMeshes`; UGRID is an interchange format handled by explicit IO functions (reaffirms locate spec §3) |
| Location representation | Flat named axis `:cell`/`:node`/`:edge` + trailing | Grid-native multi-axis (`:lon,:lat`) | Uniform across grids; composes with flat `interpolation_weights` node IDs; matches UXarray's `n_face` model |
| Data substrate | `YAXArrays` hard dep from v0 | weak dep / deferred | User: cubes are the substrate (NetCDF/Zarr, multi-dim, labeled axes) |
| Arithmetic inheritance | DD broadcast (free) | hand-rolled operators | Scales to any element-wise op; result stays wrapped; mesh-identity check in one place |
| Mesh-identity check | `===` (object identity) | structural `==` | O(1); catches the real bug; users share the mesh object |
| `cell ↔ face` mapping location | Inside `to_ugrid`/`from_ugrid` only | field-layer rename | Keeps the field internally consistent; UGRID vocab only at the IO boundary |
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
6. **Out-of-core `mapCube` compute** operating directly on fields.
7. **`ManifoldRegrid.jl`** — the consumer: `remap(src_field, dest_mesh; method)`
   → new `DiscreteField` on the destination grid, using `locate_cell` +
   `interpolation_weights` + gather.
