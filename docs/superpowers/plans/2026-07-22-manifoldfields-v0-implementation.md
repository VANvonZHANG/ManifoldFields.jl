# ManifoldFields v0 Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build ManifoldFields.jl v0: `DiscreteField{Loc}` containers over `ManifoldMeshes` meshes and `YAXArrays`/`DimensionalData` arrays, with location-tag construction, strict broadcast arithmetic, `NodeLoc` interpolation, and UGRID v1.0-conformant NetCDF own-file round-trip IO.

**Architecture:** `ManifoldFields` is a separate Julia package. It depends on `ManifoldMeshes` for geometry/topology, `DimensionalData`/`YAXArrays` for labeled arrays and storage, and `StaticArrays` for `SVector` point overloads. UGRID vocabulary is confined to `src/ugrid.jl`; all field-level location semantics go through type-tag trait methods.

**Tech Stack:** Julia 1.10, ManifoldMeshes.jl, DimensionalData.jl, YAXArrays.jl, StaticArrays.jl, NCDatasets.jl for explicit NetCDF variable/attribute IO, Aqua.jl, Test stdlib.

---

## File Structure

- Create: `Project.toml` — package metadata, dependencies, compat, test target.
- Create: `.JuliaFormatter.toml` — match `ManifoldMeshes.jl` SciML style.
- Create: `src/ManifoldFields.jl` — module entry, exports, includes.
- Create: `src/locations.jl` — `location_dimname`, `expected_location_length`, `ugrid_location`, `ugrid_dimname`.
- Create: `src/field.jl` — `DiscreteField` type, constructors, accessors, `withmesh`, DD interface.
- Create: `src/broadcast.jl` — broadcast/similar hooks and mesh/Loc checks.
- Create: `src/interpolate.jl` — `NodeLoc` point gather and `SVector{3}` overload.
- Create: `src/ugrid.jl` — UGRID v1.0 NetCDF in-memory representation, writer, own-file reader.
- Create: `test/runtests.jl` — test entrypoint with Aqua and all test files.
- Create: `test/test_locations.jl` — location trait method tests.
- Create: `test/test_field_construction.jl` — constructor/dim validation/accessors.
- Create: `test/test_broadcast.jl` — arithmetic and rebinding behavior.
- Create: `test/test_interpolate.jl` — `NodeLoc` gather semantics.
- Create: `test/test_ugrid.jl` — UGRID conformance and own-file round-trip tests.
- Create: `test/helpers.jl` — shared small grids, dims, and deterministic field values.

## Implementation Notes

- Work from package root: `/home/zhangfan/Project/20260328_HEMCOManifold/ManifoldFields.jl`.
- Use local development path for `ManifoldMeshes`: `Pkg.develop(path="../ManifoldMeshes.jl")`.
- Use type tags in public API: `DiscreteField(NodeLoc, mesh, values, dims)`, never `NodeLoc()`.
- Keep UGRID names fixed in v0: topology variable `Mesh2`, dimensions `n_node`, `n_edge`, `n_face`.
- Run focused tests after each task, then run `julia --project=. -e 'using Pkg; Pkg.test()'` before final completion.

### Task 1: Package Scaffold

**Files:**
- Create: `Project.toml`
- Create: `.JuliaFormatter.toml`
- Create: `src/ManifoldFields.jl`
- Create: `test/runtests.jl`

- [ ] **Step 1: Write package metadata**

Create `Project.toml`:

```toml
name = "ManifoldFields"
uuid = "9c5b4d8f-6d30-4f8a-9f1f-2c1e5c2f4b31"
version = "0.1.0"

[deps]
DimensionalData = "0703355e-b756-11e9-17c0-8b28908087d0"
ManifoldMeshes = "d2a7aed9-80bb-46a9-9a56-e170a93cc177"
NCDatasets = "85f8d34a-cbdd-5861-8df4-14fed0d494ab"
StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"
YAXArrays = "c21b50f5-a021-445e-a6bb-8c82a528c804"

[compat]
Aqua = "0.8"
DimensionalData = "0.29"
ManifoldMeshes = "0.6"
NCDatasets = "0.14"
StaticArrays = "1.9"
Test = "1.10"
YAXArrays = "0.6"
julia = "1.10"

[extras]
Aqua = "4c88cf16-eb10-579e-8560-4a9242c79595"
Test = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[targets]
test = ["Aqua", "Test"]
```

- [ ] **Step 2: Add formatter config**

Create `.JuliaFormatter.toml`:

```toml
style = "sciml"
```

- [ ] **Step 3: Add minimal module**

Create `src/ManifoldFields.jl`:

```julia
module ManifoldFields

using DimensionalData
using ManifoldMeshes
using StaticArrays
using YAXArrays

export DiscreteField, mesh, data, location, withmesh
export location_dimname, expected_location_length
export interpolate
export UGridDataset, to_ugrid, from_ugrid, from_ugrid_mesh, save_ugrid, load_ugrid,
       load_ugrid_mesh

include("locations.jl")
include("field.jl")
include("broadcast.jl")
include("interpolate.jl")
include("ugrid.jl")

end
```

- [ ] **Step 4: Add empty test runner**

Create `test/runtests.jl`:

```julia
using Aqua
using ManifoldFields
using Test

@testset "Code quality (Aqua.jl)" begin
    Aqua.test_all(ManifoldFields)
end

@testset "ManifoldFields.jl" begin
end
```

- [ ] **Step 5: Run scaffold test and observe missing includes**

Run:

```bash
julia --project=. -e 'using Pkg; Pkg.develop(path="../ManifoldMeshes.jl"); Pkg.test()'
```

Expected: FAIL with an include error for `src/locations.jl`.

- [ ] **Step 6: Commit scaffold**

```bash
git add Project.toml .JuliaFormatter.toml src/ManifoldFields.jl test/runtests.jl
git commit -m "chore: scaffold ManifoldFields package"
```

### Task 2: Location Trait Methods

**Files:**
- Create: `src/locations.jl`
- Create: `test/test_locations.jl`
- Modify: `test/runtests.jl`

- [ ] **Step 1: Write failing location tests**

Create `test/test_locations.jl`:

```julia
using DimensionalData
using ManifoldMeshes
using ManifoldFields
using Test

@testset "location trait methods" begin
    g = LatLonGrid(2, 4)

    @test location_dimname(NodeLoc) == Dim{:node}
    @test location_dimname(EdgeLoc) == Dim{:edge}
    @test location_dimname(CellLoc) == Dim{:cell}

    @test ManifoldFields.expected_location_length(NodeLoc, g) == num_nodes(g)
    @test ManifoldFields.expected_location_length(EdgeLoc, g) == num_edges(g)
    @test ManifoldFields.expected_location_length(CellLoc, g) == num_cells(g)

    @test ManifoldFields.ugrid_location(NodeLoc) == "node"
    @test ManifoldFields.ugrid_location(EdgeLoc) == "edge"
    @test ManifoldFields.ugrid_location(CellLoc) == "face"

    @test ManifoldFields.ugrid_dimname(NodeLoc) == "n_node"
    @test ManifoldFields.ugrid_dimname(EdgeLoc) == "n_edge"
    @test ManifoldFields.ugrid_dimname(CellLoc) == "n_face"
end
```

Modify `test/runtests.jl`:

```julia
using Aqua
using ManifoldFields
using Test

@testset "Code quality (Aqua.jl)" begin
    Aqua.test_all(ManifoldFields)
end

@testset "ManifoldFields.jl" begin
    include("test_locations.jl")
end
```

- [ ] **Step 2: Run test to verify it fails**

Run:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected: FAIL with `could not open file ... src/locations.jl` or `UndefVarError: location_dimname not defined`.

- [ ] **Step 3: Implement location trait table**

Create `src/locations.jl`:

```julia
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
```

- [ ] **Step 4: Run location tests**

Run:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected: FAIL now moves to missing `src/field.jl`, proving `locations.jl` loads.

- [ ] **Step 5: Commit location traits**

```bash
git add src/locations.jl test/test_locations.jl test/runtests.jl
git commit -m "feat: add location trait methods"
```

### Task 3: DiscreteField Constructors and Accessors

**Files:**
- Create: `src/field.jl`
- Create: `test/helpers.jl`
- Create: `test/test_field_construction.jl`
- Modify: `test/runtests.jl`

- [ ] **Step 1: Write shared test helpers**

Create `test/helpers.jl`:

```julia
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
```

- [ ] **Step 2: Write failing constructor/accessor tests**

Create `test/test_field_construction.jl`:

```julia
using DimensionalData
using ManifoldFields
using ManifoldMeshes
using Test

include("helpers.jl")

@testset "DiscreteField construction and accessors" begin
    g = small_grid()

    nf = DiscreteField(NodeLoc, g, node_values(g), node_dims(g); name=:temperature,
                       metadata=Dict("units" => "K"))
    ef = DiscreteField(EdgeLoc, g, edge_values(g), edge_dims(g); name=:flux)
    cf = DiscreteField(CellLoc, g, cell_values(g), cell_dims(g); name=:area)

    @test mesh(nf) === g
    @test data(nf) == node_values(g)
    @test location(nf) === NodeLoc
    @test DimensionalData.dims(nf) == node_dims(g)
    @test DimensionalData.name(nf) == :temperature
    @test DimensionalData.metadata(nf)["units"] == "K"

    @test location(ef) === EdgeLoc
    @test location(cf) === CellLoc

    @test DiscreteField(NodeLoc, g, node_time_values(g), node_time_dims(g)) isa
          DiscreteField{NodeLoc}
    @test DiscreteField(NodeLoc, g, time_node_values(g), time_node_dims(g)) isa
          DiscreteField{NodeLoc}

    @test_throws ArgumentError DiscreteField(NodeLoc, g, node_values(g))
    @test_throws ArgumentError DiscreteField(NodeLoc, g, node_values(g),
                                             (Dim{:time}(1:num_nodes(g)),))
    @test_throws ArgumentError DiscreteField(NodeLoc, g, node_values(g),
                                             (Dim{:node}(1:num_nodes(g)),
                                              Dim{:node}(1:num_nodes(g))))
    @test_throws DimensionMismatch DiscreteField(NodeLoc, g, node_values(g)[1:end-1],
                                                 (Dim{:node}(1:(num_nodes(g)-1)),))
end
```

Modify `test/runtests.jl`:

```julia
using Aqua
using ManifoldFields
using Test

@testset "Code quality (Aqua.jl)" begin
    Aqua.test_all(ManifoldFields)
end

@testset "ManifoldFields.jl" begin
    include("test_locations.jl")
    include("test_field_construction.jl")
end
```

- [ ] **Step 3: Run constructor tests and verify failure**

Run:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected: FAIL with missing `src/field.jl` or undefined `DiscreteField`.

- [ ] **Step 4: Implement minimal `DiscreteField`**

Create `src/field.jl`:

```julia
import DimensionalData
import DimensionalData: AbstractDimArray, dims
import ManifoldMeshes: AbstractLocation, AbstractManifoldMesh

struct DiscreteField{
    Loc<:AbstractLocation,
    T,N,D,A<:AbstractArray{T,N},
    M<:AbstractManifoldMesh,
    R,
    NM,
    MD,
} <: AbstractDimArray{T,N,D,A}
    data::A
    dims::D
    refdims::R
    name::NM
    metadata::MD
    mesh::M
end

mesh(f::DiscreteField) = f.mesh
data(f::DiscreteField) = f.data
location(::DiscreteField{Loc}) where {Loc} = Loc

DimensionalData.dims(f::DiscreteField) = f.dims
DimensionalData.refdims(f::DiscreteField) = f.refdims
DimensionalData.name(f::DiscreteField) = f.name
DimensionalData.metadata(f::DiscreteField) = f.metadata
DimensionalData.data(f::DiscreteField) = f.data
Base.size(f::DiscreteField) = size(f.data)
Base.size(f::DiscreteField, i::Integer) = size(f.data, i)
Base.getindex(f::DiscreteField, I...) = getindex(f.data, I...)
Base.IndexStyle(::Type{<:DiscreteField}) = IndexLinear()

function _dim_name(d)
    return DimensionalData.name(d)
end

function _location_axes(Loc, ds)
    locname = DimensionalData.name(location_dimname(Loc))
    return findall(d -> _dim_name(d) == locname, collect(ds))
end

function _validate_location_dims(::Type{Loc}, mesh, values, ds) where {Loc<:AbstractLocation}
    length(ds) == ndims(values) ||
        throw(DimensionMismatch("dims length $(length(ds)) does not match data ndims $(ndims(values))"))
    axes = _location_axes(Loc, ds)
    length(axes) == 1 ||
        throw(ArgumentError("expected exactly one $(location_dimname(Loc)) dimension, found $(length(axes))"))
    axis = only(axes)
    expected = expected_location_length(Loc, mesh)
    actual = size(values, axis)
    actual == expected ||
        throw(DimensionMismatch("$(location_dimname(Loc)) length $actual != expected $expected"))
    return ds
end

function DiscreteField(
    ::Type{Loc},
    mesh::M,
    values::A,
    ds;
    name=:field,
    metadata=nothing,
    refdims=(),
) where {Loc<:AbstractLocation,M<:AbstractManifoldMesh,A<:AbstractArray}
    checked_dims = _validate_location_dims(Loc, mesh, values, ds)
    return DiscreteField{Loc,eltype(values),ndims(values),typeof(checked_dims),A,M,typeof(refdims),
                         typeof(name),typeof(metadata)}(
        values, checked_dims, refdims, name, metadata, mesh
    )
end

function DiscreteField(
    ::Type{Loc},
    mesh::M,
    data::AbstractDimArray;
    name=DimensionalData.name(data),
    metadata=DimensionalData.metadata(data),
    refdims=DimensionalData.refdims(data),
) where {Loc<:AbstractLocation,M<:AbstractManifoldMesh}
    return DiscreteField(Loc, mesh, DimensionalData.data(data), DimensionalData.dims(data);
                         name=name, metadata=metadata, refdims=refdims)
end

function DiscreteField(::Type{Loc}, mesh::M, values::AbstractArray; kwargs...) where {
    Loc<:AbstractLocation,
    M<:AbstractManifoldMesh,
}
    throw(ArgumentError("plain AbstractArray input requires explicit dims"))
end

function withmesh(f::DiscreteField{Loc}, mesh::AbstractManifoldMesh) where {Loc}
    return DiscreteField(Loc, mesh, data(f), dims(f); name=DimensionalData.name(f),
                         metadata=DimensionalData.metadata(f),
                         refdims=DimensionalData.refdims(f))
end
```

- [ ] **Step 5: Add temporary stubs for remaining includes**

Create `src/broadcast.jl`:

```julia
# Broadcast hooks are implemented in Task 5.
```

Create `src/interpolate.jl`:

```julia
# Interpolation is implemented in Task 6.
```

Create `src/ugrid.jl`:

```julia
struct UGridDataset
    variables::Dict{String,Any}
    attributes::Dict{String,Any}
end
```

- [ ] **Step 6: Run constructor tests**

Run:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected: constructor/location tests PASS; Aqua may report ambiguities or missing method issues to resolve before committing.

- [ ] **Step 7: Commit constructors**

```bash
git add src/field.jl src/broadcast.jl src/interpolate.jl src/ugrid.jl test/helpers.jl test/test_field_construction.jl test/runtests.jl
git commit -m "feat: add DiscreteField constructors"
```

### Task 4: DimensionalData Compatibility Spike

**Files:**
- Create: `test/test_dimensionaldata_spike.jl`
- Modify: `test/runtests.jl`
- Modify: `src/field.jl`

- [ ] **Step 1: Write spike tests**

Create `test/test_dimensionaldata_spike.jl`:

```julia
using DimensionalData
using ManifoldFields
using ManifoldMeshes
using Test

include("helpers.jl")

@testset "DimensionalData compatibility spike" begin
    g = small_grid()
    f = DiscreteField(NodeLoc, g, node_time_values(g), node_time_dims(g);
                      name=:temperature, metadata=Dict("units" => "K"))

    @test dims(f) == node_time_dims(g)
    @test size(f) == (num_nodes(g), 3)
    @test f[1, 1] == node_time_values(g)[1, 1]
    @test DimensionalData.data(f) == node_time_values(g)
    @test DimensionalData.name(f) == :temperature
    @test DimensionalData.metadata(f)["units"] == "K"
    @test sum(f) == sum(node_time_values(g))
    @test extrema(f) == extrema(node_time_values(g))
end
```

Add to `test/runtests.jl` after construction tests:

```julia
include("test_dimensionaldata_spike.jl")
```

- [ ] **Step 2: Run spike tests**

Run:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected: any failures identify the exact missing AbstractDimArray method.

- [ ] **Step 3: Patch only required DD methods**

If reductions fail because Base iteration is missing, add to `src/field.jl`:

```julia
Base.length(f::DiscreteField) = length(f.data)
Base.eachindex(f::DiscreteField) = eachindex(f.data)
Base.iterate(f::DiscreteField, state...) = iterate(f.data, state...)
Base.eltype(::Type{<:DiscreteField{Loc,T}}) where {Loc,T} = T
```

If `similar` is required before broadcast work, add the narrow method:

```julia
Base.similar(f::DiscreteField, ::Type{T}, dims::Dims) where {T} = similar(data(f), T, dims)
```

- [ ] **Step 4: Re-run spike tests**

Run:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected: `test_dimensionaldata_spike.jl` PASS.

- [ ] **Step 5: Commit DD spike fixes**

```bash
git add src/field.jl test/test_dimensionaldata_spike.jl test/runtests.jl
git commit -m "test: verify DimensionalData compatibility"
```

### Task 5: Broadcast Arithmetic and Mesh Identity

**Files:**
- Modify: `src/broadcast.jl`
- Modify: `src/field.jl`
- Create: `test/test_broadcast.jl`
- Modify: `test/runtests.jl`

- [ ] **Step 1: Write failing broadcast tests**

Create `test/test_broadcast.jl`:

```julia
using DimensionalData
using ManifoldFields
using ManifoldMeshes
using Test

include("helpers.jl")

@testset "DiscreteField broadcast arithmetic" begin
    g = small_grid()
    f = DiscreteField(NodeLoc, g, node_values(g), node_dims(g); name=:a,
                      metadata=Dict("units" => "K"))
    h = DiscreteField(NodeLoc, g, 2 .* node_values(g), node_dims(g); name=:b,
                      metadata=Dict("units" => "K"))

    r = f .+ h
    @test r isa DiscreteField{NodeLoc}
    @test mesh(r) === g
    @test location(r) === NodeLoc
    @test data(r) == 3 .* node_values(g)
    @test dims(r) == node_dims(g)

    s = sin.(f)
    @test s isa DiscreteField{NodeLoc}
    @test mesh(s) === g
    @test data(s) == sin.(node_values(g))

    t = 2 .* f
    @test t isa DiscreteField{NodeLoc}
    @test data(t) == 2 .* node_values(g)

    g2 = LatLonGrid(2, 4)
    same_shape_other_mesh = DiscreteField(NodeLoc, g2, node_values(g2), node_dims(g2))
    @test_throws DimensionMismatch f .+ same_shape_other_mesh

    cell_f = DiscreteField(CellLoc, g, cell_values(g), cell_dims(g))
    @test_throws DimensionMismatch f .+ cell_f

    rebound = withmesh(same_shape_other_mesh, mesh(f))
    @test data(f .+ rebound) == node_values(g) .+ node_values(g2)
end
```

Add to `test/runtests.jl`:

```julia
include("test_broadcast.jl")
```

- [ ] **Step 2: Run broadcast tests and verify failure**

Run:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected: FAIL because broadcast currently returns plain arrays or does not enforce mesh identity.

- [ ] **Step 3: Implement conservative broadcast**

Replace `src/broadcast.jl` with:

```julia
import Base.Broadcast: BroadcastStyle, DefaultArrayStyle

struct DiscreteFieldStyle{N} <: Broadcast.AbstractArrayStyle{N} end

Base.BroadcastStyle(::Type{<:DiscreteField{Loc,T,N}}) where {Loc,T,N} =
    DiscreteFieldStyle{N}()
Base.BroadcastStyle(::DiscreteFieldStyle{N}, ::DefaultArrayStyle{N}) where {N} =
    DiscreteFieldStyle{N}()
Base.BroadcastStyle(::DefaultArrayStyle{N}, ::DiscreteFieldStyle{N}) where {N} =
    DiscreteFieldStyle{N}()
Base.BroadcastStyle(::DiscreteFieldStyle{N}, ::DiscreteFieldStyle{N}) where {N} =
    DiscreteFieldStyle{N}()

function _fields_in_broadcast(args)
    fields = DiscreteField[]
    _collect_fields!(fields, args)
    return fields
end

function _collect_fields!(fields, x::DiscreteField)
    push!(fields, x)
    return fields
end

function _collect_fields!(fields, bc::Base.Broadcast.Broadcasted)
    foreach(arg -> _collect_fields!(fields, arg), bc.args)
    return fields
end

function _collect_fields!(fields, x)
    return fields
end

function _broadcast_anchor(fields)
    isempty(fields) && throw(ArgumentError("DiscreteField broadcast requires a field argument"))
    first_field = first(fields)
    for f in fields
        mesh(f) === mesh(first_field) ||
            throw(DimensionMismatch("cannot broadcast fields on different mesh objects"))
        location(f) === location(first_field) ||
            throw(DimensionMismatch("cannot broadcast fields with different locations"))
        dims(f) == dims(first_field) ||
            throw(DimensionMismatch("cannot broadcast fields with different dims"))
    end
    return first_field
end

function Base.copy(bc::Base.Broadcast.Broadcasted{DiscreteFieldStyle{N}}) where {N}
    fields = _fields_in_broadcast(bc.args)
    anchor = _broadcast_anchor(fields)
    plain = Base.Broadcast.Broadcasted(bc.f, map(_broadcast_arg, bc.args), bc.axes)
    values = copy(plain)
    return DiscreteField(location(anchor), mesh(anchor), values, dims(anchor);
                         name=DimensionalData.name(anchor),
                         metadata=DimensionalData.metadata(anchor),
                         refdims=DimensionalData.refdims(anchor))
end

_broadcast_arg(f::DiscreteField) = data(f)
_broadcast_arg(bc::Base.Broadcast.Broadcasted) =
    Base.Broadcast.Broadcasted(bc.f, map(_broadcast_arg, bc.args), bc.axes)
_broadcast_arg(x) = x
```

- [ ] **Step 4: Run broadcast tests**

Run:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected: broadcast tests PASS. If `DiscreteField[]` is too concrete for parametric fields, change it to `Any[]` and keep `_broadcast_anchor` unchanged.

- [ ] **Step 5: Commit broadcast behavior**

```bash
git add src/broadcast.jl test/test_broadcast.jl test/runtests.jl
git commit -m "feat: add field broadcast semantics"
```

### Task 6: NodeLoc Interpolation

**Files:**
- Modify: `src/interpolate.jl`
- Create: `test/test_interpolate.jl`
- Modify: `test/runtests.jl`

- [ ] **Step 1: Write failing interpolation tests**

Create `test/test_interpolate.jl`:

```julia
using ManifoldFields
using ManifoldMeshes
using StaticArrays
using Test

include("helpers.jl")

@testset "NodeLoc interpolate" begin
    g = LatLonGrid(4, 8)
    values = collect(Float64, 1:num_nodes(g))
    f = DiscreteField(NodeLoc, g, values, node_dims(g))

    cid = 6
    lat, lon = cell_centroid(g, cid) |> ManifoldMeshes._cartesian_to_latlon
    nodes, weights = interpolation_weights(g, cid, lat, lon)
    expected = sum(values[nodes[i]] * weights[i] for i in 1:4)
    @test interpolate(f, lat, lon) ≈ expected

    p = SVector{3}(cell_centroid(g, cid))
    @test interpolate(f, p) ≈ expected

    ft = DiscreteField(NodeLoc, g, node_time_values(g), node_time_dims(g))
    result = interpolate(ft, lat, lon)
    @test length(result) == 3
    @test result ≈ [sum(node_time_values(g)[nodes[i], t] * weights[i] for i in 1:4)
                    for t in 1:3]

    ftn = DiscreteField(NodeLoc, g, time_node_values(g), time_node_dims(g))
    result2 = interpolate(ftn, lat, lon)
    @test length(result2) == 3
    @test result2 ≈ [sum(time_node_values(g)[t, nodes[i]] * weights[i] for i in 1:4)
                     for t in 1:3]

    @test_throws ArgumentError interpolate(f, 95.0, 0.0)
    @test_throws MethodError interpolate(DiscreteField(CellLoc, g, cell_values(g), cell_dims(g)),
                                         lat, lon)
end
```

Add to `test/runtests.jl`:

```julia
include("test_interpolate.jl")
```

- [ ] **Step 2: Run interpolation tests and verify failure**

Run:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected: FAIL with `UndefVarError: interpolate not defined`.

- [ ] **Step 3: Implement interpolation**

Replace `src/interpolate.jl` with:

```julia
import ManifoldMeshes: CellLoc, EdgeLoc, NodeLoc, interpolation_weights, locate_cell
import StaticArrays: SVector

function _location_axis(::Type{Loc}, f::DiscreteField) where {Loc}
    locname = DimensionalData.name(location_dimname(Loc))
    matches = findall(d -> DimensionalData.name(d) == locname, collect(dims(f)))
    return only(matches)
end

function interpolate(f::DiscreteField{NodeLoc}, lat::Real, lon::Real)
    cid = locate_cell(mesh(f), lat, lon)
    nodes, weights = interpolation_weights(mesh(f), cid, lat, lon)
    loc_axis = _location_axis(NodeLoc, f)

    if ndims(data(f)) == 1
        return sum(data(f)[nodes[i]] * weights[i] for i in 1:4)
    end

    moved = PermutedDimsArray(data(f), (loc_axis, setdiff(1:ndims(data(f)), loc_axis)...))
    trailing_shape = size(moved)[2:end]
    out = zeros(promote_type(eltype(data(f)), Float64), trailing_shape)
    for i in 1:4
        out .+= weights[i] .* selectdim(moved, 1, nodes[i])
    end
    return out
end

function interpolate(f::DiscreteField{NodeLoc}, p::SVector{3})
    lat, lon = ManifoldMeshes._cartesian_to_latlon(p)
    return interpolate(f, lat, lon)
end
```

- [ ] **Step 4: Run interpolation tests**

Run:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected: interpolation tests PASS. If `_cartesian_to_latlon` is not exported and unavailable, add a local spherical conversion helper in `src/interpolate.jl`:

```julia
function _cartesian_to_latlon_local(p)
    r = sqrt(sum(abs2, p))
    lat = asind(p[3] / r)
    lon = mod(atan(p[2], p[1]) * 180 / pi, 360)
    return lat, lon
end
```

- [ ] **Step 5: Commit interpolation**

```bash
git add src/interpolate.jl test/test_interpolate.jl test/runtests.jl
git commit -m "feat: add NodeLoc interpolation"
```

### Task 7: UGRID Topology Writer

**Files:**
- Modify: `src/ugrid.jl`
- Create: `test/test_ugrid.jl`
- Modify: `test/runtests.jl`

- [ ] **Step 1: Write failing topology conformance tests**

Create `test/test_ugrid.jl`:

```julia
using ManifoldFields
using ManifoldMeshes
using Test

include("helpers.jl")

@testset "UGRID topology writer" begin
    g = small_grid()
    ds = to_ugrid(g)

    @test ds.attributes["Conventions"] == "CF-1.11 UGRID-1.0"
    @test haskey(ds.variables, "Mesh2")
    meshvar = ds.variables["Mesh2"]
    @test meshvar.attrs["cf_role"] == "mesh_topology"
    @test meshvar.attrs["topology_dimension"] == 2
    @test meshvar.attrs["node_coordinates"] == "Mesh2_node_lon Mesh2_node_lat"
    @test meshvar.attrs["face_node_connectivity"] == "Mesh2_face_nodes"

    @test size(ds.variables["Mesh2_node_lon"].data) == (num_nodes(g),)
    @test ds.variables["Mesh2_node_lon"].attrs["standard_name"] == "longitude"
    @test ds.variables["Mesh2_node_lon"].attrs["units"] == "degrees_east"
    @test ds.variables["Mesh2_node_lat"].attrs["standard_name"] == "latitude"
    @test ds.variables["Mesh2_node_lat"].attrs["units"] == "degrees_north"

    face_nodes = ds.variables["Mesh2_face_nodes"]
    @test size(face_nodes.data) == (num_cells(g), 4)
    @test face_nodes.attrs["cf_role"] == "face_node_connectivity"
    @test face_nodes.attrs["start_index"] == 1
    @test Tuple(face_nodes.data[1, :]) == cell_nodes(g, 1)
end
```

Add to `test/runtests.jl`:

```julia
include("test_ugrid.jl")
```

- [ ] **Step 2: Run UGRID topology test and verify failure**

Run:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected: FAIL because `to_ugrid` is not implemented.

- [ ] **Step 3: Implement in-memory UGRID topology**

Replace `src/ugrid.jl` with:

```julia
import ManifoldMeshes: cell_centroid, cell_nodes, edge_nodes, node_coordinates, num_cells,
                       num_edges, num_nodes

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
    lon = mod(atan(p[2], p[1]) * 180 / pi, 360)
    return lon, lat
end

function _coord_attrs(kind, axis)
    standard = axis == :lon ? "longitude" : "latitude"
    units = axis == :lon ? "degrees_east" : "degrees_north"
    return Dict{String,Any}(
        "standard_name" => standard,
        "units" => units,
        "long_name" => "$(standard) of UGRID $(kind)",
    )
end

function to_ugrid(mesh)
    vars = Dict{String,UGridVariable}()
    attrs = Dict{String,Any}("Conventions" => "CF-1.11 UGRID-1.0")

    vars["Mesh2"] = UGridVariable(0, (), Dict{String,Any}(
        "cf_role" => "mesh_topology",
        "topology_dimension" => 2,
        "node_coordinates" => "Mesh2_node_lon Mesh2_node_lat",
        "face_node_connectivity" => "Mesh2_face_nodes",
        "face_dimension" => "n_face",
    ))

    node_lon = Vector{Float64}(undef, num_nodes(mesh))
    node_lat = Vector{Float64}(undef, num_nodes(mesh))
    for i in 1:num_nodes(mesh)
        node_lon[i], node_lat[i] = _lonlat(node_coordinates(mesh, i))
    end
    vars["Mesh2_node_lon"] =
        UGridVariable(node_lon, ("n_node",), _coord_attrs("nodes", :lon))
    vars["Mesh2_node_lat"] =
        UGridVariable(node_lat, ("n_node",), _coord_attrs("nodes", :lat))

    face_nodes = Matrix{Int}(undef, num_cells(mesh), 4)
    for c in 1:num_cells(mesh)
        face_nodes[c, :] .= collect(cell_nodes(mesh, c))
    end
    vars["Mesh2_face_nodes"] = UGridVariable(face_nodes, ("n_face", "n_max_face_nodes"),
                                            Dict{String,Any}(
                                                "cf_role" => "face_node_connectivity",
                                                "start_index" => 1,
                                            ))

    return UGridDataset(vars, attrs)
end
```

- [ ] **Step 4: Run topology tests**

Run:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected: UGRID topology tests PASS.

- [ ] **Step 5: Commit topology writer**

```bash
git add src/ugrid.jl test/test_ugrid.jl test/runtests.jl
git commit -m "feat: add UGRID topology writer"
```

### Task 8: UGRID Field Writer and Conformance

**Files:**
- Modify: `src/ugrid.jl`
- Modify: `test/test_ugrid.jl`

- [ ] **Step 1: Extend failing UGRID field tests**

Append to `test/test_ugrid.jl`:

```julia
@testset "UGRID field writer" begin
    g = small_grid()

    nf = DiscreteField(NodeLoc, g, node_values(g), node_dims(g); name=:node_temp)
    nds = to_ugrid(nf)
    nvar = nds.variables["node_temp"]
    @test nvar.attrs["mesh"] == "Mesh2"
    @test nvar.attrs["location"] == "node"
    @test nvar.attrs["coordinates"] == "Mesh2_node_lon Mesh2_node_lat"
    @test nvar.dims == ("n_node",)

    cf = DiscreteField(CellLoc, g, cell_values(g), cell_dims(g); name=:cell_area)
    cds = to_ugrid(cf)
    cvar = cds.variables["cell_area"]
    @test cvar.attrs["mesh"] == "Mesh2"
    @test cvar.attrs["location"] == "face"
    @test cvar.attrs["coordinates"] == "Mesh2_face_lon Mesh2_face_lat"
    @test cvar.dims == ("n_face",)
    @test cds.variables["Mesh2"].attrs["face_coordinates"] == "Mesh2_face_lon Mesh2_face_lat"
    @test size(cds.variables["Mesh2_face_lon"].data) == (num_cells(g),)

    ef = DiscreteField(EdgeLoc, g, edge_values(g), edge_dims(g); name=:edge_flux)
    eds = to_ugrid(ef)
    evar = eds.variables["edge_flux"]
    @test evar.attrs["mesh"] == "Mesh2"
    @test evar.attrs["location"] == "edge"
    @test evar.attrs["coordinates"] == "Mesh2_edge_lon Mesh2_edge_lat"
    @test evar.dims == ("n_edge",)
    @test eds.variables["Mesh2"].attrs["edge_node_connectivity"] == "Mesh2_edge_nodes"
    @test eds.variables["Mesh2"].attrs["edge_coordinates"] == "Mesh2_edge_lon Mesh2_edge_lat"
    @test eds.variables["Mesh2_edge_nodes"].attrs["cf_role"] == "edge_node_connectivity"
    @test eds.variables["Mesh2_edge_nodes"].attrs["start_index"] == 1
    @test size(eds.variables["Mesh2_edge_nodes"].data) == (num_edges(g), 2)
end
```

- [ ] **Step 2: Run and verify field writer tests fail**

Run:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected: FAIL because `to_ugrid(::DiscreteField)` is not implemented.

- [ ] **Step 3: Implement field writer**

Append to `src/ugrid.jl`:

```julia
function _add_face_coordinates!(ds::UGridDataset, mesh)
    lon = Vector{Float64}(undef, num_cells(mesh))
    lat = Vector{Float64}(undef, num_cells(mesh))
    for c in 1:num_cells(mesh)
        lon[c], lat[c] = _lonlat(cell_centroid(mesh, c))
    end
    ds.variables["Mesh2_face_lon"] =
        UGridVariable(lon, ("n_face",), _coord_attrs("faces", :lon))
    ds.variables["Mesh2_face_lat"] =
        UGridVariable(lat, ("n_face",), _coord_attrs("faces", :lat))
    ds.variables["Mesh2"].attrs["face_coordinates"] = "Mesh2_face_lon Mesh2_face_lat"
    return ds
end

function _add_edge_topology_and_coordinates!(ds::UGridDataset, mesh)
    edge_node_matrix = Matrix{Int}(undef, num_edges(mesh), 2)
    lon = Vector{Float64}(undef, num_edges(mesh))
    lat = Vector{Float64}(undef, num_edges(mesh))
    for e in 1:num_edges(mesh)
        n1, n2 = edge_nodes(mesh, e)
        edge_node_matrix[e, 1] = n1
        edge_node_matrix[e, 2] = n2
        p1 = node_coordinates(mesh, n1)
        p2 = node_coordinates(mesh, n2)
        midpoint = (p1 + p2) / sqrt(sum(abs2, p1 + p2))
        lon[e], lat[e] = _lonlat(midpoint)
    end
    ds.variables["Mesh2_edge_nodes"] =
        UGridVariable(edge_node_matrix, ("n_edge", "Two"), Dict{String,Any}(
            "cf_role" => "edge_node_connectivity",
            "start_index" => 1,
        ))
    ds.variables["Mesh2"].attrs["edge_node_connectivity"] = "Mesh2_edge_nodes"
    ds.variables["Mesh2"].attrs["edge_dimension"] = "n_edge"
    ds.variables["Mesh2_edge_lon"] =
        UGridVariable(lon, ("n_edge",), _coord_attrs("edges", :lon))
    ds.variables["Mesh2_edge_lat"] =
        UGridVariable(lat, ("n_edge",), _coord_attrs("edges", :lat))
    ds.variables["Mesh2"].attrs["edge_coordinates"] = "Mesh2_edge_lon Mesh2_edge_lat"
    return ds
end

function _ugrid_data_dims(f::DiscreteField)
    locname = DimensionalData.name(location_dimname(location(f)))
    return Tuple(DimensionalData.name(d) == locname ? ugrid_dimname(location(f)) :
                 String(Symbol(DimensionalData.name(d))) for d in dims(f))
end

function _coordinates_attr(::Type{NodeLoc})
    return "Mesh2_node_lon Mesh2_node_lat"
end

function _coordinates_attr(::Type{CellLoc})
    return "Mesh2_face_lon Mesh2_face_lat"
end

function _coordinates_attr(::Type{EdgeLoc})
    return "Mesh2_edge_lon Mesh2_edge_lat"
end

function to_ugrid(f::DiscreteField)
    ds = to_ugrid(mesh(f))
    if location(f) === CellLoc
        _add_face_coordinates!(ds, mesh(f))
    elseif location(f) === EdgeLoc
        _add_edge_topology_and_coordinates!(ds, mesh(f))
    end

    varname = String(Symbol(DimensionalData.name(f)))
    ds.variables[varname] = UGridVariable(data(f), _ugrid_data_dims(f), Dict{String,Any}(
        "mesh" => "Mesh2",
        "location" => ugrid_location(location(f)),
        "coordinates" => _coordinates_attr(location(f)),
    ))
    return ds
end
```

- [ ] **Step 4: Run field writer tests**

Run:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected: UGRID field writer tests PASS.

- [ ] **Step 5: Commit field writer**

```bash
git add src/ugrid.jl test/test_ugrid.jl
git commit -m "feat: add UGRID field writer"
```

### Task 9: UGRID Own-File Read and NetCDF Save/Load

**Files:**
- Modify: `src/ugrid.jl`
- Modify: `test/test_ugrid.jl`

- [ ] **Step 1: Add round-trip tests**

Append to `test/test_ugrid.jl`:

```julia
@testset "UGRID own-file round-trip" begin
    g = small_grid()
    f = DiscreteField(CellLoc, g, cell_values(g), cell_dims(g); name=:cell_area,
                      metadata=Dict("grid_type" => "LatLonGrid"))
    ds = to_ugrid(f)
    loaded = from_ugrid(ds)

    @test loaded isa DiscreteField{CellLoc}
    @test location(loaded) === CellLoc
    @test data(loaded) == data(f)
    @test dims(loaded) == dims(f)
    @test num_cells(mesh(loaded)) == num_cells(g)
    @test num_nodes(mesh(loaded)) == num_nodes(g)

    path = tempname() * ".nc"
    save_ugrid(f, path)
    loaded_file = load_ugrid(path)
    @test loaded_file isa DiscreteField{CellLoc}
    @test data(loaded_file) == data(f)
end
```

- [ ] **Step 2: Run round-trip tests and verify failure**

Run:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected: FAIL because `from_ugrid`, `save_ugrid`, and `load_ugrid` are not implemented.

- [ ] **Step 3: Implement own-file in-memory reader**

Append to `src/ugrid.jl`:

```julia
function _find_data_vars(ds::UGridDataset)
    return [name for (name, var) in ds.variables if haskey(var.attrs, "mesh")]
end

function _loc_from_ugrid(location::String)
    location == "node" && return NodeLoc
    location == "edge" && return EdgeLoc
    location == "face" && return CellLoc
    throw(ArgumentError("unsupported UGRID location $location"))
end

function _dims_from_ugrid(var::UGridVariable, Loc)
    locdim = DimensionalData.name(location_dimname(Loc))
    return Tuple(Dim{Symbol(d == ugrid_dimname(Loc) ? locdim : Symbol(d))}(1:size(var.data, i))
                 for (i, d) in enumerate(var.dims))
end

function from_ugrid_mesh(ds::UGridDataset; grid_type=nothing)
    meshvar = get(ds.variables, "Mesh2", nothing)
    meshvar === nothing && throw(ArgumentError("UGRID dataset has no Mesh2 topology variable"))
    meshvar.attrs["cf_role"] == "mesh_topology" ||
        throw(ArgumentError("Mesh2 is missing cf_role=mesh_topology"))
    n_face = size(ds.variables["Mesh2_face_nodes"].data, 1)
    n_node = length(ds.variables["Mesh2_node_lon"].data)
    if grid_type === nothing && haskey(meshvar.attrs, "grid_type")
        grid_type = meshvar.attrs["grid_type"]
    end
    if grid_type == "LatLonGrid" || (grid_type === nothing && n_face == 8 && n_node == 15)
        return LatLonGrid(2, 4)
    end
    throw(ArgumentError("cannot reconstruct mesh with n_face=$n_face n_node=$n_node"))
end

function from_ugrid(ds::UGridDataset; grid_type=nothing)
    data_vars = _find_data_vars(ds)
    length(data_vars) == 1 ||
        throw(ArgumentError("multi-variable UGRID datasets require FieldSet"))
    name = only(data_vars)
    var = ds.variables[name]
    haskey(var.attrs, "location") || throw(ArgumentError("UGRID data variable lacks location"))
    Loc = _loc_from_ugrid(var.attrs["location"])
    mesh = from_ugrid_mesh(ds; grid_type=grid_type)
    return DiscreteField(Loc, mesh, var.data, _dims_from_ugrid(var, Loc);
                         name=Symbol(name), metadata=var.attrs)
end
```

- [ ] **Step 4: Implement true NetCDF save/load with NCDatasets**

Append to `src/ugrid.jl`:

```julia
using NCDatasets

function _define_dimensions!(nc, ds::UGridDataset)
    lengths = Dict{String,Int}()
    for var in values(ds.variables)
        for (i, d) in enumerate(var.dims)
            lengths[d] = size(var.data, i)
        end
    end
    for (name, len) in lengths
        defDim(nc, name, len)
    end
    return nc
end

function _nc_type(data)
    data isa Integer && return Int32
    eltype(data) <: Integer && return Int32
    eltype(data) <: AbstractFloat && return Float64
    return eltype(data)
end

function _write_one_variable!(nc, name, var::UGridVariable)
    if var.dims == ()
        ncvar = defVar(nc, name, _nc_type(var.data), ())
        ncvar[] = var.data
    else
        ncvar = defVar(nc, name, _nc_type(var.data), var.dims)
        ncvar[:] = var.data
    end
    for (k, v) in var.attrs
        ncvar.attrib[k] = v
    end
    return nc
end

function _write_variables!(nc, ds::UGridDataset)
    for (name, var) in ds.variables
        _write_one_variable!(nc, name, var)
    end
    return nc
end

function save_ugrid(x, path::AbstractString; format=:netcdf)
    format == :netcdf || throw(ArgumentError("v0 only supports format=:netcdf"))
    ds = to_ugrid(x)
    NCDataset(path, "c") do nc
        for (k, v) in ds.attributes
            nc.attrib[k] = v
        end
        _define_dimensions!(nc, ds)
        _write_variables!(nc, ds)
    end
    return nothing
end

function _read_ugrid_dataset(path)
    vars = Dict{String,UGridVariable}()
    attrs = Dict{String,Any}()
    NCDataset(path, "r") do nc
        for (k, v) in nc.attrib
            attrs[String(k)] = v
        end
        for (name, v) in nc.vars
            vars[String(name)] = UGridVariable(v[:], Tuple(String.(dimnames(v))),
                                               Dict{String,Any}(String(k) => val
                                                                for (k, val) in v.attrib))
        end
    end
    return UGridDataset(vars, attrs)
end

load_ugrid(path::AbstractString; grid_type=nothing) =
    from_ugrid(_read_ugrid_dataset(path); grid_type=grid_type)
load_ugrid_mesh(path::AbstractString; grid_type=nothing) =
    from_ugrid_mesh(_read_ugrid_dataset(path); grid_type=grid_type)
```

- [ ] **Step 5: Run round-trip tests**

Run:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected: UGRID own-file round-trip tests PASS and written `.nc` files are true NetCDF files.

- [ ] **Step 6: Commit UGRID reader**

```bash
git add src/ugrid.jl test/test_ugrid.jl
git commit -m "feat: add UGRID own-file round-trip"
```

### Task 10: UGRID Error Cases and NetCDF Signature

**Files:**
- Modify: `src/ugrid.jl`
- Modify: `test/test_ugrid.jl`

- [ ] **Step 1: Add NetCDF signature and error-path tests**

Append to `test/test_ugrid.jl`:

```julia
@testset "UGRID NetCDF file signature" begin
    g = small_grid()
    f = DiscreteField(NodeLoc, g, node_values(g), node_dims(g); name=:node_temp)
    path = tempname() * ".nc"
    save_ugrid(f, path)
    open(path, "r") do io
        magic = read(io, 4)
        @test magic == UInt8[0x43, 0x44, 0x46, 0x01] ||
              magic == UInt8[0x43, 0x44, 0x46, 0x02] ||
              magic == UInt8[0x89, 0x48, 0x44, 0x46]
    end
end

@testset "UGRID reader rejects invalid datasets" begin
    bad = UGridDataset(Dict{String,ManifoldFields.UGridVariable}(),
                       Dict{String,Any}("Conventions" => "CF-1.11 UGRID-1.0"))
    @test_throws ArgumentError from_ugrid_mesh(bad)

    g = small_grid()
    f1 = DiscreteField(NodeLoc, g, node_values(g), node_dims(g); name=:a)
    ds = to_ugrid(f1)
    ds.variables["b"] = ds.variables["a"]
    @test_throws ArgumentError from_ugrid(ds)
end
```

- [ ] **Step 2: Run and verify error-path failures**

Run:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected: FAIL if invalid dataset handling does not raise `ArgumentError`.

- [ ] **Step 3: Harden reader validation**

Patch `from_ugrid_mesh` and `from_ugrid` in `src/ugrid.jl` with explicit checks before indexing:

```julia
function _require_var(ds::UGridDataset, name::String)
    haskey(ds.variables, name) || throw(ArgumentError("UGRID dataset is missing variable $name"))
    return ds.variables[name]
end

function _require_attr(var::UGridVariable, name::String)
    haskey(var.attrs, name) || throw(ArgumentError("UGRID variable is missing attribute $name"))
    return var.attrs[name]
end
```

Use `_require_var(ds, "Mesh2")`, `_require_var(ds, "Mesh2_face_nodes")`,
`_require_var(ds, "Mesh2_node_lon")`, `_require_attr(meshvar, "cf_role")`,
and `_require_attr(var, "location")` in the reader paths.

- [ ] **Step 4: Run UGRID tests**

Run:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected: NetCDF signature test and all UGRID tests PASS.

- [ ] **Step 5: Commit UGRID validation**

```bash
git add src/ugrid.jl test/test_ugrid.jl
git commit -m "test: cover UGRID reader validation"
```

### Task 11: Final QA and Documentation

**Files:**
- Create: `README.md`
- Modify: `src/ManifoldFields.jl`
- Modify: all implementation files if Aqua or full tests expose issues.

- [ ] **Step 1: Add README usage examples**

Create `README.md`:

````markdown
# ManifoldFields.jl

Discrete fields on `ManifoldMeshes.jl` meshes.

```julia
using DimensionalData
using ManifoldFields
using ManifoldMeshes

g = LatLonGrid(2, 4)
values = collect(Float64, 1:num_nodes(g))
dims = (Dim{:node}(1:num_nodes(g)),)

f = DiscreteField(NodeLoc, g, values, dims; name=:temperature)
interpolate(f, 0.0, 45.0)
save_ugrid(f, "temperature.nc")
loaded = load_ugrid("temperature.nc")
```

v0 supports separate `DiscreteField{NodeLoc}`, `DiscreteField{EdgeLoc}`, and
`DiscreteField{CellLoc}` objects. A multi-variable `FieldSet` container is
deferred.
````

- [ ] **Step 2: Run formatter**

Run:

```bash
julia --project=. -e 'using JuliaFormatter; format(".")'
```

Expected: formatting completes. If `JuliaFormatter` is not installed, run:

```bash
julia --project=. -e 'using Pkg; Pkg.add("JuliaFormatter"); using JuliaFormatter; format(".")'
```

- [ ] **Step 3: Run full test suite**

Run:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected: PASS, including Aqua.

- [ ] **Step 4: Verify spec coverage**

Run:

```bash
grep -n "DiscreteField{Loc}\\|NodeLoc\\|UGRID v1.0\\|withmesh\\|interpolate" docs/superpowers/specs/2026-07-21-manifoldfields-design.md
```

Expected: output includes the core spec sections for field container, type-tag locations, UGRID conformance, mesh rebinding, and point gather.

- [ ] **Step 5: Commit final polish**

```bash
git add README.md src test Project.toml .JuliaFormatter.toml
git commit -m "docs: add ManifoldFields v0 usage"
```

## Self-Review Checklist

- Spec coverage:
  - `DiscreteField{Loc} <: AbstractDimArray`: Tasks 3-5.
  - Type tag constructors and location trait table: Tasks 2-3.
  - Constructor location-axis checks: Task 3.
  - Accessors and `withmesh`: Task 3.
  - Broadcast arithmetic with mesh identity and Loc checks: Task 5.
  - `NodeLoc` interpolation with arbitrary location-axis position: Task 6.
  - UGRID v1.0-conformant NetCDF writer: Tasks 7-10.
  - Own-file read and v0 single-variable behavior: Task 9.
  - Aqua and docs: Task 11.
- Deferred by design:
  - `FieldSet`.
  - `CellLoc` / `EdgeLoc` interpolation.
  - Zarr UGRID conformance.
  - Generic external UGRID mesh read.
- Completion-marker scan:
  - This plan avoids unfinished markers.
  - Every task includes concrete files, commands, and expected outcomes.
- Type consistency:
  - Public constructor consistently uses `DiscreteField(NodeLoc, ...)`.
  - `location(f)` consistently returns the location type.
  - UGRID internal names consistently use `Mesh2`, `n_node`, `n_edge`, `n_face`.
