### A Pluto.jl notebook ###
# ╔═╡ 00000000-0000-4000-8000-000000000021
md"""
# 03 · Fields & FieldSet

Values live at a *location* of the mesh — nodes, cells, or edges — and that
fact is part of the type. A `FieldSet` holds several fields on one shared
mesh, like a `UxDataset`. Twin notebook: `python/03-fields-and-fieldset.py`.
"""

# ╔═╡ 00000000-0000-4000-8000-000000000022
begin
    import Pkg
    Pkg.activate(@__DIR__; io = devnull)
    if isnothing(Base.find_package("ManifoldFields")) || isnothing(Base.find_package("ManifoldMeshes"))
        # Bootstrap on a fresh clone (Pluto's nbpkg cannot use [sources] on Julia 1.10):
        Pkg.develop(path = joinpath(@__DIR__, "..", ".."); io = devnull)
        Pkg.add(url = "https://github.com/VANvonZHANG/ManifoldMeshes.jl"; io = devnull)
    end
end

# ╔═╡ 00000000-0000-4000-8000-000000000023
begin
    using ManifoldFields
    using ManifoldMeshes
    using DimensionalData
    using Statistics
    using CairoMakie
end

# ╔═╡ 00000000-0000-4000-8000-000000000024
DATA = Pkg.ensure_artifact_installed("example-data", joinpath(@__DIR__, "..", "Artifacts.toml"))

# ╔═╡ 00000000-0000-4000-8000-000000000025
begin
    g = LatLonGrid(lat_edges = collect(-90.0:15.0:90.0), lon_edges = collect(0.0:15.0:360.0))
    (num_nodes(g), num_cells(g), num_edges(g))
end

# ╔═╡ 00000000-0000-4000-8000-000000000026
begin
    node_t = DiscreteField(NodeLoc, g, randn(num_nodes(g)),
                           (Dim{:node}(1:num_nodes(g)),); name = :node_t)
    cell_t = DiscreteField(CellLoc, g, randn(num_cells(g)),
                           (Dim{:cell}(1:num_cells(g)),); name = :cell_t)
end

# ╔═╡ 00000000-0000-4000-8000-000000000027
md"""
✦ **Python equivalent** — location is implied by the dimension name:

```python
uxds["psi"].dims  # ("n_face",) — cell data; node data would carry ("n_node",)
```

In Julia the location is a type parameter (`NodeLoc`, `CellLoc`, `EdgeLoc`),
so `interpolate` can accept only node fields, plotting can pick the right
coordinates, and mixed-location `FieldSet`s stay valid.
"""

# ╔═╡ 00000000-0000-4000-8000-000000000028
begin
    fs = FieldSet(g, :node_t => node_t, :cell_t => cell_t)
    (field_names(fs), location.(collect(fields(fs))))
end

# ╔═╡ 00000000-0000-4000-8000-000000000029
md"""
## Broadcasting keeps the mesh

Arithmetic broadcasts field-by-field and preserves the shared mesh — compare
UXarray's `2 * uxds + 1`.
"""

# ╔═╡ 00000000-0000-4000-8000-00000000002a
begin
    fs2 = 2 .* fs .+ 1
    ManifoldFields.mesh(fs2) === ManifoldFields.mesh(fs)
end

# ╔═╡ 00000000-0000-4000-8000-00000000002b
md"""
## Reductions

`mean` (also `sum`, `std`, …) reduces over `dims` — and **drops** the reduced
dimensions by default; `keepdims = true` retains them as length-1 dimensions.
xarray made the same choice (`keepdims=True`), so this default is one the
twins share — numpy is the odd one out, keeping length-1 axes. Where they
genuinely disagree: reducing a *location* dimension (`Dim{:node}`,
`Dim{:cell}`) throws, because the mesh extent is part of the field's
contract, while UXarray reduces `n_face` like any other dimension.
"""

# ╔═╡ 00000000-0000-4000-8000-00000000002c
mean(fs)  # NamedTuple of scalars, one entry per field

# ╔═╡ 00000000-0000-4000-8000-00000000002d
begin
    node_t4 = DiscreteField(NodeLoc, g, randn(4, num_nodes(g)),
                            (Dim{:time}(1:4), Dim{:node}(1:num_nodes(g))); name = :node_t)
    fs_t = FieldSet(g, :node_t => node_t4)
    dims(node_t4)
end

# ╔═╡ 00000000-0000-4000-8000-00000000002e
begin
    dropped = mean(fs_t; dims = Dim{:time})                  # Dim{:time} dropped
    kept    = mean(fs_t; dims = Dim{:time}, keepdims = true) # retained, length 1
    (dims(dropped[:node_t]), dims(kept[:node_t]))
end

# ╔═╡ 00000000-0000-4000-8000-00000000002f
try
    mean(fs; dims = Dim{:node}, keepdims = true)
catch e
    showerror(stdout, e)
end

# ╔═╡ 00000000-0000-4000-8000-000000000030
md"""
## On real data

The analytic HEALPix field from chapter 01, minus its mean, colored by cell:
"""

# ╔═╡ 00000000-0000-4000-8000-000000000031
begin
    fs_h = load_ugrid(joinpath(DATA, "healpix_analytic.nc"))
    anomaly = fs_h[:psi] .- mean(fs_h[:psi])
    v = data(anomaly)
    ManifoldMeshes.plot_mesh_filled(ManifoldFields.mesh(fs_h); color_by = i -> v[i])
end

# ╔═╡ 00000000-0000-4000-8000-000000000032
md"""
## Recap

Location-aware fields, one mesh shared by many variables, broadcasts that
preserve metadata, reductions that drop dims by default — and location dims
that refuse to be reduced at all. Next: **04 · Point Interpolation**.
"""

