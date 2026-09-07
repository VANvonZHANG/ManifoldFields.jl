### A Pluto.jl notebook ###
# ╔═╡ 00000000-0000-4000-8000-000000000001
md"""
# 01 · Introduction

The same dataset, opened and explored in ManifoldFields.jl (Julia) and UXarray
(Python). Twin notebook: `python/01-introduction.py`.
"""

# ╔═╡ 00000000-0000-4000-8000-00000000000d
begin
    import Pkg
    # Use the environment pinned by examples/julia/Project.toml (the directory of
    # this notebook) instead of Pluto's built-in package manager: it contains
    # unregistered [sources] packages that Pluto cannot resolve on its own.
    Pkg.activate(@__DIR__; io = devnull)
    # NB: test loadability, not Pkg.project().dependencies — the committed Project.toml
    # already lists the [deps], so only the Manifest (gitignored) tells fresh clones apart.
    if isnothing(Base.find_package("ManifoldFields")) ||
       isnothing(Base.find_package("ManifoldMeshes"))
        # Bootstrap on a fresh clone (Pluto's nbpkg cannot use [sources] on Julia 1.10):
        Pkg.develop(path = joinpath(@__DIR__, "..", ".."); io = devnull)
        Pkg.add(url = "https://github.com/VANvonZHANG/ManifoldMeshes.jl"; io = devnull)
    end
end

# ╔═╡ 00000000-0000-4000-8000-000000000002
begin
    using ManifoldFields
    using ManifoldMeshes
    using DimensionalData
    using CairoMakie
end

# ╔═╡ 00000000-0000-4000-8000-000000000003
DATA = Pkg.ensure_artifact_installed("example-data", joinpath(@__DIR__, "..", "Artifacts.toml"))

# ╔═╡ 00000000-0000-4000-8000-000000000004
md"""
✦ **Python equivalent** (runs in the twin notebook):

```python
uxds = ux.open_dataset(f"{DATA}/healpix_analytic.nc")
```
"""

# ╔═╡ 00000000-0000-4000-8000-000000000005
fs = load_ugrid(joinpath(DATA, "healpix_analytic.nc"))

# ╔═╡ 00000000-0000-4000-8000-000000000006
md"""
`fs` is a `FieldSet` — multiple variables on one shared mesh. UXarray's twin
concept is `UxDataset` with its `.uxgrid`.
"""

# ╔═╡ 00000000-0000-4000-8000-000000000007
field_names(fs)

# ╔═╡ 00000000-0000-4000-8000-000000000008
f = fs[:psi]

# ╔═╡ 00000000-0000-4000-8000-000000000009
begin
    v = data(f)
    extrema(v)
end

# ╔═╡ 00000000-0000-4000-8000-00000000000a
md"""
`psi` is an analytic field, `psi(lat, lon) = cos(lat)·sin(lon) + 0.3·sin(2·lat)`,
stored at cell centers of a HEALPix grid with nside = 32 (12 288 cells) — written
by `save_ugrid` in Julia, read natively by UXarray.
"""

# ╔═╡ 00000000-0000-4000-8000-00000000000b
ManifoldMeshes.plot_mesh_filled(ManifoldFields.mesh(fs); color_by = i -> v[i])

# ╔═╡ 00000000-0000-4000-8000-00000000000c
md"""
## Recap

| here (ManifoldFields) | twin (UXarray) |
|---|---|
| `load_ugrid(path)` | `ux.open_dataset(path)` |
| `FieldSet` | `UxDataset` |
| `DiscreteField` | `UxDataArray` |
| `ManifoldFields.mesh(fs)` | `uxds.uxgrid` |

Next: **02 · Grids & UGRID IO**.
"""
