### A Pluto.jl notebook ###
# ╔═╡ 00000000-0000-4000-8000-000000000011
md"""
# 02 · Grids & UGRID IO

ManifoldMeshes ships four manifold-native S² grid types. UXarray instead
represents *any* UGRID topology (MPAS, ICON, SCRIP, …) — this chapter shows
both sides of that boundary. Twin notebook: `python/02-grids-and-ugrid-io.py`.
"""

# ╔═╡ 00000000-0000-4000-8000-000000000012
begin
    import Pkg
    Pkg.activate(@__DIR__; io = devnull)
    if isnothing(Base.find_package("ManifoldFields")) || isnothing(Base.find_package("ManifoldMeshes"))
        # Bootstrap on a fresh clone (Pluto's nbpkg cannot use [sources] on Julia 1.10):
        Pkg.develop(path = joinpath(@__DIR__, "..", ".."); io = devnull)
        Pkg.add(url = "https://github.com/VANvonZHANG/ManifoldMeshes.jl"; io = devnull)
    end
end

# ╔═╡ 00000000-0000-4000-8000-000000000013
begin
    using ManifoldFields
    using ManifoldMeshes
    using DimensionalData
    using NCDatasets
end

# ╔═╡ 00000000-0000-4000-8000-000000000014
DATA = Pkg.ensure_artifact_installed("example-data", joinpath(@__DIR__, "..", "Artifacts.toml"))

# ╔═╡ 00000000-0000-4000-8000-000000000015
md"""
✦ **Python equivalent** — there is no grid zoo to tour; any UGRID/MPAS/… file
opens through one entry point:

```python
uxds   = ux.open_dataset(grid_path, data_path)  # MPAS pair
uxgrid = ux.open_grid(grid_path)                # topology only
```
"""

# ╔═╡ 00000000-0000-4000-8000-000000000016
begin
    meshes = (
        ("LatLonGrid", LatLonGrid(lat_edges = collect(-90.0:10.0:90.0),
                                  lon_edges = collect(0.0:10.0:360.0))),
        ("CubedSphereGrid", CubedSphereGrid(n = 16)),
        ("ReducedGaussianGrid", ReducedGaussianGrid(nlat = 64)),
        ("HEALPixGrid", HEALPixGrid(nside = 32)),
    )
    [(name, num_cells(m)) for (name, m) in meshes]
end

# ╔═╡ 00000000-0000-4000-8000-000000000017
md"""
## Round-trip on our own files

`save_ugrid` / `load_ugrid` round-trip any `FieldSet`; the file carries the
mesh metadata needed to reconstruct the exact grid type.
"""

# ╔═╡ 00000000-0000-4000-8000-000000000018
begin
    fs = load_ugrid(joinpath(DATA, "cubedsphere_analytic.nc"))
    save_ugrid(fs, joinpath(@__DIR__, "cubedsphere_copy.nc"))
    fs2 = load_ugrid(joinpath(@__DIR__, "cubedsphere_copy.nc"))
    isapprox(data(fs2[:psi]), data(fs[:psi]); atol = 1e-12)
end

# ╔═╡ 00000000-0000-4000-8000-000000000019
md"""
## The other direction: a foreign UGRID file

UXarray wrote `oQU480.ugrid.nc` from an MPAS mesh. The load stops at the very
first stage: ManifoldFields looks up the UGRID topology container under its
own writer name — `_require_var(ds, "Mesh2")` in `src/ugrid.jl` — while
uxarray names the container `grid_topology`. And even past that stage, an
arbitrary MPAS quadrilateral topology would match none of the four grid types
`load_ugrid` can reconstruct from own-file metadata. When a foreign file's
geometry *does* match, pass `load_ugrid(path; mesh = m)`.
"""

# ╔═╡ 00000000-0000-4000-8000-00000000001a
try
    load_ugrid(joinpath(DATA, "oQU480.ugrid.nc"))
catch e
    showerror(stdout, e)
end

# ╔═╡ 00000000-0000-4000-8000-00000000001b
md"""
## HEALPix values without UGRID topology

`psi_healpix.nc` is a bare value vector (12 288 cells, no coordinates, no
connectivity). UXarray reads it via `UxDataset.from_healpix`. On the Julia side we
wrap the values on a matching `HEALPixGrid` — the cell ordering that agrees
with UXarray's reading is `:nested` (see `DATA.md` for the measurement);
values correspond cell-by-cell by index, and the positional deviations on the
equatorial belt cells are documented in `DATA.md`.
"""

# ╔═╡ 00000000-0000-4000-8000-00000000001c
begin
    ds = Dataset(joinpath(DATA, "psi_healpix.nc"))
    psi_raw = Array(ds["psi"][:])
    close(ds)
    hp = HEALPixGrid(nside = 32, ordering = :nested)
    psi_hp = DiscreteField(CellLoc, hp, psi_raw, (Dim{:cell}(1:num_cells(hp)),); name = :psi)
    length(psi_raw) == num_cells(hp)
end

# ╔═╡ 00000000-0000-4000-8000-00000000001d
md"""
## Recap

Construct grids from parameters; persist any `FieldSet` as UGRID and read it
back exactly. Foreign topologies are UXarray's home turf — ManifoldFields
excels when the sphere grid itself is first-class. Next: **03 · Fields & FieldSet**.
"""
