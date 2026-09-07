# ManifoldFields.jl × UXarray — a paired tutorial

Side-by-side notebooks for UXarray (Python) users adopting ManifoldFields.jl
(Julia). Each chapter is a pair — `julia/NN-name.jl` (Pluto) and
`python/NN-name.py` (marimo) — with identical section titles and one-to-one
cells, so you can run both and compare. Every "✦ equivalent" snippet in one
notebook runs for real in its twin.

## Concept map

| UXarray (Python) | ManifoldFields.jl (Julia) |
|---|---|
| `ux.open_dataset(grid, data)` | `load_ugrid(path)` → `FieldSet` |
| `UxDataArray` | `DiscreteField` (location-aware: `NodeLoc` / `CellLoc` / `EdgeLoc`) |
| `UxDataset` / `.uxgrid` | `FieldSet` / `mesh(fs)` |
| any UGRID topology (MPAS, ICON, SCRIP, Exodus, …) | four manifold-native grid types: `LatLonGrid`, `CubedSphereGrid`, `ReducedGaussianGrid`, `HEALPixGrid` |
| `uxgrid.to_xarray(grid_format="ugrid")` | `to_ugrid(mesh)` / `save_ugrid(fs, path)` |
| `da.mean(dim=...)` (keeps length-1 dims) | `mean(fs; dims = ..., keepdims = false)` (dims dropped by default) |
| `get_ball_tree()` / `remap.nearest_neighbor(...)` | `interpolate(field, lat, lon)` (true interpolation) |

See [CAPABILITY_MAP.md](CAPABILITY_MAP.md) for full coverage and
[DATA.md](DATA.md) for the datasets and measured interop results.

## Running the notebooks

### Julia (Pluto)

```bash
julia -e 'import Pkg; Pkg.add("Pluto")'
julia -e 'import Pluto; Pluto.run()'
```

Open `examples/julia/01-introduction.jl` (this README's directory). The
committed `julia/Project.toml` pins the notebook environment: `ManifoldFields`
resolves to the repository you are in via `[sources]`, `ManifoldMeshes` is
fetched from its GitHub URL.

### Python (marimo)

```bash
python -m venv .venv && source .venv/bin/activate
pip install -r python/requirements.txt
marimo edit python/01-introduction.py
```

### Data

Notebooks download the shared datasets on first run (≈2.4 MB, sha256-verified)
into the gitignored `data/` directory. Julia resolves `Artifacts.toml` with
`Pkg.ensure_artifact_installed`; Python reads the same manifest through
`python/_data.py`. Nothing is bundled with `Pkg.add("ManifoldFields")`.

## Why paired notebooks?

Pluto and marimo are both reactive: change a cell and everything downstream
recomputes — the exploration loop feels the same in both languages.
