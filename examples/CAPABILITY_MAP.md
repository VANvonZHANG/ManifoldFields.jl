# UXarray ↔ ManifoldFields capability map

Rows follow the UXarray user guide (https://uxarray.readthedocs.io/).
Remapping and batch-interpolation roadmap items live in [../ROADMAP.md](../ROADMAP.md).

| UXarray capability | UXarray API | ManifoldFields today | Julia-side status |
|---|---|---|---|
| Open MPAS / ICON / SCRIP / Exodus / … | `ux.open_dataset(grid, data)` | manifold-native grids instead (`LatLonGrid`, `CubedSphereGrid`, `ReducedGaussianGrid`, `HEALPixGrid`) | by design |
| UGRID read / write | `uxgrid.to_xarray(grid_format="ugrid")` | `load_ugrid` / `save_ugrid` / `to_ugrid`, round-trips all four grid types | ✅ shipped |
| HEALPix | `ux.UxDataset.from_healpix(path)` / `ux.Grid.from_healpix(zoom=...)` | `HEALPixGrid(nside, ordering)` | ✅ shipped |
| Subsetting | `uxds.subset.bounding_box(...)` | — | roadmap |
| Cross-sections | `uxds.cross_section(...)` | — | roadmap |
| Zonal / azimuthal means | `da.zonal_mean(...)` / `da.azimuthal_mean(...)` | — | roadmap (conservative ops) |
| Weighted / plain means | `da.mean(dim=...)` / `da.weighted(...)` | `mean(fs; dims, keepdims)` | ✅ shipped |
| Remapping | `da.remap.nearest_neighbor(...)` / `.bilinear(...)` / `.inverse_distance_weighted(...)` | — | roadmap (NodeLoc bilinear prototype) |
| Point interpolation | `uxgrid.get_ball_tree()` (nearest) | `interpolate(field, lat, lon)` (true interpolation) | ✅ shipped; batched queries roadmap |
| Face areas / integrals | `uxgrid.calculate_total_face_area(...)` / `da.integrate(...)` | `cell_volume` (ManifoldMeshes) | ✅ shipped |
| Vector calculus | `da.gradient(...)` etc. | — | roadmap |
| Dual mesh | `uxgrid.get_dual()` | dual-mesh support in ManifoldMeshes | partial |
| Topological aggregation | `da.topological_mean(...)` | — | roadmap |
| Plotting | `da.plot()`, HoloViz | `ManifoldMeshes.plot_mesh` / `plot_mesh_filled` (CairoMakie) | ✅ shipped |

> Note: chapter 04 measured an upstream bilinear-weights transposition bug in ManifoldMeshes 0.6.0 (LatLon/ReducedGaussian grids) — point interpolation works but with elevated error until the upstream fix; see `examples/DATA.md` and chapter 04.
