# Changelog

All notable changes to this project will be documented in this file.

The format is based on Keep a Changelog, and this project follows Semantic
Versioning. While the package is in the `0.x` series, backward-compatible
features and fixes should use patch releases; breaking changes should use minor
releases.

## [Unreleased]

### Added

- Add GitHub Actions CI, documentation, formatting, CompatHelper, TagBot, and
  Dependabot configuration.
- Add Documenter.jl documentation skeleton.
- Add project roadmap for FieldSet, UGRID IO, interpolation, remapping, and
  lazy/out-of-core work.
- Add `FieldSet`: multi-field container of `DiscreteField`s sharing one mesh,
  subtyping `AbstractDimStack`; mixed `NodeLoc`/`EdgeLoc`/`CellLoc` members;
  `fs[:key]` materializes a `DiscreteField`; per-field dot-broadcast via
  `FieldSetStyle`; `merge` with mesh validation; dimensional slicing and
  reductions.
- Add UGRID multi-field IO: `to_ugrid(::FieldSet)`, `from_ugrid` returns a
  `FieldSet` for any number of data variables, `save_ugrid`/`load_ugrid`
  multi-field round-trips, full topology validation with `start_index`
  normalization, and external-file reading via `load_ugrid(path; mesh = ...)`.
- Add `fields(fs)` and `field_names(fs)` accessors.
- Reductions on `FieldSet` drop fully-reduced dimensions (unlike DimensionalData's length-1 retention).
- Own-file UGRID round-trips for `CubedSphereGrid`, `ReducedGaussianGrid`, and
  `HEALPixGrid` (constructor-parameter serialization including `rotation` for
  `CubedSphereGrid`/`HEALPixGrid`, validated against both connectivity and node
  coordinates).
- `keepdims` keyword for `FieldSet` reductions (retain length-1 reduced dims).
- `from_ugrid` reads UGRID global attributes back into FieldSet-level metadata.
- Add `examples/` — a paired Pluto/marimo tutorial for UXarray users, with a
  release-hosted example-data artifact (nothing bundled with `Pkg.add`).
- Foreign UGRID topologies now reconstruct generically as
  `ManifoldMeshes.UnstructuredMesh` (any cell arity ≥ 3 — quads, triangles,
  ICON hexagons) instead of erroring — `load_ugrid` is a general UGRID reader
  (requires ManifoldMeshes 0.7). Mixed-arity meshes round-trip with
  `_FillValue`-padded connectivity; parametric-grid output is unchanged.
  Degenerate connectivity keeps an actionable `mesh=` error.
  `load_ugrid(path; mesh = m)` still forces a specific grid.

### Changed

- `from_ugrid`/`load_ugrid` now return a `FieldSet` (previously a single
  `DiscreteField`; multi-variable datasets were rejected).

### Fixed

- `from_ugrid`/`load_ugrid` now discover the UGRID topology container via
  `cf_role = "mesh_topology"` instead of requiring the variable name `Mesh2`,
  resolve node coordinates and face-node connectivity through the container's
  own attributes, exclude coordinate/connectivity variables from data-variable
  discovery, and infer a data variable's location from its dimensions when the
  `location` attribute is absent. Third-party UGRID files (e.g. written by
  UXarray, which names the container `grid_topology`) now reach grid
  reconstruction; foreign geometries fail there with an actionable error
  pointing at `mesh =`, and UXarray round-trips of ManifoldFields-written files
  load natively (topology metadata survives the round trip).
- Node-major `(n_max_face_nodes, n_face)` face-node connectivity layouts (as
  emitted in this repo's earlier tutorial artifacts) are oriented on read by
  comparing the connectivity's second dimension to the declared
  `face_dimension`.
- `_FillValue` is now written before the data definition (and cast to the
  on-disk type), matching netCDF define-mode practice.
- `Base.cat` on `FieldSet`s now validates mesh identity (`===`) before
  concatenating.

## [0.1.0] - 2026-07-23

### Added

- Initial `DiscreteField` container for `ManifoldMeshes.jl` meshes.
- Location tags for node, edge, and cell fields.
- Node-location interpolation.
- UGRID own-file NetCDF round-trip support.
- Broadcast support preserving field metadata.

[unreleased]: https://github.com/VANvonZHANG/ManifoldFields.jl/compare/v0.1.0...HEAD
[0.1.0]: https://github.com/VANvonZHANG/ManifoldFields.jl/releases/tag/v0.1.0

