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
  `HEALPixGrid` (constructor-parameter serialization plus node-coordinate and
  connectivity validation; rotated `CubedSphereGrid` fails loudly — use
  `mesh=` injection).
- `keepdims` keyword for `FieldSet` reductions (retain length-1 reduced dims).
- `from_ugrid` reads UGRID global attributes back into FieldSet-level metadata.

### Changed

- `from_ugrid`/`load_ugrid` now return a `FieldSet` (previously a single
  `DiscreteField`; multi-variable datasets were rejected).

### Fixed

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

