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

## [0.1.0] - 2026-07-23

### Added

- Initial `DiscreteField` container for `ManifoldMeshes.jl` meshes.
- Location tags for node, edge, and cell fields.
- Node-location interpolation.
- UGRID own-file NetCDF round-trip support.
- Broadcast support preserving field metadata.

[unreleased]: https://github.com/VANvonZHANG/ManifoldFields.jl/compare/v0.1.0...HEAD
[0.1.0]: https://github.com/VANvonZHANG/ManifoldFields.jl/releases/tag/v0.1.0

