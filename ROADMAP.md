# Roadmap

This file records the next engineering steps after the initial
`DiscreteField`-focused package surface.

## Current Decision

Use `ManifoldFields.jl` as a private/development package for now. Do not submit
to the Julia General Registry until the `FieldSet` API and UGRID multi-field IO
contracts have stabilized, local tests and GitHub Actions are green, and release
notes are ready in `CHANGELOG.md`.

Registration prerequisites:

- GitHub Actions CI is green on Julia 1.10 and 1.11.
- `LICENSE` exists.
- `README.md` and Documenter docs describe the public API.
- `Project.toml` version is bumped according to Julia `0.x` SemVer practice.
- JuliaRegistrator is installed on the repository.
- TagBot and `DOCUMENTER_KEY`/deploy keys are configured if docs deployment is
  expected.

## Priority Order

1. CI and docs cleanup.
2. `FieldSet`.
3. UGRID multi-field IO.
4. Batch interpolation.
5. Remapping.
6. Lazy/out-of-core hardening.

## Compatibility Audit

Current direct dependency bounds are intentionally conservative.

- `julia = "1.10"` is appropriate for the planned GitHub Actions matrix because
  Julia compat `"1.10"` allows Julia 1.10 and later 1.x releases, including
  Julia 1.11.
- `NCDatasets = "0.14"` is not over-tight; it permits compatible 0.14 patch
  releases.
- `DimensionalData = "0.29"` and `YAXArrays = "0.6"` should stay pinned to
  those minor series for now. A temporary solve against `DimensionalData 0.30.1`
  and `YAXArrays 0.7.7` reaches the test suite but fails in broadcast with
  `MethodError: no method matching format(::Nothing, ::OffsetArrays.OffsetVector...)`.

Before widening those two bounds, fix or adapt `src/broadcast.jl` against the
new `DimensionalData` broadcast behavior and rerun the full test suite.

## FieldSet

Add `FieldSet(mesh, fields)` as the multi-variable container for fields that
share the same mesh. A typical constructor should support
`FieldSet(mesh, Dict(:u => u, :v => v))`.

Design constraints:

- Values may be `NodeLoc`, `EdgeLoc`, or `CellLoc`.
- Every member field must reference the same mesh as the set.
- UGRID multi-variable read/write should attach to `FieldSet`, not overload
  `DiscreteField` with multi-variable semantics.
- Keep `DiscreteField` as the single-variable building block.

## UGRID IO

The current IO scope is own-file round trips. The next reader/writer pass should
handle more external UGRID files while staying aligned with concrete mesh types
available in `ManifoldMeshes.jl`.

Planned coverage:

- Multiple data variables in one file.
- Different dimension orders.
- External mesh topology validation.
- More complete CF/UGRID attributes.

Do not broaden to arbitrary unstructured mesh support until `ManifoldMeshes.jl`
has a corresponding concrete mesh representation.

## Interpolation

Current interpolation supports `NodeLoc`.

Planned extensions:

- Batch node interpolation for many `(lat, lon)` query points in one call.
- `CellLoc` nearest or containing-cell gather.
- Defer `EdgeLoc` interpolation until there is a clear physical use case.

## Remapping

Introduce a remapping interface, either as a module or a package such as
`ManifoldRegrid.jl`.

Target API:

```julia
remap(field, dest_mesh; method = :bilinear)
```

Implementation should reuse existing cell location, interpolation weight, and
`DiscreteField` machinery.

## Lazy and Out-of-Core

Correctness remains the priority. Before using large atmospheric chemistry or
HEMCO datasets heavily, verify:

- `YAXArrays` lazy data remains lazy through field construction.
- Broadcast does not unexpectedly materialize large arrays.
- UGRID save can write large arrays in chunks.
