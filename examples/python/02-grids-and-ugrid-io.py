import marimo

__generated_with = "0.9.0"
app = marimo.App(width="medium")


@app.cell
def _():
    import marimo as mo
    return (mo,)


@app.cell
def _(mo):
    mo.md(
        r"""
    # 02 · Grids & UGRID IO

    ManifoldMeshes ships four manifold-native S² grid types. UXarray instead
    represents *any* UGRID topology (MPAS, ICON, SCRIP, …) — this chapter shows
    both sides of that boundary. Twin notebook: `julia/02-grids-and-ugrid-io.jl`.
    """
    )
    return


@app.cell
def _():
    import sys
    from pathlib import Path

    sys.path.insert(0, str(Path(__file__).parent))
    from _data import ensure_data

    DATA = ensure_data()
    return DATA, Path


@app.cell
def _():
    import uxarray as ux
    return (ux,)


@app.cell
def _(mo):
    mo.md(
        r"""
    ## A file the Julia side wrote

    `healpix_analytic.nc` / `cubedsphere_analytic.nc` were produced by
    `ManifoldFields.save_ugrid`. UGRID is UXarray's native format, so they
    open directly.
    """
    )
    return


@app.cell
def _(DATA, ux):
    uxds = ux.open_dataset(f"{DATA}/cubedsphere_analytic.nc")
    uxds
    return (uxds,)


@app.cell
def _(mo):
    mo.md(
        r"""
    ✦ **Julia equivalent**:

    ```julia
    fs = load_ugrid(joinpath(DATA, "cubedsphere_analytic.nc"))
    fs[:psi]
    ```
    """
    )
    return


@app.cell
def _():
    import matplotlib.pyplot as plt
    return (plt,)


@app.cell
def _(Path, uxds):
    out_path = Path(__file__).parent / "cubedsphere_copy.nc"
    out = uxds.uxgrid.to_xarray()
    for var in uxds.data_vars:
        out[var] = uxds[var].variable
        out[var].attrs.update(uxds[var].attrs)
        out[var].attrs["mesh"] = "Mesh2"
    out.to_netcdf(out_path)
    print("wrote", out_path)
    return (out_path,)


@app.cell
def _(out_path, ux, uxds):
    import numpy as np

    uxds2 = ux.open_dataset(out_path)
    print("round-trip equal:", np.allclose(uxds2["psi"].values, uxds["psi"].values))
    return


@app.cell
def _(mo):
    mo.md(
        r"""
    ## The UXarray home turf: MPAS natively

    ✦ **Julia side** of this story: `load_ugrid` on the converted UGRID file
    raises `ArgumentError: UGRID dataset is missing variable Mesh2` —
    ManifoldFields expects the topology container under its own writer name
    `Mesh2`, while uxarray names it `grid_topology`; and even past that, an
    MPAS quadrilateral topology matches none of the four manifold-native grid
    types (see `DATA.md`). That boundary is exactly where UXarray's generic
    UGRID model wins.
    """
    )
    return


@app.cell
def _(DATA, plt, ux):
    uxds_mpas = ux.open_dataset(f"{DATA}/oQU480.grid.nc", f"{DATA}/oQU480.data.nc")
    name = list(uxds_mpas.data_vars)[0]
    print("first data var:", name)
    uxds_mpas[name].plot()
    plt.gca()
    return


@app.cell
def _(DATA, ux):
    uxds_hp = ux.UxDataset.from_healpix(f"{DATA}/psi_healpix.nc")
    print("healpix cells:", int(uxds_hp.uxgrid.n_face))
    return


@app.cell
def _(mo):
    mo.md(
        r"""
    ## Recap

    One entry point opens any topology; UGRID is the interchange format both
    tools agree on. Next: **03 · Fields & FieldSet**.
    """
    )
    return


if __name__ == "__main__":
    app.run()
