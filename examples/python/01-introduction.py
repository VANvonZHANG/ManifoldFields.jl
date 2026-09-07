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
    # 01 · Introduction

    The same dataset, opened and explored in UXarray (Python) and ManifoldFields.jl
    (Julia). Twin notebook: `julia/01-introduction.jl`.
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
    print("data root:", DATA)
    return (DATA,)


@app.cell
def _():
    import uxarray as ux
    return (ux,)


@app.cell
def _():
    import matplotlib.pyplot as plt
    return (plt,)


@app.cell
def _(DATA, ux):
    uxds = ux.open_dataset(f"{DATA}/healpix_analytic.nc")
    uxds
    return (uxds,)


@app.cell
def _(mo):
    mo.md(
        r"""
    ✦ **Julia equivalent** (runs in the twin notebook):

    ```julia
    fs = load_ugrid(joinpath(DATA, "healpix_analytic.nc"))  # FieldSet
    f = fs[:psi]                                           # DiscreteField
    ManifoldFields.mesh(fs)                                # the HEALPixGrid
    ```
    """
    )
    return


@app.cell
def _(uxds):
    psi = uxds["psi"]
    print("shape:", psi.shape)
    print("range:", float(psi.min()), float(psi.max()))
    return (psi,)


@app.cell
def _(mo):
    mo.md(
        r"""
    `psi` is an analytic field, `psi(lat, lon) = cos(lat)·sin(lon) + 0.3·sin(2·lat)`,
    stored at cell centers of a HEALPix grid with nside = 32 (12 288 cells) — written
    by `save_ugrid` in Julia, read natively by UXarray.
    """
    )
    return


@app.cell
def _(plt, psi):
    psi.plot()
    plt.gca()
    return


@app.cell
def _(mo):
    mo.md(
        r"""
    ## Recap

    | here (UXarray) | twin (ManifoldFields) |
    |---|---|
    | `ux.open_dataset(path)` | `load_ugrid(path)` |
    | `UxDataset` | `FieldSet` |
    | `UxDataArray` | `DiscreteField` |
    | `uxds.uxgrid` | `mesh(fs)` |

    Next: **02 · Grids & UGRID IO**.
    """
    )
    return


if __name__ == "__main__":
    app.run()
