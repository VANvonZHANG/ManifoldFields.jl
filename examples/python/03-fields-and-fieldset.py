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
    # 03 · Fields & FieldSet

    Values live at a *location* of the mesh — nodes, cells, or edges — and that
    fact is part of the type. A `FieldSet` holds several fields on one shared
    mesh, like a `UxDataset`. Twin notebook: `julia/03-fields-and-fieldset.jl`.
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
    import matplotlib.pyplot as plt
    import uxarray as ux
    return plt, ux


@app.cell
def _(DATA, ux):
    uxds = ux.open_dataset(f"{DATA}/healpix_analytic.nc")
    psi = uxds["psi"]
    print("dims:", psi.dims)
    psi
    return (psi, uxds)


@app.cell
def _(mo):
    mo.md(
        r"""
    Location is implied by the dimension name — face data carries
    `("n_face",)`; node data would carry `("n_node",)`.

    ✦ **Julia equivalent**: `fs_h[:psi]` is a `DiscreteField{CellLoc}` —
    location is a *type parameter* there, not a naming convention.
    """
    )
    return


@app.cell
def _(psi):
    anomaly = psi - psi.mean()
    print("anomaly range:", float(anomaly.min()), float(anomaly.max()))
    return (anomaly,)


@app.cell
def _(mo):
    mo.md(
        r"""
    ## Reductions

    xarray reductions **drop** the reduced dimension by default;
    `keepdims=True` retains it as size 1 — ManifoldFields made the same
    choice (`keepdims = true`), so this default is one the twins share
    (numpy is the odd one out). Note *which* dimension is being reduced
    below: `n_face` is the location dimension, and UXarray reduces it like
    any other. The Julia twin cannot — `mean(fs; dims = Dim{:node})` throws,
    because the mesh extent is part of the field's contract:

    ```julia
    mean(fs)                                        # NamedTuple of scalars
    mean(fs_t; dims = Dim{:time})                   # Dim{:time} dropped
    mean(fs_t; dims = Dim{:time}, keepdims = true)  # retained as length 1
    ```
    """
    )
    return


@app.cell
def _(psi):
    print("dropped by default:", psi.mean(dim=psi.dims[0]).sizes)
    print("kept with keepdims=True:", psi.mean(dim=psi.dims[0], keepdims=True).sizes)
    return


@app.cell
def _(anomaly, plt):
    anomaly.plot()
    plt.gca()
    return


@app.cell
def _(mo):
    mo.md(
        r"""
    ## Recap

    Same operations, different metadata contracts: both drop reduced dims by
    default, but UXarray keys location off dimension names while ManifoldFields
    types it — and refuses to reduce it. Next: **04 · Point Interpolation**.
    """
    )
    return


if __name__ == "__main__":
    app.run()
