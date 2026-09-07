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
    # 04 · Point Interpolation

    Ask for the value *between* grid points. `interpolate` gives true on-manifold
    interpolation at any (lat, lon) — today for `NodeLoc` fields, with batched
    queries on the roadmap. Twin notebook: `julia/04-interpolation.jl`.
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
    return (DATA,)


@app.cell
def _():
    import numpy as np
    import uxarray as ux
    return np, ux


@app.cell
def _(DATA, ux):
    uxds = ux.open_dataset(f"{DATA}/healpix_analytic.nc")
    return (uxds,)


@app.cell
def _(mo):
    mo.md(
        r"""
    ✦ **Julia equivalent** — one call, and the result is *interpolated*, not
    nearest-cell:

    ```julia
    f = DiscreteField(NodeLoc, g, node_values, dims; name = :psi)
    interpolate(f, 12.3, 45.6)   # ≈ analytic(12.3, 45.6)
    ```
    """
    )
    return


@app.cell
def _(np):
    def analytic(lat, lon):
        return np.cos(np.deg2rad(lat)) * np.sin(np.deg2rad(lon)) + 0.3 * np.sin(
            np.deg2rad(2 * lat)
        )

    def to_xyz(lat, lon):
        la, lo = np.deg2rad(lat), np.deg2rad(lon)
        return [np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo), np.sin(la)]

    return analytic, to_xyz


@app.cell
def _(uxds):
    bt = uxds.uxgrid.get_ball_tree(
        coordinate_system="cartesian", distance_metric="euclidean"
    )
    return (bt,)


@app.cell
def _(bt, to_xyz, uxds):
    dist, idx = bt.query([to_xyz(12.3, 45.6)])
    nearest = float(uxds["psi"].values[idx.item()])
    print("nearest cell value:", nearest)
    return (nearest,)


@app.cell
def _(analytic, nearest):
    print("analytic truth:     ", analytic(12.3, 45.6))
    return


@app.cell
def _(mo):
    mo.md(
        r"""
    ## Error sweep

    Nearest-cell against the analytic truth on 200 random points (HEALPix
    nside = 32 ⇒ ~1.8° cells). The Julia twin runs a comparable sweep with true
    interpolation — currently `0.105`, dominated by the LatLonGrid weight bug
    documented in its final cell, versus `0.026` here.
    """
    )
    return


@app.cell
def _(analytic, bt, np, to_xyz, uxds):
    rng = np.random.default_rng(0)
    lats = rng.uniform(-80, 80, 200)
    lons = rng.uniform(0, 360, 200)
    _, ind = bt.query(np.array([to_xyz(a, b) for a, b in zip(lats, lons)]))
    vals = uxds["psi"].values[ind].ravel()
    truth = np.array([analytic(a, b) for a, b in zip(lats, lons)])
    print("nearest-cell max error:", float(np.abs(vals - truth).max()))
    return


@app.cell
def _(mo):
    mo.md(
        r"""
    ## Recap

    For grid-to-grid work UXarray offers `remap.nearest_neighbor` /
    `remap.bilinear` / `remap.inverse_distance_weighted`; for point queries the
    ball tree is the tool. ManifoldFields' `interpolate` answers point queries
    with interpolated values directly. Next: `CAPABILITY_MAP.md`.
    """
    )
    return


if __name__ == "__main__":
    app.run()
