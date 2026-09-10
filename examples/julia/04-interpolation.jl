### A Pluto.jl notebook ###
# ╔═╡ 00000000-0000-4000-8000-000000000031
md"""
# 04 · Point Interpolation

Ask for the value *between* grid points. `interpolate` gives true on-manifold
interpolation at any (lat, lon) — today for `NodeLoc` fields, with batched
queries on the roadmap. Twin notebook: `python/04-interpolation.py`.
"""

# ╔═╡ 00000000-0000-4000-8000-00000000003b
begin
    import Pkg
    Pkg.activate(@__DIR__; io = devnull)
    if isnothing(Base.find_package("ManifoldFields")) ||
       isnothing(Base.find_package("ManifoldMeshes"))
        # Bootstrap on a fresh clone (Pluto's nbpkg cannot use [sources] on Julia 1.10):
        Pkg.develop(path = joinpath(@__DIR__, "..", ".."); io = devnull)
        Pkg.add(url = "https://github.com/VANvonZHANG/ManifoldMeshes.jl"; io = devnull)
    end
end

# ╔═╡ 00000000-0000-4000-8000-000000000032
begin
    using ManifoldFields
    using ManifoldMeshes
    using DimensionalData
end

# ╔═╡ 00000000-0000-4000-8000-000000000033
begin
    analytic(lat, lon) = cosd(lat) * sind(lon) + 0.3 * sind(2lat)

    function lonlat(p)
        lat = asind(p[3])
        lon = atand(p[2], p[1])
        lon < 0 && (lon += 360.0)
        return (lat, lon)
    end
end

# ╔═╡ 00000000-0000-4000-8000-000000000034
md"""
Sample the analytic field at every **node** of a 5° lat-lon grid (interpolation
is defined for `NodeLoc` fields), then forget the formula and ask the field.
"""

# ╔═╡ 00000000-0000-4000-8000-000000000035
begin
    g = LatLonGrid(lat_edges = collect(-90.0:5.0:90.0), lon_edges = collect(0.0:5.0:360.0))
    vals = [analytic(lonlat(node_coordinates(g, i))...) for i in 1:num_nodes(g)]
    f = DiscreteField(NodeLoc, g, vals, (Dim{:node}(1:num_nodes(g)),); name = :psi)
end

# ╔═╡ 00000000-0000-4000-8000-000000000036
interpolate(f, 12.3, 45.6)

# ╔═╡ 00000000-0000-4000-8000-000000000037
analytic(12.3, 45.6)

# ╔═╡ 00000000-0000-4000-8000-000000000038
md"""
✦ **Python equivalent** — UXarray's point query is *nearest cell* via a ball
tree (true interpolation between points is offline remapping, `remap.bilinear`):

```python
bt = uxds.uxgrid.get_ball_tree(coordinate_system="cartesian",
                               distance_metric="euclidean")
dist, idx = bt.query([to_xyz(lat, lon)])
nearest = uxds["psi"].values[idx.item()]
```
"""

# ╔═╡ 00000000-0000-4000-8000-000000000039
begin
    errs = [abs(interpolate(f, lat, lon) - analytic(lat, lon))
            for lat in -80.0:8.0:80.0, lon in 10.0:17.0:350.0]
    maximum(errs)
end

# ╔═╡ 00000000-0000-4000-8000-00000000003a
md"""
Maximum interpolation error over the 441-point sweep against the analytic
truth — measured **0.0025** on ManifoldMeshes v0.7.1. (On v0.6.0/v0.7.0 the
same sweep measured **0.105**, which is how it earned its keep: it exposed a
real bug, not a too-coarse grid. `LatLonGrid` fed `_bilinear_weights` the
latitude and longitude fractions swapped, so the SE and NW corner weights
landed on each other's nodes — at (-80, 299) the weights came back
`(0.2, 0.0, 0.0, 0.8)` where correct bilinear gives `(0.2, 0.8, 0.0, 0.0)`.
The error `|latfrac − lonfrac| · |V_SE − V_NW|` vanished at nodes, along
cell diagonals, and at centroids, hiding it from the test suite.) With the
transposition fixed upstream (`s` is the SW→SE longitude fraction, `t` the
SW→NW latitude fraction), true bilinear on a 5° grid now sits an order of
magnitude below UXarray's nearest-cell `0.026` — as it should.

## Recap

| here (ManifoldFields) | twin (UXarray) |
|---|---|
| `interpolate(f, lat, lon)` — interpolated value | ball-tree nearest cell (`get_ball_tree()`) |
| point queries today; batched queries roadmap | grid-to-grid via `remap.nearest_neighbor` / `.bilinear` |

The full UXarray ↔ ManifoldFields picture: [`CAPABILITY_MAP.md`](../CAPABILITY_MAP.md).
"""
