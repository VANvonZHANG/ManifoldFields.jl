# Authoring script: write analytic test fields on three ManifoldMeshes grids to
# examples/data/*.nc (UGRID via ManifoldFields.save_ugrid), so that both the
# Julia and Python notebooks can check values against ground truth.
#
# Usage: julia --project=<repo root> examples/scripts/make_julia_data.jl
using ManifoldFields
using ManifoldMeshes
using DimensionalData

const HERE = @__DIR__
const DATA = joinpath(dirname(HERE), "data")
mkpath(DATA)

# Ground truth: smooth in lat/lon, cheap to evaluate in both languages.
analytic(lat, lon) = cosd(lat) * sind(lon) + 0.3 * sind(2lat)

function lonlat(p)
    lat = asind(p[3])
    lon = atand(p[2], p[1])
    lon < 0 && (lon += 360.0)
    return (lat, lon)
end

meshes = (
    ("latlon_analytic",
        LatLonGrid(lat_edges = collect(-90.0:5.0:90.0),
            lon_edges = collect(0.0:5.0:360.0))),
    ("cubedsphere_analytic", CubedSphereGrid(n = 24)),
    ("healpix_analytic", HEALPixGrid(nside = 32))
)

for (name, m) in meshes
    vals = [analytic(lonlat(cell_centroid(m, i))...) for i in 1:num_cells(m)]
    field = DiscreteField(CellLoc, m, vals, (Dim{:cell}(1:num_cells(m)),); name = :psi)
    fs = FieldSet(m, :psi => field)
    path = joinpath(DATA, "$name.nc")
    save_ugrid(fs, path)
    println("wrote $path  (cells = $(num_cells(m)), psi range = $(extrema(vals)))")
end
println("done")
