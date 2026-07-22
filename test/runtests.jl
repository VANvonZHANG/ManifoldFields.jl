using Aqua
using ManifoldFields
using Test

@testset "Code quality (Aqua.jl)" begin
    Aqua.test_all(ManifoldFields)
end

import ManifoldMeshes
import ManifoldMeshes: LatLonGrid

function LatLonGrid(nlat::Integer, nlon::Integer)
    return ManifoldMeshes.LatLonGrid(
        lat_edges=collect(range(-90.0, 90.0; length=nlat + 1)),
        lon_edges=collect(range(0.0, 360.0; length=nlon + 1)),
    )
end

@testset "ManifoldFields.jl" begin
    include("test_locations.jl")
    include("test_field_construction.jl")
end
