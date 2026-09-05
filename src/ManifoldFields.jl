module ManifoldFields

using DimensionalData
using ManifoldMeshes
using StaticArrays
using YAXArrays

export DiscreteField, mesh, data, location, withmesh
export location_dimname, expected_location_length
export interpolate
export UGridDataset, to_ugrid, from_ugrid, from_ugrid_mesh, save_ugrid, load_ugrid,
       load_ugrid_mesh

include("locations.jl")
include("field.jl")
include("fieldset.jl")
include("broadcast.jl")
include("fieldset_broadcast.jl")
include("interpolate.jl")
include("ugrid.jl")

end
