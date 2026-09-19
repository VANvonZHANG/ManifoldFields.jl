module ManifoldFields

using DimensionalData
using ManifoldMeshes
using StaticArrays
using YAXArrays

export DiscreteField, mesh, data, location, withmesh
export FieldSet, fields, field_names
export location_dimname, location_axis, expected_location_length
export interpolate
export UGridDataset, to_ugrid, from_ugrid, from_ugrid_mesh, save_ugrid, load_ugrid,
       load_ugrid_mesh

# `field.jl` defines `DiscreteField`, which `locations.jl`'s public
# `location_axis(f::DiscreteField{Loc})` names in its signature; signatures are
# resolved when the method is defined, so the type must come first.
include("field.jl")
include("locations.jl")
include("fieldset.jl")
include("broadcast.jl")
include("fieldset_broadcast.jl")
include("interpolate.jl")
include("ugrid.jl")

end
