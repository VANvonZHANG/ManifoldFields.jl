import DimensionalData
import DimensionalData: AbstractDimStack
import ManifoldMeshes: AbstractLocation, AbstractManifoldMesh

const _LOCATION_DIM_NAMES = (:node, :edge, :cell)

"""Return the location type implied by a layer's dimension names.

Throws `ArgumentError` unless exactly one dimension is named `:node`, `:edge`,
or `:cell`.
"""
function _loc_from_layerdims(ld::Tuple)
    hit = Symbol[]
    for d in ld
        n = Symbol(DimensionalData.name(d))
        n in _LOCATION_DIM_NAMES && push!(hit, n)
    end
    length(hit) == 1 || throw(ArgumentError(
        "expected exactly one location dimension ($(join(_LOCATION_DIM_NAMES, "/"))), found $(length(hit))"))
    n = only(hit)
    return n === :node ? NodeLoc : (n === :edge ? EdgeLoc : CellLoc)
end
