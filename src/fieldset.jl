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

import DimensionalData: basedims, combinedims, data_eltype

"""
    FieldSet(mesh, fields) <: AbstractDimStack

A multi-variable container for `DiscreteField`s that share one mesh. Fields
may live at different locations (`NodeLoc`, `EdgeLoc`, `CellLoc`).

`fields` may be a `NamedTuple`, a `Pair{Symbol,<:DiscreteField}` varargs, or a
`Dict{Symbol,<:DiscreteField}`. `FieldSet(:u => u, :v => v)` uses the first
member's mesh as the authority. Every member must reference the *same* mesh
object as the authority (`===`); use `withmesh` to rebind fields first.

Indexing with a `Symbol` returns the corresponding `DiscreteField`.
"""
struct FieldSet{
    K,
    T,
    N,
    L,
    D,
    R,
    LD,
    M,
    LM,
    MT <: AbstractManifoldMesh
} <: AbstractDimStack{K, T, N, L}
    data::L
    dims::D
    refdims::R
    layerdims::NamedTuple{K, LD}
    metadata::M
    layermetadata::NamedTuple{K, LM}
    mesh::MT
    function FieldSet(
            data, dims, refdims, layerdims::NamedTuple{K, LD}, metadata,
            layermetadata::NamedTuple{K}, mesh::MT
    ) where {K, LD, MT <: AbstractManifoldMesh}
        for key in K
            ld = Tuple(layerdims[key])
            Loc = _loc_from_layerdims(ld)
            locdim = only(d
            for d in ld if Symbol(DimensionalData.name(d)) in _LOCATION_DIM_NAMES)
            # layerdims entries are basedims (no index values), so resolve the
            # location dim against the combined stack dims before taking length
            locdim = DimensionalData.dims(dims, locdim)
            expected = expected_location_length(Loc, mesh)
            length(locdim) == expected || throw(DimensionMismatch(
                "field :$key location dimension length $(length(locdim)) != expected $expected"))
        end
        T = data_eltype(data)
        N = length(dims)
        return FieldSet{K, T, N}(
            data, dims, refdims, layerdims, metadata, layermetadata, mesh)
    end
    function FieldSet{K, T, N}(
            data::L, dims::D, refdims::R, layerdims::NamedTuple{K, LD}, metadata::M,
            layermetadata::NamedTuple{K, LM}, mesh::MT
    ) where {K, T, N, L, D, R, LD, M, LM, MT <: AbstractManifoldMesh}
        return new{K, T, N, L, D, R, LD, M, LM, MT}(
            data, dims, refdims, layerdims, metadata, layermetadata, mesh)
    end
end

function FieldSet(
        mesh0::M, fields::NamedTuple{K, <:Tuple{Vararg{<:DiscreteField}}}
) where {M <: AbstractManifoldMesh, K}
    isempty(fields) && throw(ArgumentError("FieldSet requires at least one field"))
    for key in K
        mesh(fields[key]) === mesh0 || throw(DimensionMismatch(
            "field :$key mesh does not match the FieldSet mesh (=== required); use withmesh to rebind"))
    end
    arrays = map(parent, fields)
    ds = combinedims(collect(Tuple(fields)))
    lds = map(basedims, fields)
    md = map(DimensionalData.metadata, fields)
    return FieldSet(arrays, ds, (), lds, DimensionalData.NoMetadata(), md, mesh0)
end

function FieldSet(mesh0::AbstractManifoldMesh, fields::Pair{Symbol, <:DiscreteField}...)
    FieldSet(mesh0, (; fields...))
end

function FieldSet(fields::Pair{Symbol, <:DiscreteField}...)
    nt = (; fields...)
    isempty(nt) && throw(ArgumentError("FieldSet requires at least one field"))
    return FieldSet(mesh(first(nt)), nt)
end

function FieldSet(mesh0::AbstractManifoldMesh, fields::Dict{Symbol, <:DiscreteField})
    FieldSet(mesh0, (; sort!(collect(fields); by = first)...))
end

mesh(fs::FieldSet) = getfield(fs, :mesh)

"""    field_names(fs::FieldSet)

Return the field name tuple `(Symbol, ...)` of a `FieldSet`."""
field_names(fs::FieldSet) = keys(fs)

import DimensionalData: layerdims, layermetadata, refdims

for f in (:getindex, :view, :dotview)
    @eval function Base.$f(fs::FieldSet, key::Symbol)
        # NamedTuple Symbol indexing throws ErrorException, not KeyError
        key in keys(fs) || throw(KeyError(key))
        ld = layerdims(fs, key)
        return DiscreteField(
            _loc_from_layerdims(Tuple(ld)), mesh(fs), DimensionalData.data(fs)[key], ld;
            name = key,
            metadata = layermetadata(fs, key),
            refdims = refdims(fs)
        )
    end
end

"""    fields(fs::FieldSet) -> NamedTuple

Materialize every member as a `DiscreteField`."""
function fields(fs::FieldSet{K}) where {K}
    return NamedTuple{K}(map(key -> fs[key], K))
end

function Base.show(io::IO, fs::FieldSet)
    print(io, "FieldSet on ", typeof(mesh(fs)), " with fields ", field_names(fs))
end
