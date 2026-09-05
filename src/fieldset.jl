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
    # format wraps raw-range lookups (e.g. `Dim{:time}(1:3)`) into proper
    # `Sampled` lookups, as DimArray construction does; without this, inherited
    # reductions (`sum(fs; dims=...)`) fail in `reducelookup` on layer dims.
    # Formatted dims stay `==` to the raw input dims.
    ds = DimensionalData.format(combinedims(collect(Tuple(fields))))
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

import DimensionalData: rebuild, rebuild_from_arrays
import DimensionalData: AbstractBasicDimArray

function DimensionalData.rebuild(
        s::FieldSet;
        data = DimensionalData.data(s),
        dims = DimensionalData.dims(s),
        refdims = DimensionalData.refdims(s),
        layerdims = DimensionalData.layerdims(s),
        metadata = DimensionalData.metadata(s),
        layermetadata = DimensionalData.layermetadata(s)
)
    # never use `s.mesh` here: AbstractDimStack getproperty redirects to getindex
    return FieldSet(
        data, dims, refdims, layerdims, metadata, layermetadata, getfield(s, :mesh))
end
function DimensionalData.rebuild(
        s::FieldSet, data, dims = DimensionalData.dims(s),
        refdims = DimensionalData.refdims(s),
        layerdims = DimensionalData.layerdims(s),
        metadata = DimensionalData.metadata(s),
        layermetadata = DimensionalData.layermetadata(s)
)
    return FieldSet(
        data, dims, refdims, layerdims, metadata, layermetadata, getfield(s, :mesh))
end

function DimensionalData.rebuild_from_arrays(
        s::FieldSet, das::Tuple{Vararg{<:AbstractBasicDimArray}}; kw...)
    return rebuild_from_arrays(s, NamedTuple{field_names(s)}(das); kw...)
end
function DimensionalData.rebuild_from_arrays(
        s::FieldSet,
        das::NamedTuple{<:Any, <:Tuple{Vararg{<:AbstractBasicDimArray}}};
        data = map(parent, das),
        refdims = DimensionalData.refdims(s),
        metadata = DimensionalData.metadata(s),
        dims = nothing,
        layerdims = map(basedims, das),
        layermetadata = map(DimensionalData.metadata, das)
)
    if isnothing(dims)
        # layerdims are basedims (Colon lookups); the inner constructor resolves
        # location-dim lengths against these stack-level combinedims
        dims = combinedims(collect(Tuple(das)))
    end
    return FieldSet(
        data, dims, refdims, layerdims, metadata, layermetadata, getfield(s, :mesh))
end

"""
    getindex(fs::FieldSet, d::Dimension, ds::Dimension...)

Slice every field of the `FieldSet` along the given dimension selector(s),
returning a new `FieldSet` on the same mesh; fields that contain none of the
selected dimensions are passed through unchanged (a warning is emitted when
no field contains any of them). Slicing a location dimension (`node`, `edge`,
or `cell`) throws an `ArgumentError` because fields keep their mesh-aligned
location extent.
"""
function Base.getindex(
        fs::FieldSet, d::DimensionalData.Dimension, ds::DimensionalData.Dimension...; kw...)
    sel = (d, ds...)
    selnames = Set(Symbol(DimensionalData.name(s)) for s in sel)
    isempty(intersect(selnames, Set(_LOCATION_DIM_NAMES))) || throw(ArgumentError(
        "cannot slice a location dimension ($(join(_LOCATION_DIM_NAMES, "/"))); " *
        "FieldSet layers keep their mesh-aligned location extent"))
    all_layer_dims = union!(
        Set{Symbol}(),
        (Set(Symbol(DimensionalData.name(x)) for x in DimensionalData.dims(fs[key]))
        for key in field_names(fs))...)
    if isempty(intersect(selnames, all_layer_dims))
        @warn "FieldSet getindex: dimension(s) $(sort(collect(selnames))) not found in any field; returning an unchanged copy"
    end

    layers = map(field_names(fs)) do key
        f = fs[key]
        flagnames = Set(Symbol(DimensionalData.name(x)) for x in DimensionalData.dims(f))
        if isempty(intersect(selnames, flagnames))
            return f   # layer untouched by this selection
        end
        return DimensionalData.getindex(f, d, ds...; kw...)
    end
    return FieldSet(mesh(fs), NamedTuple{field_names(fs)}(Tuple(layers)))
end

function Base.merge(s::FieldSet, pairs::Pair{Symbol, <:DiscreteField}...)
    for (key, f) in pairs
        mesh(f) === mesh(s) || throw(DimensionMismatch(
            "field :$key mesh does not match the FieldSet mesh (=== required); use withmesh to rebind"))
    end
    merged = merge(fields(s), (; pairs...))
    return FieldSet(mesh(s), merged)
end
function Base.merge(s::FieldSet, nt::NamedTuple{<:Any, <:Tuple{Vararg{<:DiscreteField}}})
    for key in keys(nt)
        mesh(nt[key]) === mesh(s) || throw(DimensionMismatch(
            "field :$key mesh does not match the FieldSet mesh (=== required); use withmesh to rebind"))
    end
    return FieldSet(mesh(s), merge(fields(s), nt))
end
Base.merge(s::FieldSet) = s

function Base.:(==)(s1::FieldSet, s2::FieldSet)
    mesh(s1) === mesh(s2) || return false
    field_names(s1) == field_names(s2) || return false
    return DimensionalData.data(s1) == DimensionalData.data(s2) &&
           DimensionalData.layerdims(s1) == DimensionalData.layerdims(s2)
end

# Reductions. DimensionalData's inherited stack reductions keep each
# fully-reduced dimension as a length-1 dimension in every affected layer; for
# FieldSet those dimensions are dropped instead, so `sum(fs; dims = Dim{:time})`
# yields fields whose dims no longer contain the reduced dimension. Reducing
# along a location dimension throws: the per-layer rebuild goes through
# DiscreteField validation, which rejects a mesh location extent of length 1
# (layers that have the location dim) or a layer left with no location
# dimension at all.
import Statistics

for (mod,
    fnames) in (:Base => (:sum, :prod, :maximum, :minimum, :extrema),
    :Statistics => (:mean, :median, :std, :var))
    for fname in fnames
        @eval function $(mod).$(fname)(s::FieldSet; dims = :, kw...)
            return DimensionalData.maplayers(s) do A
                if dims isa Colon
                    return $(mod).$(fname)(A; dims = :, kw...)
                end
                ld = DimensionalData.commondims(A, dims)
                isempty(ld) && return A   # layer untouched by this reduction
                return Base.dropdims($(mod).$(fname)(A; dims = ld, kw...); dims = ld)
            end
        end
    end
end
for (mod,
    fnames) in (:Base => (:reduce, :sum, :prod, :maximum, :minimum, :extrema),
    :Statistics => (:mean,))
    for fname in fnames
        @eval function $(mod).$(fname)(f::Function, s::FieldSet; dims = Colon())
            return DimensionalData.maplayers(s) do A
                if dims isa Colon
                    return $(mod).$(fname)(f, A; dims = :)
                end
                ld = DimensionalData.commondims(A, dims)
                isempty(ld) && return A   # layer untouched by this reduction
                return Base.dropdims($(mod).$(fname)(f, A; dims = ld); dims = ld)
            end
        end
    end
end
