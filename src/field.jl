import DimensionalData
import DimensionalData: AbstractDimArray, dims
import ManifoldMeshes: AbstractLocation, AbstractManifoldMesh

struct DiscreteField{
    Loc <: AbstractLocation,
    T, N, D <: Tuple, A <: AbstractArray{T, N},
    M <: AbstractManifoldMesh,
    R,
    NM,
    MD
} <: AbstractDimArray{T, N, D, A}
    data::A
    dims::D
    refdims::R
    name::NM
    metadata::MD
    mesh::M
    function DiscreteField{Loc}(
            data::A,
            dims::D,
            refdims::R,
            name::NM,
            metadata::MD,
            mesh::M
    ) where {Loc <: AbstractLocation, D <: Tuple, A <: AbstractArray{T, N},
            M <: AbstractManifoldMesh, R, NM, MD} where {T, N}
        return new{Loc, T, N, D, A, M, R, NM, MD}(data, dims, refdims, name, metadata, mesh)
    end
end

mesh(f::DiscreteField) = f.mesh
data(f::DiscreteField) = f.data
location(::DiscreteField{Loc}) where {Loc} = Loc

DimensionalData.dims(f::DiscreteField) = f.dims
DimensionalData.refdims(f::DiscreteField) = f.refdims
DimensionalData.name(f::DiscreteField) = f.name
DimensionalData.metadata(f::DiscreteField) = f.metadata
DimensionalData.data(f::DiscreteField) = f.data
Base.size(f::DiscreteField) = size(f.data)
Base.size(f::DiscreteField, i::Integer) = size(f.data, i)
function Base.IndexStyle(::Type{<:DiscreteField{<:Any, <:Any, <:Any, <:Any, A}}) where {A}
    Base.IndexStyle(A)
end
Base.parent(f::DiscreteField) = f.data

function _valid_field_name(name)
    return isempty(String(Symbol(name))) ? :field : name
end

function _dim_name(d)
    return DimensionalData.name(d)
end

function _location_axes(Loc, ds)
    locname = DimensionalData.name(location_dimname(Loc))
    return findall(d -> _dim_name(d) == locname, collect(ds))
end

function _validate_location_dims(::Type{Loc}, mesh, values, ds) where {Loc <:
                                                                       AbstractLocation}
    axes = _location_axes(Loc, ds)
    length(axes) == 1 ||
        throw(ArgumentError("expected exactly one $(location_dimname(Loc)) dimension, found $(length(axes))"))
    length(ds) == ndims(values) ||
        throw(DimensionMismatch("dims length $(length(ds)) does not match data ndims $(ndims(values))"))
    for axis in 1:ndims(values)
        actual_axis = size(values, axis)
        dim_axis = length(ds[axis])
        actual_axis == dim_axis ||
            throw(DimensionMismatch("dimension $axis length $dim_axis != data axis length $actual_axis"))
    end
    axis = only(axes)
    expected = expected_location_length(Loc, mesh)
    actual = size(values, axis)
    actual == expected ||
        throw(DimensionMismatch("$(location_dimname(Loc)) length $actual != expected $expected"))
    return ds
end

function DiscreteField(
        ::Type{Loc},
        mesh::M,
        values::A,
        ds;
        name = :field,
        metadata = nothing,
        refdims = ()
) where {Loc <: AbstractLocation, M <: AbstractManifoldMesh, A <: AbstractArray}
    checked_dims = _validate_location_dims(Loc, mesh, values, ds)
    return DiscreteField{Loc}(
        values, checked_dims, refdims, _valid_field_name(name), metadata, mesh
    )
end

function DiscreteField(
        ::Type{Loc},
        mesh::M,
        data::AbstractDimArray;
        name = DimensionalData.name(data),
        metadata = DimensionalData.metadata(data),
        refdims = DimensionalData.refdims(data)
) where {Loc <: AbstractLocation, M <: AbstractManifoldMesh}
    return DiscreteField(Loc, mesh, DimensionalData.data(data), DimensionalData.dims(data);
        name = name, metadata = metadata, refdims = refdims)
end

function _dimarray(f::DiscreteField)
    return DimensionalData.DimArray(
        parent(f),
        DimensionalData.dims(f);
        refdims = DimensionalData.refdims(f),
        name = DimensionalData.name(f),
        metadata = DimensionalData.metadata(f)
    )
end

function _maybe_discrete_field(f::DiscreteField{Loc}, result) where {Loc}
    result isa AbstractDimArray || return result
    isempty(_location_axes(Loc, DimensionalData.dims(result))) && return result
    _validate_location_dims(Loc, mesh(f), parent(result), DimensionalData.dims(result))
    return DiscreteField(Loc, mesh(f), parent(result), DimensionalData.dims(result);
        refdims = DimensionalData.refdims(result),
        name = DimensionalData.name(result),
        metadata = DimensionalData.metadata(result))
end

function Base.getindex(f::DiscreteField, d::DimensionalData.Dimension,
        ds::DimensionalData.Dimension...; kw...)
    result = getindex(_dimarray(f), d, ds...; kw...)
    return _maybe_discrete_field(f, result)
end

function Base.getindex(f::DiscreteField; kw...)
    result = getindex(_dimarray(f); kw...)
    return _maybe_discrete_field(f, result)
end

function DiscreteField(::Type{Loc}, mesh::M, values::AbstractArray;
        kwargs...) where {
        Loc <: AbstractLocation,
        M <: AbstractManifoldMesh
}
    throw(ArgumentError("plain AbstractArray input requires explicit dims"))
end

function withmesh(f::DiscreteField{Loc}, mesh::AbstractManifoldMesh) where {Loc}
    return DiscreteField(Loc, mesh, data(f), dims(f); name = DimensionalData.name(f),
        metadata = DimensionalData.metadata(f),
        refdims = DimensionalData.refdims(f))
end

function DimensionalData.rebuild(
        f::DiscreteField{Loc},
        values::AbstractArray,
        ds::Tuple = DimensionalData.dims(f),
        refdims = DimensionalData.refdims(f),
        name = DimensionalData.name(f),
        metadata = DimensionalData.metadata(f)
) where {Loc}
    return DiscreteField(Loc, mesh(f), values, ds; name = name, metadata = metadata,
        refdims = refdims)
end

function DimensionalData.rebuild(
        f::DiscreteField{Loc};
        data = parent(f),
        dims = DimensionalData.dims(f),
        refdims = DimensionalData.refdims(f),
        name = DimensionalData.name(f),
        metadata = DimensionalData.metadata(f)
) where {Loc}
    return DiscreteField(Loc, mesh(f), data, dims; name = name, metadata = metadata,
        refdims = refdims)
end
