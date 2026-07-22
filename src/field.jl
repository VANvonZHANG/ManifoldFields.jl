import DimensionalData
import DimensionalData: AbstractDimArray, dims
import ManifoldMeshes: AbstractLocation, AbstractManifoldMesh

struct DiscreteField{
    Loc<:AbstractLocation,
    T,N,D<:Tuple,A<:AbstractArray{T,N},
    M<:AbstractManifoldMesh,
    R,
    NM,
    MD,
} <: AbstractDimArray{T,N,D,A}
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
        mesh::M,
    ) where {Loc<:AbstractLocation,D<:Tuple,A<:AbstractArray{T,N},
             M<:AbstractManifoldMesh,R,NM,MD} where {T,N}
        return new{Loc,T,N,D,A,M,R,NM,MD}(data, dims, refdims, name, metadata, mesh)
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
Base.IndexStyle(::Type{<:DiscreteField{<:Any,<:Any,<:Any,<:Any,A}}) where {A} =
    Base.IndexStyle(A)
Base.parent(f::DiscreteField) = f.data

function _dim_name(d)
    return DimensionalData.name(d)
end

function _location_axes(Loc, ds)
    locname = DimensionalData.name(location_dimname(Loc))
    return findall(d -> _dim_name(d) == locname, collect(ds))
end

function _validate_location_dims(::Type{Loc}, mesh, values, ds) where {Loc<:AbstractLocation}
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
    name=:field,
    metadata=nothing,
    refdims=(),
) where {Loc<:AbstractLocation,M<:AbstractManifoldMesh,A<:AbstractArray}
    checked_dims = _validate_location_dims(Loc, mesh, values, ds)
    return DiscreteField{Loc}(
        values, checked_dims, refdims, name, metadata, mesh
    )
end

function DiscreteField(
    ::Type{Loc},
    mesh::M,
    data::AbstractDimArray;
    name=DimensionalData.name(data),
    metadata=DimensionalData.metadata(data),
    refdims=DimensionalData.refdims(data),
) where {Loc<:AbstractLocation,M<:AbstractManifoldMesh}
    return DiscreteField(Loc, mesh, DimensionalData.data(data), DimensionalData.dims(data);
                         name=name, metadata=metadata, refdims=refdims)
end

function DiscreteField(::Type{Loc}, mesh::M, values::AbstractArray; kwargs...) where {
    Loc<:AbstractLocation,
    M<:AbstractManifoldMesh,
}
    throw(ArgumentError("plain AbstractArray input requires explicit dims"))
end

function withmesh(f::DiscreteField{Loc}, mesh::AbstractManifoldMesh) where {Loc}
    return DiscreteField(Loc, mesh, data(f), dims(f); name=DimensionalData.name(f),
                         metadata=DimensionalData.metadata(f),
                         refdims=DimensionalData.refdims(f))
end

function DimensionalData.rebuild(
    f::DiscreteField{Loc},
    values::AbstractArray,
    ds::Tuple=DimensionalData.dims(f),
    refdims=DimensionalData.refdims(f),
    name=DimensionalData.name(f),
    metadata=DimensionalData.metadata(f),
) where {Loc}
    return DiscreteField(Loc, mesh(f), values, ds; name=name, metadata=metadata,
                         refdims=refdims)
end

function DimensionalData.rebuild(
    f::DiscreteField{Loc};
    data=parent(f),
    dims=DimensionalData.dims(f),
    refdims=DimensionalData.refdims(f),
    name=DimensionalData.name(f),
    metadata=DimensionalData.metadata(f),
) where {Loc}
    return DiscreteField(Loc, mesh(f), data, dims; name=name, metadata=metadata,
                         refdims=refdims)
end
