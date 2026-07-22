import Base.Broadcast: AbstractArrayStyle, BroadcastStyle, Broadcasted

struct DiscreteFieldStyle{N} <: AbstractArrayStyle{N} end

DiscreteFieldStyle{N}(::Val{N}) where {N} = DiscreteFieldStyle{N}()

BroadcastStyle(::Type{<:DiscreteField{<:Any,<:Any,N}}) where {N} = DiscreteFieldStyle{N}()

function Base.copy(bc::Broadcasted{DiscreteFieldStyle{N}}) where {N}
    fields = Any[]
    _collect_fields!(fields, bc)
    isempty(fields) && return copy(_unwrap_fields(bc))

    anchor = first(fields)
    _validate_broadcast_fields(anchor, fields)

    values = copy(_unwrap_fields(bc))
    return DiscreteField(
        location(anchor),
        mesh(anchor),
        values,
        dims(anchor);
        name=DimensionalData.name(anchor),
        metadata=DimensionalData.metadata(anchor),
        refdims=DimensionalData.refdims(anchor),
    )
end

function _validate_broadcast_fields(anchor::DiscreteField, fields)
    for f in fields
        mesh(f) === mesh(anchor) ||
            throw(DimensionMismatch("broadcast requires identical mesh objects"))
        location(f) === location(anchor) ||
            throw(DimensionMismatch("broadcast requires identical field locations"))
        dims(f) == dims(anchor) ||
            throw(DimensionMismatch("broadcast requires identical field dimensions"))
    end
    return nothing
end

function _collect_fields!(fields, bc::Broadcasted)
    foreach(arg -> _collect_fields!(fields, arg), bc.args)
    return fields
end
function _collect_fields!(fields, ex::Base.Broadcast.Extruded)
    _collect_fields!(fields, ex.x)
    return fields
end
function _collect_fields!(fields, f::DiscreteField)
    push!(fields, f)
    return fields
end
_collect_fields!(fields, _) = fields

function _unwrap_fields(bc::Broadcasted)
    return Broadcasted(bc.f, map(_unwrap_fields, bc.args))
end
function _unwrap_fields(ex::Base.Broadcast.Extruded)
    return Base.Broadcast.Extruded(_unwrap_fields(ex.x), ex.keeps, ex.defaults)
end
_unwrap_fields(f::DiscreteField) = parent(f)
_unwrap_fields(x) = x
