import Base.Broadcast: AbstractArrayStyle, BroadcastStyle, Broadcasted, DefaultArrayStyle

struct FieldSetStyle{N} <: AbstractArrayStyle{N} end

FieldSetStyle{N}(::Val{N}) where {N} = FieldSetStyle{N}()
function FieldSetStyle{N}(::Val{M}) where {N, M}
    throw(DimensionMismatch(
        "FieldSet broadcast must preserve rank/dims; set rank $N cannot produce rank $M",
    ))
end

BroadcastStyle(::Type{<:FieldSet{<:Any, <:Any, N}}) where {N} = FieldSetStyle{N}()
# FieldSet is iterable but not an AbstractArray, so Base's generic
# `broadcastable(x) = collect(x)` fallback would materialize it into an
# Array of NamedTuples before style inference. Return the stack itself
# (axes/ndims/getindex are defined by AbstractDimStack) so FieldSetStyle
# drives the broadcast.
Base.Broadcast.broadcastable(fs::FieldSet) = fs
BroadcastStyle(::FieldSetStyle{N}, ::DefaultArrayStyle{N}) where {N} = FieldSetStyle{N}()
BroadcastStyle(::DefaultArrayStyle{N}, ::FieldSetStyle{N}) where {N} = FieldSetStyle{N}()
BroadcastStyle(::FieldSetStyle{N}, ::FieldSetStyle{N}) where {N} = FieldSetStyle{N}()
function BroadcastStyle(::FieldSetStyle{N}, ::DimensionalData.DimensionalStyle) where {N}
    return FieldSetStyle{N}()
end
function BroadcastStyle(::DimensionalData.DimensionalStyle, ::FieldSetStyle{N}) where {N}
    return FieldSetStyle{N}()
end
function BroadcastStyle(::FieldSetStyle{N}, ::DiscreteFieldStyle{M}) where {N, M}
    return FieldSetStyle{N}()
end
function BroadcastStyle(::DiscreteFieldStyle{M}, ::FieldSetStyle{N}) where {N, M}
    return FieldSetStyle{N}()
end

function Base.copy(bc::Broadcasted{FieldSetStyle{N}}) where {N}
    fss = Any[]
    _collect_fieldsets!(fss, bc)
    anchor = first(fss)
    _validate_broadcast_fieldsets(anchor, fss)
    out = map(field_names(anchor)) do key
        copy(_replace_fieldsets(bc, key))
    end
    return FieldSet(mesh(anchor), NamedTuple{field_names(anchor)}(Tuple(out)))
end

function _validate_broadcast_fieldsets(anchor::FieldSet, fss)
    akeys = field_names(anchor)
    for fs in fss
        field_names(fs) == akeys || throw(DimensionMismatch(
            "FieldSet broadcast requires identical field names; got $(field_names(fs)) and $akeys"))
        mesh(fs) === mesh(anchor) ||
            throw(DimensionMismatch("FieldSet broadcast requires identical mesh objects"))
    end
    return nothing
end

function _collect_fieldsets!(fss, bc::Broadcasted)
    foreach(arg -> _collect_fieldsets!(fss, arg), bc.args)
    return fss
end
_collect_fieldsets!(fss, ex::Base.Broadcast.Extruded) = _collect_fieldsets!(fss, ex.x)
function _collect_fieldsets!(fss, fs::FieldSet)
    push!(fss, fs)
    return fss
end
_collect_fieldsets!(fss, _) = fss

# Pass `nothing` as axes (not the stack-level `bc.axes`, which is the union of
# all layer dims): the untyped `Broadcasted` constructor re-infers the style
# from the replaced args, and `instantiate` recomputes per-field axes that
# match the layer's own dimensionality.
function _replace_fieldsets(bc::Broadcasted, key::Symbol)
    Broadcasted(bc.f, map(a -> _replace_fieldsets(a, key), bc.args), nothing)
end
function _replace_fieldsets(ex::Base.Broadcast.Extruded, key::Symbol)
    return Base.Broadcast.Extruded(_replace_fieldsets(ex.x, key), ex.keeps, ex.defaults)
end
_replace_fieldsets(fs::FieldSet, key::Symbol) = fs[key]
_replace_fieldsets(x, key::Symbol) = x
