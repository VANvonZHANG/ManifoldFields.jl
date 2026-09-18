import ManifoldMeshes: NodeLoc, interpolation_weights, locate_cell
import StaticArrays: SVector

function _location_axis(::Type{Loc}, f::DiscreteField) where {Loc}
    locname = DimensionalData.name(location_dimname(Loc))
    matches = findall(d -> DimensionalData.name(d) == locname, collect(dims(f)))
    return only(matches)
end

function _cartesian_to_latlon(p::SVector{3})
    r = sqrt(sum(abs2, p))
    lat = asind(p[3] / r)
    lon = mod(rad2deg(atan(p[2], p[1])), 360)
    return lat, lon
end

function interpolate(f::DiscreteField{NodeLoc}, lat::Real, lon::Real)
    cid = locate_cell(mesh(f), lat, lon)
    nodes, weights = interpolation_weights(mesh(f), cid, lat, lon)
    values = data(f)
    loc_axis = _location_axis(NodeLoc, f)

    if ndims(values) == 1
        return sum(values[nodes[i]] * weights[i] for i in eachindex(nodes))
    end

    other_axes = filter(!=(loc_axis), 1:ndims(values))
    moved = PermutedDimsArray(values, (loc_axis, other_axes...))
    trailing_shape = size(moved)[2:end]
    out = zeros(promote_type(eltype(values), eltype(weights)), trailing_shape)
    for i in eachindex(nodes)
        out .+= weights[i] .* selectdim(moved, 1, nodes[i])
    end
    return out
end

function interpolate(f::DiscreteField{NodeLoc}, p::SVector{3})
    lat, lon = _cartesian_to_latlon(p)
    return interpolate(f, lat, lon)
end

"""
    interpolate(f::DiscreteField{NodeLoc}, lats, lons) -> DimArray

Bilinear node interpolation at many query points in one call. `lats` and `lons`
are equal-length vectors in degrees (`lon` is taken mod 360); the result is a
`DimArray` carrying the field name, a leading `Dim{:point}` axis and the
field's trailing axes unchanged.
"""
function interpolate(f::DiscreteField{NodeLoc}, lats::AbstractVector{<:Real},
        lons::AbstractVector{<:Real})
    length(lats) == length(lons) || throw(DimensionMismatch(
        "lats and lons must have equal length, got $(length(lats)) and $(length(lons))"
    ))
    n = length(lats)
    for (i, lat) in enumerate(lats)
        -90 <= lat <= 90 || throw(ArgumentError(
            "latitude $lat at index $i is outside [-90, 90]"
        ))
    end

    values = data(f)
    loc_axis = _location_axis(NodeLoc, f)
    trailing_dims = Tuple(d for (i, d) in enumerate(dims(f)) if i != loc_axis)
    T = promote_type(eltype(values), Float64)

    if ndims(values) == 1
        out = Vector{T}(undef, n)
        for (i, (lat, lon)) in enumerate(zip(lats, lons))
            out[i] = interpolate(f, lat, lon)
        end
    else
        out = Array{T}(undef, (n, map(length, trailing_dims)...))
        for (i, (lat, lon)) in enumerate(zip(lats, lons))
            copyto!(selectdim(out, 1, i), interpolate(f, lat, lon))
        end
    end

    return DimensionalData.DimArray(out, (Dim{:point}(1:n), trailing_dims...);
        name = DimensionalData.name(f))
end
