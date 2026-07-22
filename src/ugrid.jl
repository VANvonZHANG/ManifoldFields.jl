import NCDatasets

struct UGridDataset
    variables::Dict{String,Any}
    attributes::Dict{String,Any}
end

function to_ugrid(args...; kwargs...)
    throw(ErrorException("to_ugrid is implemented in a later task"))
end

function from_ugrid(args...; kwargs...)
    throw(ErrorException("from_ugrid is implemented in a later task"))
end

function from_ugrid_mesh(args...; kwargs...)
    throw(ErrorException("from_ugrid_mesh is implemented in a later task"))
end

function save_ugrid(args...; kwargs...)
    throw(ErrorException("save_ugrid is implemented in a later task"))
end

function load_ugrid(args...; kwargs...)
    throw(ErrorException("load_ugrid is implemented in a later task"))
end

function load_ugrid_mesh(args...; kwargs...)
    throw(ErrorException("load_ugrid_mesh is implemented in a later task"))
end
