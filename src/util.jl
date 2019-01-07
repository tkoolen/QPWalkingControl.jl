Base.@propagate_inbounds function horizontal_projection(v::AbstractVector)
    @boundscheck length(v) == 3 || throw(ArgumentError())
    SVector(v[1], v[2])
end
