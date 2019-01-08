Base.@propagate_inbounds function horizontal_projection(v::AbstractVector)
    @boundscheck length(v) == 3 || throw(ArgumentError())
    SVector(v[1], v[2])
end

critically_damped_gains(k::Number) = PDGains(k, 2 * sqrt(k))

function icp(c::Point3D, ċ::FreeVector3D, ω::Number)
    @framecheck c.frame ċ.frame
    c + ċ / ω
end
