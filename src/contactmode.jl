struct ContactMode{T}
    frame::CartesianFrame3D
    points_body::Dict{BodyID, Vector{SPoint3D{T}}} # TODO: consider removing (really only there to circumvent allocation)
    points_frame::Dict{BodyID, Vector{SPoint3D{T}}}
    projected_points_frame::Vector{SVector{2, T}}
    hull::FlexibleConvexHull{T}
end

function ContactMode{T}(frame::CartesianFrame3D, bodies) where T
    points_body = Dict{BodyID, Vector{SPoint3D{T}}}()
    points_frame = Dict{BodyID, Vector{SPoint3D{T}}}()
    for body in bodies
        points_body[BodyID(body)] = SPoint3D{T}[]
        points_frame[BodyID(body)] = SPoint3D{T}[]
    end
    projected_points_frame = SVector{2, T}[]
    hull = ConvexHull{CCW, T}()
    ContactMode{T}(frame, points_body, points_frame, projected_points_frame, hull)
end

_transform_to_root(state::MechanismState, bodyid::BodyID) = transform_to_root(state, bodyid)
_transform_to_root(transforms_to_root::AbstractDict{BodyID, <:Transform3D}, bodyid::BodyID) = transforms_to_root[bodyid]

function update!(mode::ContactMode{T}, transformsource::Union{MechanismState, AbstractDict{BodyID, <:Transform3D}};
        root_to_frame::Transform3D=one(Transform3D{T}, mode.frame)) where T
    empty!(mode.projected_points_frame)
    for (bodyid, points_body) in mode.points_body
        points_frame = mode.points_frame[bodyid]
        resize!(points_frame, length(points_body))
        to_root = _transform_to_root(transformsource, bodyid)
        to_frame = root_to_frame * to_root
        i = 1
        @inbounds for point_body in points_body
            points_frame[i] = to_frame * point_body
            push!(mode.projected_points_frame, horizontal_projection(points_frame[i].v))
            i += 1
        end
    end
    jarvis_march!(mode.hull, mode.projected_points_frame)
    mode
end
