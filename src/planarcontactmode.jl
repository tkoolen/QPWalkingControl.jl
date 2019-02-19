# TODO: remove
struct PlanarContactMode{T}
    frame::CartesianFrame3D
    projected_points::Vector{Vec2{T}}
    hull::DConvexHull{T}
end

function PlanarContactMode{T}(frame::CartesianFrame3D) where T
    projected_points = Vec2{T}[]
    hull = DConvexHull{T}()
    PlanarContactMode{T}(frame, projected_points, hull)
end

_transform_to_root(state::MechanismState, bodyid::BodyID) = transform_to_root(state, bodyid)
_transform_to_root(transforms_to_root::AbstractDict{BodyID, <:Transform3D}, bodyid::BodyID) = transforms_to_root[bodyid]

function update!(mode::PlanarContactMode{T}, active_points::Dict{BodyID, Vector{SPoint3D{T}}},
        transformsource::Union{MechanismState, AbstractDict{BodyID, <:Transform3D}};
        root_to_frame::Transform3D=one(Transform3D{T}, mode.frame)) where T
    empty!(mode.projected_points)
    for (bodyid, points_body) in active_points
        to_root = _transform_to_root(transformsource, bodyid)
        to_frame = root_to_frame * to_root
        for point_body in points_body
            point_frame = to_frame * point_body
            push!(mode.projected_points, horizontal_projection(point_frame.v))
        end
    end
    update!(mode)
    mode
end

function update!(mode::PlanarContactMode)
    jarvis_march!(mode.hull, mode.projected_points)
end
