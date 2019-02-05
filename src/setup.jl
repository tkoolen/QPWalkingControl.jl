function add_sole_frames!(mechanism::Mechanism)
    soleframes = Dict{BodyID, CartesianFrame3D}()
    for body in bodies(mechanism)
        bodypoints = RigidBodyDynamics.contact_points(body)
        if !isempty(bodypoints)
            bodyid = BodyID(body)
            bodyframe = default_frame(body)
            soleheight = missing
            for point in bodypoints
                position = location(point)
                @framecheck position.frame bodyframe
                if soleheight === missing
                    soleheight = position.v[3]
                else
                    position.v[3] == soleheight || throw(ArgumentError("Can't currently handle non-z-normal sole planes."))
                end
            end
            soleframe = soleframes[bodyid] = CartesianFrame3D("$(body)_sole")
            sole_to_body = Transform3D(soleframe, bodyframe, SVector(0, 0, soleheight))
            add_frame!(body, sole_to_body)
        end
    end
    soleframes
end

function foot_polygons(mechanism::Mechanism{T}, soleframes::AbstractDict{BodyID, CartesianFrame3D};
        num_extreme_points::Int) where T
    footpolygons = Dict{BodyID, SConvexHull{num_extreme_points, Float64}}()
    for body in bodies(mechanism)
        bodypoints = RigidBodyDynamics.contact_points(body)
        if !isempty(bodypoints)
            bodyid = BodyID(body)
            to_sole = inv(frame_definition(body, soleframes[bodyid]))
            projected_points = Vector{SVector{2, T}}(undef, length(bodypoints))
            for (i, point) in enumerate(bodypoints)
                point_sole_frame = transform(location(point), to_sole)
                @assert point_sole_frame.v[3] == 0
                projected_points[i] = horizontal_projection(point_sole_frame.v)
            end
            flexiblehull = DConvexHull{T}()
            jarvis_march!(flexiblehull, projected_points)
            @assert num_vertices(flexiblehull) == num_extreme_points
            footpolygons[bodyid] = SConvexHull{num_extreme_points, T}(flexiblehull)
        end
    end
    footpolygons
end
