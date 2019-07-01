function make_foot_polygons(mechanism::Mechanism{T}, soleframes::AbstractDict{BodyID, CartesianFrame3D},
        contact_points::AbstractDict{BodyID, <:Vector{<:Point3D}};
        num_extreme_points::Int) where T
    footpolygons = Dict{BodyID, Framed{SConvexHull{num_extreme_points, Float64}}}()
    for (bodyid, bodypoints) in contact_points
        body = findbody(mechanism, bodyid)
        soleframe = soleframes[bodyid]
        to_sole = inv(frame_definition(body, soleframe))
        projected_points = Vector{SVector{2, T}}(undef, length(bodypoints))
        for (i, point) in enumerate(bodypoints)
            point_sole_frame = transform(point, to_sole)
            @assert point_sole_frame.v[3] == 0
            projected_points[i] = horizontal_projection(point_sole_frame.v)
        end
        flexiblehull = DConvexHull{T}()
        jarvis_march!(flexiblehull, projected_points)
        @assert num_vertices(flexiblehull) == num_extreme_points
        footpolygons[bodyid] = in_frame(soleframe, SConvexHull{num_extreme_points, T}(flexiblehull))
    end
    footpolygons
end
