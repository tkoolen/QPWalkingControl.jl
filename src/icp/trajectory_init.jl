function transfer_weight!(
        generator::ICPTrajectoryGenerator{T}, state::MechanismState,
        foot_polygons::AbstractDict{BodyID, <:Framed{<:ConvexHull}}, bodyid::BodyID;
        Δt::Number, ω::Number) where {T}
    foot_polygons_world = typeof(foot_polygons)()
    contactmode = PlanarContactMode{T}(root_frame(state.mechanism))
    for (bodyid, foot_polygon_sole) in foot_polygons
        foot_polygon_world = transform(foot_polygon_sole, transform_to_root(state, foot_polygon_sole.frame))
        foot_polygons_world[bodyid] = foot_polygon_world
        for point in vertices(unwrap(foot_polygon_world))
            push!(contactmode.projected_points, point)
        end
    end
    update!(contactmode)
    foot_polygon_world = foot_polygons_world[bodyid]
    desired_icp = @framechecked centroid(foot_polygon_world)
    @framecheck desired_icp.frame contactmode.frame
    empty!(generator)
    generator.initial_icp[] = horizontal_projection(icp(state, ω).v)
    generator.final_icp[] = unwrap(desired_icp)
    push_segment!(generator, Δt, ω, contactmode.hull, unwrap(desired_icp))
    solve!(generator)
    nothing
end
