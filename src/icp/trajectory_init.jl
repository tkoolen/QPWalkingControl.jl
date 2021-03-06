function transfer_weight!(
        generator::ICPTrajectoryGenerator{T}, state::MechanismState,
        foot_polygons::AbstractDict{BodyID, <:Framed{<:ConvexHull}}, bodyid::BodyID;
        Δt::Number) where {T}
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
    ξ0 = horizontal_projection(icp(state, generator.ω).v)
    generator.initial_cop[] = ξ0
    generator.initial_icp[] = ξ0
    generator.final_icp[] = unwrap(desired_icp)
    push_segment!(generator, Δt, contactmode.hull, unwrap(desired_icp))
    solve!(generator)
    nothing
end
