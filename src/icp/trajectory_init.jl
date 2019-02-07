# Wrong. Just need to find desired ICP for one shift, use centroid as desired CoP.

function init_back_and_forth_shift!(
        generator::ICPTrajectoryGenerator{T}, state::MechanismState,
        soleframes::AbstractDict{BodyID, CartesianFrame3D},
        foot_polygons_sole_frame::AbstractDict{BodyID, <:ConvexHull};
        num_shifts::Integer, Δt = 1.0) where T
    # Compute ω (assumes ground is at z = 0)
    mechanism = state.mechanism
    gz = norm(mechanism.gravitational_acceleration)
    z = center_of_mass(state).v[3]
    ω = sqrt(abs(gz / z))

    # Transform foot polygons to world frame and set up contact mode
    foot_polygons_world_frame = typeof(foot_polygons_sole_frame)()
    contactmode = PlanarContactMode{T}(root_frame(mechanism))
    for (bodyid, hull_sole_frame) in foot_polygons_sole_frame
        soleframe = soleframes[bodyid]
        tf = transform_to_root(state, soleframe)
        R = horizontal_projection(rotation(tf))
        p = horizontal_projection(translation(tf))
        hull_world_frame = typeof(hull_sole_frame)(map(x -> R * x + p, vertices(hull_sole_frame)))
        foot_polygons_world_frame[bodyid] = hull_world_frame
        for point in vertices(hull_world_frame)
            push!(contactmode.projected_points, point)
        end
    end
    update!(contactmode)

    # First ICP knot point is ICP at given state
    generator.initial_icp[] = horizontal_projection(icp(state, ω).v)

    # Subsequent ICPs cycle between centroids of foot polygons in world frame
    for (i, bodyid) in enumerate(cycle(keys(foot_polygons_world_frame))) # TODO: order undefined
        i > num_shifts && break
        hull_world_frame = foot_polygons_world_frame[bodyid]
        t, icp = last(icpknots)
        next_icp = centroid(hull_world_frame)
        push!(icpknots, t + Δt => next_icp)
    end

    # Initialize trajectory generator
    empty!(generator)
    for i = 1 : num_active_segments
        push_segment!(generator, Δt, ω, foot_polygon_hrep + foot_centers[i], foot_centers[i])
    end
    generator.initial_icp[] = foot_centers[1] + SVector(0.02, 0.01)
    generator.final_icp[] = foot_centers[num_active_segments]

    initialize!(generator, Constant(contactmode), icpknots, fill(ω, length(icpknots) - 1))
end
