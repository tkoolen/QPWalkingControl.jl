# TODO: this is very messy

struct PlanarContactSituation{T}
    active_body_poses::Dict{BodyID, Transform3D{T}}
    support_polygon::DConvexHull{T}
end

PlanarContactSituation{T}() where {T} = PlanarContactSituation(Dict{BodyID, Transform3D{T}}(), DConvexHull{T}())

struct ICPWalkingStateMachine{N, I<:ICPTrajectoryGenerator, S<:SE3PDController}
    contacts::Dict{BodyID, Vector{QPControl.ContactPoint{N}}}
    icp_trajectory_generator::I
    end_effector_controllers::Dict{BodyID, S}
    contact_plan::Piecewise{Constant{PlanarContactSituation{Float64}}, Float64, Vector{Constant{PlanarContactSituation{Float64}}}, Vector{Float64}}
    bodies::Vector{BodyID}
    # TODO: footstep generator
end

# TODO: not super nice:
isdone(trajectory::SE3Trajectory, t) = isdone(trajectory.angular, t) || isdone(trajectory.linear, t)
isdone(trajectory::Constant, t) = false
isdone(trajectory::BasicFootTrajectory, t) = t > trajectory.tf

function ICPWalkingStateMachine(
            mechanism::Mechanism,
            contacts::Dict{BodyID, <:Vector{<:QPControl.ContactPoint}},
            icp_trajectory_generator::ICPTrajectoryGenerator
        )
    T = Float64
    n = num_segments(icp_trajectory_generator)
    contact_situations = [PlanarContactSituation{Float64}() for i = 1 : n]
    times = zeros(T, n + 1)
    contact_plan = Piecewise(Constant.(contact_situations), times, clamp=true)
    bodies = sort!(collect(keys(contacts)), by = x -> x.value)

    # create end effector controllers
    end_effector_controllers = map(bodies) do bodyid
        body = findbody(mechanism, bodyid)
        bodyframe = default_frame(body)
        baseframe = root_frame(mechanism)
        gains = SE3PDGains(FramePDGains(bodyframe, PDGains(0.0, 20.0)), FramePDGains(bodyframe, PDGains(0.0, 0.0)))
        angulartraj = Constant(one(Quat)) # TODO
        lineartraj = convert(BasicFootTrajectory{T}, SVector(0.0, 0.0, 0.0), 0.0)
        trajectory = SE3Trajectory(bodyframe, baseframe, angulartraj, lineartraj)
        weight = Diagonal(vcat(fill(10.0, 3), fill(10.0, 3)))
        BodyID(body) => SE3PDController(BodyID(root_body(mechanism)), BodyID(body), trajectory, weight, gains)
    end |> Dict

    statemachine = ICPWalkingStateMachine(
        contacts, icp_trajectory_generator, end_effector_controllers,
        contact_plan, bodies)
end

function enable_contacts!(statemachine::ICPWalkingStateMachine, bodyid::BodyID)
    for contact in statemachine.contacts[bodyid]
        contact.maxnormalforce[] = 1e6 # TODO
    end
end

function disable_contacts!(statemachine::ICPWalkingStateMachine, bodyid::BodyID)
    for contact in statemachine.contacts[bodyid]
        disable!(contact)
    end
end

# TODO: temporary
function init_footstep_plan!(statemachine::ICPWalkingStateMachine, state::MechanismState, foot_polygons)
    T = Float64
    body_poses = Dict{BodyID, Transform3D{T}}()
    for bodyid in keys(foot_polygons)
        body_pose = transform_to_root(state, bodyid)
        body_poses[bodyid] = body_pose
    end
    Δt_first_transfer = 2.0
    Δt_transfer = 1.0
    Δt_swing = 0.8
    step_length = 0.3

    ts = breaks(statemachine.contact_plan)
    ts[1] = 0
    contact_situations = subfunctions(statemachine.contact_plan)
    n = length(contact_situations)
    support_bodies = Set(keys(foot_polygons))
    next_swing_id_iterator = Iterators.cycle(statemachine.bodies)
    next_swing_id, next_swing_id_state = iterate(next_swing_id_iterator)
    first_step = true

    for i = 1 : n
        contact_situation = contact_situations[i].value
        foot_polygons_world = valtype(foot_polygons)[]
        for bodyid in support_bodies
            foot_polygon = foot_polygons[bodyid]
            body = findbody(state.mechanism, bodyid)
            sole_to_world = body_poses[bodyid] * frame_definition(body, foot_polygon.frame)
            push!(foot_polygons_world, transform(foot_polygon, sole_to_world))
        end
        support_points = collect(Iterators.flatten(vertices(unwrap(foot_polygon_world)) for foot_polygon_world in foot_polygons_world))
        jarvis_march!(contact_situation.support_polygon, support_points)
        empty!(contact_situation.active_body_poses)
        for body in support_bodies
            push!(contact_situation.active_body_poses, body => body_poses[body])
        end
        if length(support_bodies) == length(foot_polygons)
            Δt = first_step ? Δt_first_transfer : Δt_transfer
            delete!(support_bodies, next_swing_id)
        else
            Δt = Δt_swing
            foot_displacement = if first_step
                first_step = false
                SVector(step_length, 0.0, 0.0)
            else
                SVector(2 * step_length, 0.0, 0.0)
            end
            foot_pose = body_poses[next_swing_id]
            body_poses[next_swing_id] = Transform3D(foot_pose.from, foot_pose.to, rotation(foot_pose), translation(foot_pose) + foot_displacement)
            next_swing_id, next_swing_id_state = iterate(next_swing_id_iterator, next_swing_id_state)
            for bodyid in keys(foot_polygons)
                push!(support_bodies, bodyid)
            end
        end
        ts[i + 1] = ts[i] + Δt
    end
    ξ0 = horizontal_projection(icp(state).v)
    p0 = ξ0
    ξf = centroid(last(contact_situations).value.support_polygon)
    generator = statemachine.icp_trajectory_generator
    empty!(generator)
    initialize!(generator, p0, ξ0, ξf)
    for i in eachindex(contact_situations)
        Δt = ts[i + 1] - ts[i]
        push_segment!(generator, Δt, contact_situations[i].value.support_polygon)
    end
    solve!(statemachine.icp_trajectory_generator)

    for end_effector_controller in values(statemachine.end_effector_controllers)
        init_support!(end_effector_controller, ts[2])
    end
end

function init_support!(end_effector_controller::SE3PDController, tf::Number)
    # TODO: frame lookup is kind of nasty
    bodyframe = end_effector_controller.trajectory[].body
    baseframe = end_effector_controller.trajectory[].base
    end_effector_controller.gains[] = SE3PDGains(FramePDGains(bodyframe, PDGains(0.0, 20.0)), FramePDGains(bodyframe, PDGains(0.0, 0.0)))
    angulartraj = Constant(one(Quat)) # TODO
    lineartraj = convert(BasicFootTrajectory{Float64}, SVector(0.0, 0.0, 0.0), tf)
    end_effector_controller.trajectory[] = SE3Trajectory(bodyframe, baseframe, angulartraj, lineartraj)
end

function init_swing!(end_effector_controller::SE3PDController, pose0::Transform3D, posef::Transform3D, t0::Number, tf::Number)
    # FIXME: poses are in sole frame!
    # TODO: frame lookup is kind of nasty
    bodyframe = end_effector_controller.trajectory[].body
    baseframe = end_effector_controller.trajectory[].base
    end_effector_controller.gains[] = SE3PDGains(FramePDGains(bodyframe, PDGains(100.0, 20.0)), FramePDGains(bodyframe, PDGains(100.0, 0.0)))
    angulartraj = Constant(one(Quat)) # TODO
    Δzmid = 0.1
    lineartraj = BasicFootTrajectory(t0, tf, translation(pose0), Δzmid, translation(posef), -0.1)
    end_effector_controller.trajectory[] = SE3Trajectory(bodyframe, baseframe, angulartraj, lineartraj)
end

function (statemachine::ICPWalkingStateMachine)(t, state::MechanismState)
    active_body_poses = statemachine.contact_plan(t).active_body_poses

    # Enable/disable contacts
    for bodyid in statemachine.bodies
        if bodyid in keys(active_body_poses)
            enable_contacts!(statemachine, bodyid)
        else
            disable_contacts!(statemachine, bodyid)
        end
    end

    # Gain scheduling / trajectory initialization
    contact_plan = statemachine.contact_plan
    for (bodyid, end_effector_controller) in statemachine.end_effector_controllers
        if isdone(end_effector_controller.trajectory[], t)
            t′ = clamp(t, first(breaks(contact_plan)), last(breaks(contact_plan)))
            n = length(subfunctions(contact_plan))
            index = min(searchsortedlast(breaks(contact_plan), t′), n)
            tf = breaks(contact_plan)[min(index + 1, n)]
            if bodyid in keys(active_body_poses)
                init_support!(end_effector_controller, tf)
            else
                previous_pose = subfunctions(contact_plan)[index - 1].value.active_body_poses[bodyid]
                next_pose = subfunctions(contact_plan)[index + 1].value.active_body_poses[bodyid]
                init_swing!(end_effector_controller, previous_pose, next_pose, t′, tf)
            end
        end
    end
end
