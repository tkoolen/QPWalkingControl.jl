struct PlanarContactSituation{T}
    active_bodies::Vector{BodyID}
    support_polygon::DConvexHull{T}
end

PlanarContactSituation{T}() where {T} = PlanarContactSituation(BodyID[], DConvexHull{T}())

struct ICPWalkingStateMachine{N, I<:ICPTrajectoryGenerator, S<:SE3PDController}
    contacts::Dict{BodyID, Vector{QPControl.ContactPoint{N}}}
    icp_trajectory_generator::I
    end_effector_controllers::Dict{BodyID, S}
    contact_plan::Piecewise{Constant{PlanarContactSituation{Float64}}, Float64, Vector{Constant{PlanarContactSituation{Float64}}}, Vector{Float64}}
    bodies::Vector{BodyID}
    # TODO: footstep generator
end

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
        gains = SE3PDGains(FramePDGains(bodyframe, PDGains(0.0, 0.0)), FramePDGains(bodyframe, PDGains(0.0, 0.0)))
        angulartraj = Constant(one(Quat)) # TODO
        lineartraj = BasicFootTrajectory(0.0, 1e-3, zero(SVector{3}), 1.0, zero(SVector{3}), 0.0)
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
    foot_polygons_world = typeof(foot_polygons)()
    for (bodyid, foot_polygon_sole) in foot_polygons
        foot_polygon_world = transform(foot_polygon_sole, transform_to_root(state, foot_polygon_sole.frame))
        foot_polygons_world[bodyid] = foot_polygon_world
    end
    Δt_first_transfer = 2.0
    Δt_transfer = 2.0
    Δt_swing = 1e-3#1.0
    step_length = 0.0#0.3

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
        support_points = collect(Iterators.flatten(vertices(unwrap(foot_polygons_world[bodyid])) for bodyid in support_bodies))
        jarvis_march!(contact_situation.support_polygon, support_points)
        empty!(contact_situation.active_bodies)
        append!(contact_situation.active_bodies, support_bodies)
        if length(support_bodies) == length(foot_polygons)
            Δt = first_step ? Δt_first_transfer : Δt_transfer
            delete!(support_bodies, next_swing_id)
        else
            Δt = Δt_swing
            foot_displacement = if first_step
                first_step = false
                SVector(step_length, 0.0)
            else
                SVector(2 * step_length, 0.0)
            end
            foot_polygon_world = foot_polygons_world[next_swing_id]
            foot_polygons_world[next_swing_id] = in_frame(foot_polygon_world.frame, ConvexHull{CCW}(vertices(unwrap(foot_polygon_world)) .+ Ref(foot_displacement)))
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
    initialize!(generator, p0, ξ0, ξf)
    for i in eachindex(contact_situations)
        Δt = ts[i + 1] - ts[i]
        push_segment!(generator, Δt, contact_situations[i].value.support_polygon)
    end
    solve!(statemachine.icp_trajectory_generator)
end

function (statemachine::ICPWalkingStateMachine)(t, state::MechanismState)
    active_bodies = statemachine.contact_plan(t).active_bodies
    for bodyid in statemachine.bodies
        if bodyid in active_bodies
            enable_contacts!(statemachine, bodyid)
        else
            disable_contacts!(statemachine, bodyid)
        end
    end
end
