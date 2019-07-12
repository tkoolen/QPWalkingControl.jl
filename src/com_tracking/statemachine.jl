struct CoMTrackingStateMachine{N, S<:SE3PDController}
    contacts::Dict{BodyID, Vector{QPControl.ContactPoint{N}}}
    end_effector_controllers::Dict{BodyID, S}
    pose_plans::Dict{BodyID, PosePlan{Float64}}
    bodyids::Vector{BodyID}
end

function CoMTrackingStateMachine(
            mechanism::Mechanism,
            contacts::Dict{BodyID, <:Vector{<:QPControl.ContactPoint}}
        )
    T = Float64
    bodyids = sort!(collect(keys(contacts)), by = x -> x.value)
    end_effector_controllers = map(bodyids) do bodyid
        body = findbody(mechanism, bodyid)
        bodyframe = default_frame(body)
        baseframe = root_frame(mechanism)
        gains = SE3PDGains(FramePDGains(bodyframe, PDGains(0.0, 20.0)), FramePDGains(bodyframe, PDGains(0.0, 0.0)))
        angulartraj = interpolated_orientation_trajectory(0.0, 1e-6, one(Quat{Float64}), one(Quat{Float64})) # TODO: 1e-6
        lineartraj = convert(BasicFootTrajectory{T}, SVector(0.0, 0.0, 0.0), 0.0)
        trajectory = SE3Trajectory(bodyframe, baseframe, angulartraj, lineartraj)
        weight = Diagonal(vcat(fill(10.0, 3), fill(10.0, 3)))
        bodyid => SE3PDController(BodyID(root_body(mechanism)), BodyID(body), trajectory, weight, gains)
    end |> Dict
    pose_plans = Dict(bodyid => PosePlan{Float64}() for bodyid in bodyids)
    ret = CoMTrackingStateMachine(contacts, end_effector_controllers, pose_plans, bodyids)
    for bodyid in bodyids
        enable_contacts!(ret, bodyid)
    end
    return ret
end

function set_pose_plan!(statemachine::CoMTrackingStateMachine, bodyid::BodyID, plan::PosePlan)
    statemachine.pose_plans[bodyid] = plan
    return nothing
end

function enable_contacts!(statemachine::CoMTrackingStateMachine, bodyid::BodyID)
    for contact in statemachine.contacts[bodyid]
        contact.maxnormalforce[] = 1e6 # TODO
    end
end

function disable_contacts!(statemachine::CoMTrackingStateMachine, bodyid::BodyID)
    for contact in statemachine.contacts[bodyid]
        disable!(contact)
    end
end

function (statemachine::CoMTrackingStateMachine)(t, state::MechanismState)
    end_effector_controllers = statemachine.end_effector_controllers
    pose_plans = statemachine.pose_plans
    for bodyid in statemachine.bodyids
        end_effector_controller = end_effector_controllers[bodyid]
        if isdone(end_effector_controller.trajectory[], t)
            pose_plan = pose_plans[bodyid]
            if !isempty(pose_plan) && t >= next_move_start_time(pose_plan) && !statemachine.in_contact[bodyid]
                pose0 = transform_to_root(state, bodyid)
                t0, duration, posef = popfirst!(pose_plan)
                disable_contacts!(statemachine, bodyid)
                init_swing!(end_effector_controller, pose0, posef, t0=t0, tf=t0 + duration)
            else
                enable_contacts!(statemachine, bodyid)
                init_support!(end_effector_controller; t0=t, tf=t + 1e-6) # TODO: 1e-6
            end
        end
    end
end
