struct CoMTrackingStateMachine{N, S<:SE3PDController}
    in_contact::Dict{BodyID, Bool}
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
    in_contact = Dict(bodyid => false for bodyid in keys(contacts))
    bodyids = sort!(collect(keys(contacts)), by = x -> x.value)
    end_effector_controllers = map(bodyids) do bodyid
        body = findbody(mechanism, bodyid)
        bodyframe = default_frame(body)
        baseframe = root_frame(mechanism)
        gains = SE3PDGains(FramePDGains(bodyframe, PDGains(0.0, 20.0)), FramePDGains(bodyframe, PDGains(0.0, 0.0)))
        angulartraj = interpolated_orientation_trajectory(0.0, Inf, one(Quat{Float64}), one(Quat{Float64}))
        lineartraj = convert(BasicFootTrajectory{T}, SVector(0.0, 0.0, 0.0), Inf)
        trajectory = SE3Trajectory(bodyframe, baseframe, angulartraj, lineartraj)
        weight = Diagonal(vcat(fill(10.0, 3), fill(10.0, 3)))
        bodyid => SE3PDController(BodyID(root_body(mechanism)), BodyID(body), trajectory, weight, gains)
    end |> Dict
    pose_plans = Dict(bodyid => PosePlan{Float64}() for bodyid in bodyids)
    ret = CoMTrackingStateMachine(in_contact, contacts, end_effector_controllers, pose_plans, bodyids)
    for bodyid in bodyids
        enable_contacts!(ret, bodyid)
    end
    return ret
end

function set_pose_plan!(statemachine::CoMTrackingStateMachine, bodyid::BodyID, plan::PosePlan)
    statemachine.pose_plans[bodyid] = plan
    end_effector_controller = statemachine.end_effector_controllers[bodyid]
    init_support!(end_effector_controller; t0=0.0, tf=next_move_start_time(plan))
    return nothing
end

function enable_contacts!(statemachine::CoMTrackingStateMachine, bodyid::BodyID)
    statemachine.in_contact[bodyid] = true
    for contact in statemachine.contacts[bodyid]
        contact.maxnormalforce[] = 1e6 # TODO
    end
end

function disable_contacts!(statemachine::CoMTrackingStateMachine, bodyid::BodyID)
    statemachine.in_contact[bodyid] = false
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
            if !isempty(pose_plan)
                if t >= next_move_start_time(pose_plan)
                    # starting swing
                    disable_contacts!(statemachine, bodyid)
                    @assert !isempty(pose_plan)
                    pose0 = transform_to_root(state, bodyid)
                    t0, duration, posef = popfirst!(pose_plan)
                    body_to_sole = inv(frame_definition(findbody(state.mechanism, bodyid), posef.from)) # TODO: inefficient
                    posef = posef * body_to_sole
                    println(findbody(state.mechanism, bodyid), " entering swing at $t. Goal: $posef")
                    init_swing!(end_effector_controller, pose0, posef; t0=t0, tf=t0 + duration, zdf=0.0, Δzmid=0.15)
                else
                    # coming into contact
                    let
                        Href, _, _ = end_effector_controller.trajectory[](t, Val(2))
                        H = relative_transform(state, Href.from, Href.to)
                        Herr = inv(Href) * H
                        println(findbody(state.mechanism, bodyid), " tracking error at end of step: ", Herr)
                    end
                    enable_contacts!(statemachine, bodyid)
                    println(findbody(state.mechanism, bodyid), " entering support at $t.")
                    init_support!(end_effector_controller; t0=t, tf=next_move_start_time(pose_plan))
                end
            else
                # entering final contact phase
                enable_contacts!(statemachine, bodyid)
                println(findbody(state.mechanism, bodyid), " entering final support at $t.")
                init_support!(end_effector_controller; t0=t, tf=Inf)
            end
        end
    end
end
