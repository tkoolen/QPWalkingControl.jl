struct CoMTrackingStateMachine{N, S<:SE3PDController}
    contacts::Dict{BodyID, Vector{QPControl.ContactPoint{N}}}
    end_effector_controllers::Dict{BodyID, S}
    # contact_plan::Piecewise{Constant{PlanarContactSituation{Float64}}, Float64, Vector{Constant{PlanarContactSituation{Float64}}}, Vector{Float64}}
    bodies::Vector{BodyID}
end

function CoMTrackingStateMachine(
            mechanism::Mechanism,
            contacts::Dict{BodyID, <:Vector{<:QPControl.ContactPoint}}
        )
    T = Float64
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

    CoMTrackingStateMachine(contacts, end_effector_controllers, bodies)
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
    # Enable/disable contacts
    # TODO
    for (bodyid, end_effector_controller) in statemachine.end_effector_controllers
        enable_contacts!(statemachine, bodyid)
        init_support!(end_effector_controller, 0.0) # TODO: tf?
    end

    # Gain scheduling / trajectory initialization
    # TODO
end
