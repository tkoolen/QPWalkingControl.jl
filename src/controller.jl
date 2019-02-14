struct PushRecoveryController{M<:MomentumBasedController, L}
    lowlevel::M
    robotmass::Float64
    gravitymag::Float64

    foot_angular_tasks::Dict{BodyID, AngularAccelerationTask}
    foot_linear_tasks::Dict{BodyID, LinearAccelerationTask}
    linmomtask::LinearMomentumRateTask
    pelvistask::AngularAccelerationTask
    jointtasks::Dict{JointID, JointAccelerationTask{Revolute{Float64}}}

    linear_momentum_controller::L
    pelvisgains::PDGains{Float64,Float64}
    jointgains::Dict{JointID, PDGains{Float64, Float64}}

    comref::Point3D{SVector{3, Float64}}
    jointrefs::Dict{JointID, Float64}

    active_contact_points::Dict{BodyID, Vector{SPoint3D{Float64}}}
end

function PushRecoveryController(
        lowlevel::MomentumBasedController,
        feet::Vector{<:RigidBody},
        pelvis::RigidBody,
        nominalstate::MechanismState,
        linear_momentum_controller;
        joint_regularization::Float64 = 0.05,
        linear_momentum_weight::Float64 = 1.0,
        pelvisgains::PDGains = critically_damped_gains(100.0),
        jointgains = Dict(JointID(j) => critically_damped_gains(100.0) for j in tree_joints(lowlevel.state.mechanism)),
        comref::Point3D = center_of_mass(nominalstate) - FreeVector3D(root_frame(lowlevel.state.mechanism), 0., 0., 0.05),
        jointrefs = Dict(JointID(j) => configuration(nominalstate, j)[1] for j in tree_joints(lowlevel.state.mechanism) if joint_type(j) isa Revolute))
    mechanism = lowlevel.state.mechanism
    world = root_body(mechanism)
    worldframe = root_frame(mechanism)
    m = mass(mechanism)
    g = norm(mechanism.gravitational_acceleration)

    regularize!.(Ref(lowlevel), tree_joints(mechanism), joint_regularization)

    foot_paths = Dict(BodyID(foot) => path(mechanism, world, foot) for foot in feet)

    foot_angular_tasks = Dict(BodyID(foot) => AngularAccelerationTask(mechanism, foot_paths[BodyID(foot)]) for foot in feet)
    addtask!.(Ref(lowlevel), collect(values(foot_angular_tasks)), 10.0)

    foot_linear_tasks = Dict(BodyID(foot) => LinearAccelerationTask(mechanism, foot_paths[BodyID(foot)]) for foot in feet)
    addtask!.(Ref(lowlevel), collect(values(foot_linear_tasks)))

    linmomtask = LinearMomentumRateTask(mechanism, centroidal_frame(lowlevel))
    addtask!(lowlevel, linmomtask, linear_momentum_weight)

    pelvistask = AngularAccelerationTask(mechanism, path(mechanism, world, pelvis))
    addtask!(lowlevel, pelvistask, 10.0)

    revolutejoints = filter(j -> joint_type(j) isa Revolute, tree_joints(mechanism))
    positioncontroljoints = setdiff(revolutejoints, vcat((collect(path) for path in values(foot_paths))...))
    jointtasks = Dict(JointID(j) => JointAccelerationTask(j) for j in positioncontroljoints)
    addtask!.(Ref(lowlevel), collect(values(jointtasks)), 10.0)

    active_contact_points = Dict{BodyID, Vector{SPoint3D{Float64}}}()

    PushRecoveryController(
        lowlevel, m, g,
        foot_angular_tasks, foot_linear_tasks, linmomtask, pelvistask, jointtasks,
        linear_momentum_controller, pelvisgains, jointgains,
        comref, jointrefs,
        active_contact_points)
end

function (controller::PushRecoveryController)(τ::AbstractVector, t::Number, state::MechanismState)
    # Update active contact points
    update_active_contacts!(controller.active_contact_points, controller.lowlevel.contacts)

    # Foot angular control
    for (bodyid, task) in controller.foot_angular_tasks
        # TODO: trajectory tracking
        H = transform_to_root(state, bodyid)
        T = transform(twist_wrt_world(state, bodyid), inv(H))
        ωd_des = FreeVector3D(T.frame, pd(PDGains(0.0, 20.0), rotation(H), angular(T)))
        setdesired!(task, ωd_des)
    end

    # Linear momentum control
    c = center_of_mass(state)
    h = momentum(state)
    ċ = FreeVector3D(h.frame, linear(h) / controller.robotmass)
    l̇_des = controller.linear_momentum_controller(t, c, ċ, controller.active_contact_points, state)
    world_to_centroidal = Transform3D(c.frame, centroidal_frame(controller.lowlevel), -c.v)
    l̇_des = transform(l̇_des, world_to_centroidal)
    setdesired!(controller.linmomtask, l̇_des)

    # Pelvis orientation control
    pelvis = target(controller.pelvistask.path)
    Hpelvis = transform_to_root(state, pelvis)
    Tpelvis = transform(twist_wrt_world(state, pelvis), inv(Hpelvis))
    ωd_pelvis_des = FreeVector3D(Tpelvis.frame, pd(controller.pelvisgains, rotation(Hpelvis), angular(Tpelvis)))
    setdesired!(controller.pelvistask, ωd_pelvis_des)

    # Joint position control
    for jointid in keys(controller.jointtasks)
        task = controller.jointtasks[jointid]
        gains = controller.jointgains[jointid]
        ref = controller.jointrefs[jointid]
        v̇_joint_des = pd(gains, configuration(state, jointid)[1], ref, velocity(state, jointid)[1], 0.0)
        setdesired!(task, v̇_joint_des)
    end

    controller.lowlevel(τ, t, state)
    τ
end

function update_active_contacts!(active_contact_points::AbstractDict{BodyID, <:AbstractVector{<:Point3D}}, contacts)
    for active_points in values(active_contact_points)
        empty!(active_points)
    end
    for (body, points) in contacts
        bodyid = BodyID(body)
        active_points_body = get!(valtype(active_contact_points), active_contact_points, bodyid)
        for contactpoint in points
            if QPControl.isenabled(contactpoint)
                push!(active_points_body, contactpoint.position)
            end
        end
    end
    nothing
end
