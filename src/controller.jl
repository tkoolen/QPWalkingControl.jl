struct PushRecoveryController{M<:MomentumBasedController, L}
    lowlevel::M
    robotmass::Float64
    gravitymag::Float64

    foottasks::Dict{RigidBody{Float64}, SpatialAccelerationTask}
    linmomtask::LinearMomentumRateTask
    pelvistask::AngularAccelerationTask
    jointtasks::Dict{JointID, JointAccelerationTask{Revolute{Float64}}}

    linear_momentum_controller::L
    pelvisgains::PDGains{Float64,Float64}
    jointgains::Dict{JointID, PDGains{Float64, Float64}}

    comref::Point3D{SVector{3, Float64}}
    jointrefs::Dict{JointID, Float64}

    contactmode::ContactMode{Float64}
end

function PushRecoveryController(
        lowlevel::MomentumBasedController,
        feet::Vector{<:RigidBody},
        pelvis::RigidBody,
        nominalstate::MechanismState,
        linear_momentum_controller;
        joint_regularization::Float64 = 0.05,
        linear_momentum_weight::Float64 = 1.0,
        pelvisgains::PDGains = PDGains(20., 2 * sqrt(20.0)),
        jointgains = Dict(JointID(j) => PDGains(100.0, 20.) for j in tree_joints(lowlevel.state.mechanism)),
        comref::Point3D = center_of_mass(nominalstate) - FreeVector3D(root_frame(lowlevel.state.mechanism), 0., 0., 0.05),
        jointrefs = Dict(JointID(j) => configuration(nominalstate, j)[1] for j in tree_joints(lowlevel.state.mechanism) if joint_type(j) isa Revolute))
    mechanism = lowlevel.state.mechanism
    world = root_body(mechanism)
    worldframe = root_frame(mechanism)
    m = mass(mechanism)
    g = norm(mechanism.gravitational_acceleration)

    regularize!.(Ref(lowlevel), tree_joints(mechanism), joint_regularization)

    foottasks = Dict(foot => SpatialAccelerationTask(mechanism, path(mechanism, world, foot)) for foot in feet)
    addtask!.(Ref(lowlevel), collect(values(foottasks)))

    linmomtask = LinearMomentumRateTask(mechanism, centroidal_frame(lowlevel))
    addtask!(lowlevel, linmomtask, linear_momentum_weight)

    pelvistask = AngularAccelerationTask(mechanism, path(mechanism, world, pelvis))
    addtask!(lowlevel, pelvistask, 10.0)

    revolutejoints = filter(j -> joint_type(j) isa Revolute, tree_joints(mechanism))
    positioncontroljoints = setdiff(revolutejoints, vcat((collect(task.path) for task in values(foottasks))...))
    jointtasks = Dict(JointID(j) => JointAccelerationTask(j) for j in positioncontroljoints)
    addtask!.(Ref(lowlevel), collect(values(jointtasks)), 1.0)

    contactmode = ContactMode{Float64}(worldframe, keys(lowlevel.contacts))

    PushRecoveryController(
        lowlevel, m, g,
        foottasks, linmomtask, pelvistask, jointtasks,
        linear_momentum_controller, pelvisgains, jointgains,
        comref, jointrefs,
        contactmode)
end

function (controller::PushRecoveryController)(τ::AbstractVector, t::Number, state::MechanismState)
    # Update contact mode
    contactmode = controller.contactmode
    update_active_contacts!(contactmode, controller.lowlevel.contacts, state)

    # State transitions

    # Linear momentum control
    c = center_of_mass(state)
    h = momentum(state)
    ċ = FreeVector3D(h.frame, linear(h) / controller.robotmass)
    l̇_des = controller.linear_momentum_controller(t, c, ċ, contactmode.hull)
    world_to_centroidal = Transform3D(c.frame, centroidal_frame(controller.lowlevel), -c.v)
    l̇_des = transform(l̇_des, world_to_centroidal)
    setdesired!(controller.linmomtask, l̇_des)

    # Pelvis orientation control
    pelvis = target(controller.pelvistask.path)
    Hpelvis = transform_to_root(state, pelvis)
    Tpelvis = transform(twist_wrt_world(state, pelvis), inv(Hpelvis))
    ωd_pelvis_des = FreeVector3D(Tpelvis.frame, pd(controller.pelvisgains, rotation(Hpelvis), Tpelvis.angular))
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

function update_active_contacts!(mode::ContactMode, contacts, state::MechanismState)
    for (body, contactpoints) in contacts
        bodyid = BodyID(body)
        active_points_body = mode.active_points_body[bodyid]
        empty!(active_points_body)
        for contactpoint in contactpoints
            if QPControl.isenabled(contactpoint)
                push!(active_points_body, contactpoint.position)
            end
        end
    end
    update!(mode, state)
    nothing
end
