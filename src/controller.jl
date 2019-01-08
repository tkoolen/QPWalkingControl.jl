struct PushRecoveryController{M<:MomentumBasedController, L, C<:ConvexHullProblem}
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

    support_polygon_problem::C
end

function PushRecoveryController(
        lowlevel::MomentumBasedController,
        feet::Vector{<:RigidBody},
        pelvis::RigidBody,
        nominalstate::MechanismState,
        convexhulloptimizer::MOI.AbstractOptimizer;
        joint_regularization::Float64 = 0.05,
        linear_momentum_weight::Float64 = 1.0,
        linear_momentum_controller = ICPController(lowlevel.state.mechanism),
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

    num_contacts = sum(length, values(lowlevel.contacts))
    support_polygon_problem = ConvexHullProblem{2, num_contacts, Float64}(convexhulloptimizer)

    PushRecoveryController(
        lowlevel, m, g,
        foottasks, linmomtask, pelvistask, jointtasks,
        linear_momentum_controller, pelvisgains, jointgains,
        comref, jointrefs,
        support_polygon_problem)
end

function (controller::PushRecoveryController)(τ::AbstractVector, t::Number, state::MechanismState)
    # Support polygon setup
    set_contact_points!(controller.support_polygon_problem, state, controller.lowlevel.contacts)

    # State transitions

    # Linear momentum control
    c = center_of_mass(state)
    h = momentum(state)
    ċ = FreeVector3D(h.frame, linear(h) / controller.robotmass)
    z_des = controller.comref.v[3] # TODO
    l̇_des = controller.linear_momentum_controller(c, ċ, z_des, controller.support_polygon_problem; ξ_des=controller.comref) # TODO: remove ξ_des
    centroidal = centroidal_frame(controller.lowlevel)
    world_to_centroidal = Transform3D(c.frame, centroidal, -c.v)
    l̇_des = transform(l̇_des, world_to_centroidal)
    setdesired!(controller.linmomtask, l̇_des)

    # Pelvis orientation control
    pelvistask = controller.pelvistask
    pelvisgains = controller.pelvisgains
    pelvis = target(pelvistask.path)
    Hpelvis = transform_to_root(state, pelvis)
    Tpelvis = transform(twist_wrt_world(state, pelvis), inv(Hpelvis))
    ωd_pelvis_des = FreeVector3D(Tpelvis.frame, pd(pelvisgains, rotation(Hpelvis), Tpelvis.angular))
    setdesired!(pelvistask, ωd_pelvis_des)

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

function set_contact_points!(support_polygon_problem::ConvexHullProblem, state::MechanismState, contacts)
    i = 1
    for (body, contactpoints) in contacts
        to_world = transform_to_root(state, body)
        for contactpoint in contactpoints
            position = transform(contactpoint.position::Point3D{SVector{3, Float64}}, to_world)
            @inbounds set_vertex!(support_polygon_problem, i, horizontal_projection(position.v))
            i += 1
        end
    end
    nothing
end
