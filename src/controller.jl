struct PushRecoveryController{M<:MomentumBasedController, C<:ConvexHullProblem}
    lowlevel::M
    robotmass::Float64

    foottasks::Dict{RigidBody{Float64}, SpatialAccelerationTask}
    linmomtask::LinearMomentumRateTask
    pelvistask::AngularAccelerationTask
    jointtasks::Dict{JointID, JointAccelerationTask{Revolute{Float64}}}

    icpcontroller::ICPController{Float64}
    pelvisgains::PDGains{Float64,Float64}
    jointgains::Dict{JointID, PDGains{Float64, Float64}}

    comref::Point3D{SVector{3, Float64}}
    jointrefs::Dict{JointID, Float64}

    convexhullproblem::C
end

function PushRecoveryController(
        lowlevel::MomentumBasedController,
        feet::Vector{<:RigidBody},
        pelvis::RigidBody,
        nominalstate::MechanismState,
        convexhulloptimizer::MOI.AbstractOptimizer;
        joint_regularization::Float64 = 0.05,
        linear_momentum_weight::Float64 = 1.0,
        icpcontroller::ICPController = ICPController(lowlevel.state.mechanism),
        pelvisgains::PDGains = PDGains(20., 2 * sqrt(20.0)),
        jointgains = Dict(JointID(j) => PDGains(100.0, 20.) for j in tree_joints(lowlevel.state.mechanism)),
        comref::Point3D = center_of_mass(nominalstate) - FreeVector3D(root_frame(lowlevel.state.mechanism), 0., 0., 0.05),
        jointrefs = Dict(JointID(j) => configuration(nominalstate, j)[1] for j in tree_joints(lowlevel.state.mechanism) if joint_type(j) isa Revolute))
    mechanism = lowlevel.state.mechanism
    world = root_body(mechanism)
    worldframe = root_frame(mechanism)
    m = mass(mechanism)

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
    convexhullproblem = ConvexHullProblem{2, num_contacts, Float64}(convexhulloptimizer)

    PushRecoveryController(
        lowlevel, m,
        foottasks, linmomtask, pelvistask, jointtasks,
        icpcontroller, pelvisgains, jointgains,
        comref, jointrefs,
        convexhullproblem)
end

function (controller::PushRecoveryController)(τ::AbstractVector, t::Number, state::MechanismState)
    # Support polygon setup
    convexhullproblem = controller.convexhullproblem
    update_support_polygon_problem!(convexhullproblem, state, controller.lowlevel.contacts)

    # Linear momentum control
    m = controller.robotmass
    c = center_of_mass(state)
    h = momentum(state)
    ċ = FreeVector3D(h.frame, linear(h) / m)
    zdesired = controller.comref.v[3]
    l̇desired = controller.icpcontroller(c, ċ, zdesired, convexhullproblem; ξ_des=controller.comref)
    centroidal = centroidal_frame(controller.lowlevel)
    world_to_centroidal = Transform3D(c.frame, centroidal, -c.v)
    l̇desired = transform(l̇desired, world_to_centroidal)
    setdesired!(controller.linmomtask, l̇desired)

    # Pelvis orientation control
    pelvistask = controller.pelvistask
    pelvisgains = controller.pelvisgains
    pelvis = target(pelvistask.path)
    Hpelvis = transform_to_root(state, pelvis)
    Tpelvis = transform(twist_wrt_world(state, pelvis), inv(Hpelvis))
    ωdesired = FreeVector3D(Tpelvis.frame, pd(pelvisgains, rotation(Hpelvis), Tpelvis.angular))
    setdesired!(pelvistask, ωdesired)

    # Joint position control
    for jointid in keys(controller.jointtasks)
        task = controller.jointtasks[jointid]
        gains = controller.jointgains[jointid]
        ref = controller.jointrefs[jointid]
        v̇desired = pd(gains, configuration(state, jointid)[1], ref, velocity(state, jointid)[1], 0.0)
        setdesired!(task, v̇desired)
    end

    controller.lowlevel(τ, t, state)
    τ
end

function update_support_polygon_problem!(convexhullproblem::ConvexHullProblem, state::MechanismState, contacts)
    i = 1
    for (body, contactpoints) in contacts
        to_world = transform_to_root(state, body)
        for contactpoint in contactpoints
            position = transform(contactpoint.position::Point3D{SVector{3, Float64}}, to_world)
            @inbounds set_vertex!(convexhullproblem, i, horizontal_projection(position.v))
            i += 1
        end
    end
    nothing
end
