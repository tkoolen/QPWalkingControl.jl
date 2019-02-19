struct ICPWalkingStateMachine{N, I<:ICPTrajectoryGenerator, S<:SE3PDController}
    contacts::Dict{BodyID, Vector{QPControl.ContactPoint{N}}}
    icp_trajectory_generator::I
    end_effector_controllers::Dict{BodyID, S}
    support_polygon_plan::Piecewise{Constant{DConvexHull{Float64}}, Float64, Vector{Constant{DConvexHull{Float64}}}, Vector{Float64}}
    bodies_in_contact::Vector{BodyID}
    pointsbuffer::Vector{SVector{2, Float64}}
    # TODO: footstep generator
end

function ICPWalkingStateMachine(
        contacts::Dict{BodyID, <:Vector{<:QPControl.ContactPoint}},
        icp_trajectory_generator::ICPTrajectoryGenerator,
        end_effector_controllers::Dict{BodyID, <:SE3PDController},
    )
    T = Float64
    n = num_segments(icp_trajectory_generator)
    convex_hulls = [DConvexHull{Float64}() for i = 1 : n]
    times = zeros(T, n + 1)
    support_polygon_plan = Piecewise(Constant.(convex_hulls), times, clamp=true)
    bodies_in_contact = sort!(collect(keys(contacts)), by = x -> x.value)
    pointsbuffer = SVector{2, T}[]
    ICPWalkingStateMachine(
        contacts, icp_trajectory_generator, end_effector_controllers,
        support_polygon_plan, bodies_in_contact, pointsbuffer)
end


function (statemachine::ICPWalkingStateMachine)(t, state::MechanismState)

end

