mutable struct PushApplier{T}
    path::TreePath{RigidBody{T}, Joint{T}}
    application_point::Point3D{SVector{3, T}}
    jacobian::PointJacobian{Matrix{T}}
    new_push::Bool
    t0::T
    Δt::T
    force::FreeVector3D{SVector{3, T}}
end

function PushApplier(mechanism::Mechanism{T}, application_point::Point3D;
        Δt::Number=zero(T), force::AbstractVector=zero(SVector{3, T})) where T
    world = root_body(mechanism)
    worldframe = root_frame(mechanism)
    frame = application_point.frame
    treepath = RigidBodyDynamics.path(mechanism, world, body_fixed_frame_to_body(mechanism, frame))
    jacobian = PointJacobian(worldframe, zeros(3, num_velocities(mechanism)))
    PushApplier(treepath, application_point, jacobian, false, typemax(T), T(Δt), FreeVector3D(worldframe, SVector{3, T}(force)))
end

function (controller::PushApplier)(τ, t, state::MechanismState)
    if controller.new_push
        controller.t0 = t
        controller.new_push = false
    end
    if controller.t0 <= t <= controller.t0 + controller.Δt
        jacobian = controller.jacobian
        point = transform(state, controller.application_point, jacobian.frame)
        point_jacobian!(jacobian, state, controller.path, point)
        mul!(τ, transpose(jacobian), controller.force)
    else
        τ .= 0
    end
    τ
end
