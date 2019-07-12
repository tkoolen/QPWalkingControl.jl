function init_support!(end_effector_controller::SE3PDController, tf::Number)
    # TODO: frame lookup is kind of nasty
    bodyframe = end_effector_controller.trajectory[].body
    baseframe = end_effector_controller.trajectory[].base
    end_effector_controller.gains[] = SE3PDGains(FramePDGains(bodyframe, PDGains(0.0, 15.0)), FramePDGains(bodyframe, PDGains(0.0, 0.0)))
    angulartraj = Constant(one(Quat)) # TODO
    lineartraj = convert(BasicFootTrajectory{Float64}, SVector(0.0, 0.0, 0.0), tf)
    end_effector_controller.trajectory[] = SE3Trajectory(bodyframe, baseframe, angulartraj, lineartraj)
end

function init_swing!(end_effector_controller::SE3PDController, pose0::Transform3D, posef::Transform3D,
        t0::Number, tf::Number,
        Δzmid::Number=0.15, zdf::Number=-0.1) # TODO: kwargs
    # TODO: frame lookup is kind of nasty
    bodyframe = end_effector_controller.trajectory[].body
    baseframe = end_effector_controller.trajectory[].base
    end_effector_controller.gains[] = SE3PDGains(
        FramePDGains(bodyframe, critically_damped_gains(100.)),
        FramePDGains(bodyframe, critically_damped_gains(100.))
    )
    angulartraj = Constant(one(Quat)) # TODO
    lineartraj = BasicFootTrajectory(t0, tf, translation(pose0), Δzmid, translation(posef), zdf)
    end_effector_controller.trajectory[] = SE3Trajectory(bodyframe, baseframe, angulartraj, lineartraj)
end
