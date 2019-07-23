const QUINTIC_ZERO_ONE_INTERPOLATOR = fit_quintic(x0=0.0, xf=1.0, y0=0.0, yd0=0.0, ydd0=0.0, yf=1.0, ydf=0.0, yddf=0.0)

function interpolated_orientation_trajectory(t0::Number, tf::Number, rot0::Rotation, rotf::Rotation)
    Interpolated(t0, tf, Quat(rot0), Quat(rotf), QUINTIC_ZERO_ONE_INTERPOLATOR, min_num_derivs=Val(2))
end

function init_support!(end_effector_controller::SE3PDController; t0::Number, tf::Number)
    # TODO: frame lookup is kind of nasty
    bodyframe = end_effector_controller.trajectory[].body
    baseframe = end_effector_controller.trajectory[].base
    end_effector_controller.gains[] = SE3PDGains(FramePDGains(bodyframe, PDGains(0.0, 20.0)), FramePDGains(bodyframe, PDGains(0.0, 20.0)))
    angulartraj = interpolated_orientation_trajectory(t0, tf, one(Quat{Float64}), one(Quat{Float64}))
    lineartraj = convert(BasicFootTrajectory{Float64}, SVector(0.0, 0.0, 0.0), tf)
    end_effector_controller.trajectory[] = SE3Trajectory(bodyframe, baseframe, angulartraj, lineartraj)
end

function init_swing!(end_effector_controller::SE3PDController, pose0::Transform3D, posef::Transform3D;
        t0::Number, tf::Number, Δzmid::Number=0.1, zdf::Number=-0.1)
    # TODO: frame lookup is kind of nasty
    bodyframe = end_effector_controller.trajectory[].body
    baseframe = end_effector_controller.trajectory[].base
    end_effector_controller.gains[] = SE3PDGains(
        FramePDGains(bodyframe, critically_damped_gains(100.)),
        FramePDGains(bodyframe, critically_damped_gains(100.))
    )
    angulartraj = interpolated_orientation_trajectory(t0, tf, rotation(pose0), rotation(posef))
    lineartraj = BasicFootTrajectory(t0, tf, translation(pose0), Δzmid, translation(posef), zdf)
    end_effector_controller.trajectory[] = SE3Trajectory(bodyframe, baseframe, angulartraj, lineartraj)
end
