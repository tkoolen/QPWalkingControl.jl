# TODO: consider taking contact information into account
struct PDCoMController{T, G<:PDGains, C}
    gains::G
    com_trajectory::C
    mass::T
end

function (controller::PDCoMController)(t, c::Point3D, cd::FreeVector3D, active_contact_points, state)
    mass = controller.mass
    # c_des, cd_des, cdd_des = controller.com_trajectory(t, Val(2))
    # ld_des = mass * (pd(controller.gains, c, c_des, cd, cd_des) + cdd_des)
    c_des = Point3D(c.frame, controller.com_trajectory(t))
    ld_des = mass * (pd(controller.gains, c, c_des, cd, zero(cd)))
    return ld_des
end
