# TODO: consider taking contact information into account
struct PDCoMController{T, G<:PDGains, C}
    gains::G
    com_trajectory::C
    mass::T
end

function (controller::PDCoMController)(t, c::Point3D, cd::FreeVector3D, active_contact_points, state)
    @framecheck c.frame cd.frame
    mass = controller.mass
    c_des, cd_des, cdd_des = controller.com_trajectory(t, Val(2))
    ld_des = FreeVector3D(c.frame, mass * (pd(controller.gains, c.v, c_des, cd.v, cd_des) + cdd_des.v))
    return ld_des
end
