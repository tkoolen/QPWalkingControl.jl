# TODO: consider taking contact information into account
struct PDCoMController{T, G<:PDGains, C}
    gains::G
    com_trajectory::C
    m::T
end

function (controller::PDCoMController)(t, c::Point3D, cd::FreeVector3D, active_contact_points, state)
    m = controller.m
    c_des, cd_des, cdd_des = controller.com_trajectory(t)
    ld_des = m * (pd(controller.gains, c, c_des, cd, cd_des) + cdd_des)
    return ld_des
end
