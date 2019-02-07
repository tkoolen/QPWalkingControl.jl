struct ICPController{T, I, Z}
    kxy::T
    zgains::PDGains{T, T}
    m::T
    gz::T
    icptraj::I
    ztraj::Z
    contactmode::PlanarContactMode{T}
end

function ICPController(mechanism::Mechanism{T}, icptraj, z_des::Number) where T
    ztraj = t -> (z_des, zero(z_des), zero(z_des))
    ICPController(
        T(3.0),
        critically_damped_gains(10.0),
        mass(mechanism),
        norm(mechanism.gravitational_acceleration),
        icptraj,
        ztraj,
        PlanarContactMode{T}(root_frame(mechanism))
    )
end

function (controller::ICPController)(t, c::Point3D, cd::FreeVector3D, active_contact_points, state)
    update!(controller.contactmode, active_contact_points, state)
    support_polygon = controller.contactmode.hull

    # Equation (27) in "Design of a momentum-based control framework and application to the humanoid robot atlas"
    # minus the integral term.
    kxy = controller.kxy
    zgains = controller.zgains
    m = controller.m
    gz = controller.gz
    ξ_des_2d, ξd_des_2d = controller.icptraj(t)
    ξ_des = Point3D(c.frame, ξ_des_2d[1], ξ_des_2d[2], zero(eltype(ξ_des_2d))) # TODO: kind of nasty
    ξd_des = FreeVector3D(c.frame, ξd_des_2d[1], ξd_des_2d[2], zero(eltype(ξd_des_2d))) # TODO: kind of nasty
    z_des, zd_des, zdd_des = controller.ztraj(t)
    z = c.v[3]
    zd = cd.v[3]
    ω = sqrt(abs(gz / z))

    ξ = icp(c, cd, ω)
    cmp_des = ξ - ξd_des / ω + kxy * (ξ - ξ_des)

    cmp_des_projected_xy = closest_point(horizontal_projection(cmp_des.v), support_polygon)
    cmp_des_projected = Point3D(cmp_des.frame, cmp_des_projected_xy[1], cmp_des_projected_xy[2], cmp_des.v[3])

    ld_des_xy = (c - cmp_des_projected) * (m * gz) / z
    ld_des_z = m * (pd(zgains, z, z_des, zd, zd_des) + zdd_des)
    FreeVector3D(c.frame, ld_des_xy.v[1], ld_des_xy.v[2], ld_des_z)
end
