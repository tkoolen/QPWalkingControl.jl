struct ICPController{T}
    kxy::T
    zgains::PDGains{T, T}
    m::T
    gz::T
end

function ICPController(mechanism::Mechanism{T}) where T
    ICPController(T(3.0), critically_damped_gains(10.0), mass(mechanism), norm(mechanism.gravitational_acceleration))
end

function (controller::ICPController{T})(c::Point3D, cd::FreeVector3D, z_des::Number, support_polygon_problem::ConvexHullProblem{2};
        zd_des=zero(T), zdd_des=zero(T), ξ_des=zero(c), ξd_des=zero(cd)) where T
    # Equation (27) in "Design of a momentum-based control framework and application to the humanoid robot atlas"
    # minus the integral term.
    kxy = controller.kxy
    zgains = controller.zgains
    m = controller.m
    gz = controller.gz
    z = c.v[3]
    zd = cd.v[3]
    ω = sqrt(abs(gz / z))

    ξ = icp(c, cd, ω)
    cmp_des = ξ - ξd_des / ω + kxy * (ξ - ξ_des)

    set_point!(support_polygon_problem, horizontal_projection(cmp_des.v))
    solve!(support_polygon_problem)
    cmp_des_projected_xy = closest_point(support_polygon_problem)
    cmp_des_projected = Point3D(cmp_des.frame, cmp_des_projected_xy[1], cmp_des_projected_xy[2], cmp_des.v[3])

    ld_des_xy = (c - cmp_des_projected) * (m * gz) / z
    ld_des_z = m * (pd(zgains, z, z_des, zd, zd_des) + zdd_des)
    FreeVector3D(c.frame, ld_des_xy.v[1], ld_des_xy.v[2], ld_des_z)
end
