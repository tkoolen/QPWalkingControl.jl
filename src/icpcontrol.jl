struct ICPController{T}
    kxy::T
    zgains::PDGains{T, T}
    m::T
    gz::T
    shrinkfactor::T
end

function ICPController(mechanism::Mechanism{T}) where T
    ICPController(T(3.0), PDGains(T(10.), T(2 * sqrt(10.0))), mass(mechanism), norm(mechanism.gravitational_acceleration), 0.01)
end

function (controller::ICPController{T})(c::Point3D, ċ::FreeVector3D, z_des::Number, convexhullproblem::ConvexHullProblem{2};
        zd_des=zero(T), zdd_des=zero(T), ξ_des=zero(c), ξd_des=zero(ċ)) where T
    # Equation (27) in "Design of a momentum-based control framework and application to the humanoid robot atlas"
    # minus the integral term.
    kxy = controller.kxy
    zgains = controller.zgains
    m = controller.m
    gz = controller.gz
    z = c.v[3]
    zd = ċ.v[3]
    ω = sqrt(abs(gz / z))

    ξ = c + ċ / ω
    cmp_des = ξ - ξd_des / ω + kxy * (ξ - ξ_des)

    set_point!(convexhullproblem, horizontal_projection(cmp_des.v))
    set_shrink_factor!(convexhullproblem, controller.shrinkfactor)
    solve!(convexhullproblem)
    cmp_des_projected = Point3D(cmp_des.frame, closest_point(convexhullproblem)..., cmp_des.v[3])

    ld_des_xy = (c - cmp_des_projected) * (m * gz) / z
    ld_des_z = m * (pd(zgains, z, z_des, zd, zd_des) + zdd_des)
    FreeVector3D(c.frame, ld_des_xy.v[1], ld_des_xy.v[2], ld_des_z)
end
