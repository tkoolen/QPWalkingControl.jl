function icp(c::Point3D, ċ::FreeVector3D, ω::Number)
    @framecheck c.frame ċ.frame
    Point3D(c.frame, horizontal_projection((c + ċ / ω).v)..., 0)
end

function icp(state::MechanismState, ω::Number)
    c = center_of_mass(state)
    h = momentum(state)
    ċ = FreeVector3D(h.frame, linear(momentum(state)) / mass(state.mechanism))
    icp(c, ċ, ω)
end

function icp(state::MechanismState)
    c = center_of_mass(state)
    z = c.v[3]
    gz = norm(state.mechanism.gravitational_acceleration)
    ω = sqrt(abs(gz / z))
    icp(state, ω)
end
