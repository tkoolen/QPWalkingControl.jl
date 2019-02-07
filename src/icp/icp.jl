function icp(c::Point3D, ċ::FreeVector3D, ω::Number)
    @framecheck c.frame ċ.frame
    c + ċ / ω
end

function icp(state::MechanismState, ω::Number)
    c = center_of_mass(state)
    h = momentum(state)
    ċ = FreeVector3D(h.frame, linear(momentum(state)) / mass(state.mechanism))
    icp(c, ċ, ω)
end
