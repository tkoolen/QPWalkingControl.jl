using LightXML

struct JointDamping{T}
    coeffs::Dict{JointID, Vector{T}}
end

function JointDamping{T}(mechanism::Mechanism, urdf::AbstractString) where T
    xdoc = parse_file(urdf)
    xroot = LightXML.root(xdoc)
    @assert LightXML.name(xroot) == "robot"
    coeffs = Vector{Pair{JointID, Vector{T}}}()
    xml_joints = get_elements_by_tagname(xroot, "joint")
    for xml_joint in xml_joints
        name = attribute(xml_joint, "name")
        try
            joint = findjoint(mechanism, name)
            xml_dynamics = find_element(xml_joint, "dynamics")
            if xml_dynamics != nothing
                xml_damping = attribute(xml_dynamics, "damping")
                if xml_damping != nothing
                    push!(coeffs, JointID(joint) => [parse(T, xml_damping)])
                end
            end
        catch # TODO: be more specific
        end
    end
    return JointDamping(Dict(coeffs))
end

function (damping::JointDamping)(τ::AbstractVector, t, state::MechanismState)
    τ .= 0
    v = velocity(state)
    @inbounds for (jointid, coeffs) in damping.coeffs
        τjoint = τ[jointid]
        vjoint = v[jointid]
        τjoint .= .-coeffs .* vjoint
    end
    return τ
end
