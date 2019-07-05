using LightXML
using RigidBodyDynamics: JointDict
using RigidBodyDynamics.CustomCollections: SegmentedVector

struct JointDamping{T}
    coeffs::SegmentedVector{JointID, T, Base.OneTo{JointID}, Vector{T}}
end

function JointDamping{T}(mechanism::Mechanism, urdf::AbstractString) where T
    xdoc = parse_file(urdf)
    xroot = LightXML.root(xdoc)
    @assert LightXML.name(xroot) == "robot"
    coeffs = SegmentedVector(zeros(num_velocities(mechanism)), tree_joints(mechanism), num_velocities)
    xml_joints = get_elements_by_tagname(xroot, "joint")
    for xml_joint in xml_joints
        name = attribute(xml_joint, "name")
        xml_dynamics = find_element(xml_joint, "dynamics")
        if xml_dynamics != nothing
            xml_damping = attribute(xml_dynamics, "damping")
            if xml_damping != nothing
                local joint
                try
                    joint = findjoint(mechanism, name)
                catch
                    @warn "Joint with name \"$name\" not found. Ignoring damping coefficient in URDF."
                    continue
                end
                coeffs[JointID(joint)] .= parse(T, xml_damping)
            end
        end
    end
    return JointDamping(coeffs)
end

function (damping::JointDamping)(τ::SegmentedVector{JointID}, t, state::MechanismState)
    τ .= .-damping.coeffs .* velocity(state)
    return τ
end
