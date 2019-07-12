# TODO: move to QPControl
struct PosePlan{T<:Number}
    move_start_times::Vector{T}
    move_durations::Vector{T}
    final_poses::Vector{Transform3D{T}}
end

PosePlan{T}() where {T} = PosePlan(T[], T[], Transform3D{T}[])

function Base.push!(plan::PosePlan, start_time::Number, duration::Number, final_pose::Transform3D)
    push!(plan.move_start_times, start_time)
    push!(plan.move_durations, duration)
    push!(plan.final_poses, final_pose)
    return plan
end
