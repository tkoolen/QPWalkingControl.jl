function horizontal_projection(v::StaticVector{3})
    @inbounds return SVector(v[1], v[2])
end

function horizontal_projection(r::RotMatrix{3})
    @inbounds return RotMatrix(r[SVector(1, 2), SVector(1, 2)])
end

critically_damped_gains(k::Number) = PDGains(k, 2 * sqrt(k))

"""
``A x \\le b``
"""
struct HRep{T, TA<:AbstractMatrix{T}, Tb<:AbstractVector{T}}
    A::TA
    b::Tb
end

function HRep(hull::ConvexHull)
    HRep(hrep(hull)...)
end

function Base.:*(M::StaticMatrix, hrep::HRep)
    HRep(M * hrep.A, M * hrep.b)
end

function Base.:+(hrep::HRep, c::StaticVector)
    HRep(hrep.A, hrep.b + hrep.A * c)
end

function Base.in(x::AbstractVector, hrep::HRep)
    all(hrep.A * x .<= hrep.b)
end

function RigidBodyDynamics.transform(fhull::FrameAnnotated{C}, tf::Transform3D) where C<:ConvexHull
    @framecheck fhull.frame tf.from
    R = horizontal_projection(rotation(tf))
    p = horizontal_projection(translation(tf))
    return in_frame(tf.to, C(map(x -> R * x + p, vertices(unwrap(fhull)))))
end

const SHRep{N, M, T, L} = HRep{T, SMatrix{N, M, T, L}, SVector{N, T}}
const MHRep{N, M, T, L} = HRep{T, MMatrix{N, M, T, L}, MVector{N, T}}

Base.zero(::Type{SHRep{N, M, T, L}}) where {N, M, T, L} = HRep(zero(SMatrix{N, M, T, L}), zero(SVector{N, T}))
Base.zero(::Type{MHRep{N, M, T, L}}) where {N, M, T, L} = HRep(zero(MMatrix{N, M, T, L}), zero(MVector{N, T}))
