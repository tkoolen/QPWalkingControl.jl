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
struct SHRep{N, M, T, L}
    A::SMatrix{N, M, T, L}
    b::SVector{N, T}
end

function SHRep(hull::ConvexHull)
    A, b = hrep(hull)
    SHRep(A, b)
end

function Base.zero(::Type{SHRep{N, M, T, L}}) where {N, M, T, L}
    SHRep(zero(SMatrix{N, M, T, L}), zero(SVector{N, T}))
end

function Base.:*(M::StaticMatrix, hrep::SHRep)
    SHRep(M * hrep.A, M * hrep.b)
end

function Base.:+(hrep::SHRep, c::StaticVector)
    SHRep(hrep.A, hrep.b + hrep.A * c)
end

function Base.in(x::AbstractVector, hrep::SHRep)
    all(hrep.A * x .<= hrep.b)
end
