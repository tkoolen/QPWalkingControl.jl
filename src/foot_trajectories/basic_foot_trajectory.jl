struct BasicFootTrajectory{T}
    p0::SVector{3, T}
    xaxis::SVector{2, T}
    Δx::Polynomial{4, T}
    Δz::Polynomial{3, T}
    s::Polynomial{6, T}
end

function BasicFootTrajectory(t0::Number, tf::Number, p0::SVector{3}, Δzmid::Number, pf::SVector{3}, zdf::Number)
    # Coordinates of final position relative to initial, projected onto swing plane (x-z)
    Δpf = pf - p0
    Δxf = sqrt(Δpf[1]^2 + Δpf[2]^2)
    T = eltype(Δxf)
    Δzf = convert(T, Δpf[3])
    xaxis = SVector(Δpf[1], Δpf[2]) / Δxf

    # Trajectories as function of interpolation variable `s`
    s0 = zero(T)
    sf = oneunit(T)
    smid = (s0 + sf) / 2
    Δx = fit_cubic(x0=s0, xf=sf, y0=zero(T), yd0=zero(T), yf=Δxf, ydf=zero(T))
    Δz = let
        num_coeffs = Val(3)
        A = hcat(
            SVector(SUP.coefficient_gradient(s0, num_coeffs, Val(0))),
            SVector(SUP.coefficient_gradient(smid, num_coeffs, Val(0))),
            SVector(SUP.coefficient_gradient(sf, num_coeffs, Val(0))),
        ) |> transpose
        b = SVector(tuple(zero(T), max(0, Δzf) + Δzmid, Δzf))
        Polynomial(Tuple(A \ b))
    end

    # Final value of derivative of height as a function of `s`
    dzds_f = derivative(Δz)(sf)

    # Final value of derivative of s as a function of t
    sdf = zdf / dzds_f

    # Interpolation trajectory
    s = fit_quintic(x0=t0, xf=tf, y0=s0, yd0=zero(s0), ydd0=zero(s0), yf=sf, ydf=sdf, yddf=zero(sf))

    BasicFootTrajectory(p0, xaxis, Δx, Δz, s)
end

function (traj::BasicFootTrajectory)(t::Number, ::Val{2})
    s = traj.s(t)
    sd = derivative(traj.s)(t)
    sdd = derivative(derivative(traj.s))(t)
    Δx = traj.Δx(s)
    Δxd = derivative(traj.Δx)(s) * sd
    Δxdd = derivative(traj.Δx)(s) * sdd + derivative(derivative(traj.Δx))(s) * sd
    Δz = traj.Δz(s)
    Δzd = derivative(traj.Δz)(s) * sd
    Δzdd = derivative(traj.Δz)(s) * sdd + derivative(derivative(traj.Δz))(s) * sd
    p = traj.p0 + SVector(Δx * traj.xaxis[1], Δx * traj.xaxis[2], Δz)
    pd = SVector(Δxd * traj.xaxis[1], Δxd * traj.xaxis[2], Δzd)
    pdd = SVector(Δxdd * traj.xaxis[1], Δxdd * traj.xaxis[2], Δzdd)
    p, pd, pdd
end

(traj::BasicFootTrajectory)(t::Number) = first(traj(t, Val(2)))
