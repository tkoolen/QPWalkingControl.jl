module BezierCurves

export
    BezierCurve,
    derivative

import StaticUnivariatePolynomials

using StaticUnivariatePolynomials: derivative
using Base: tail

struct BezierCurve{N, T}
    points::NTuple{N, T}
end

BezierCurve(coeffs::Tuple) = BezierCurve(promote(coeffs...))
BezierCurve(coeffs...) = BezierCurve(coeffs)

@inline (b::BezierCurve{1})(t) = b.points[1]
@inline function (b::BezierCurve)(t)
    b1 = BezierCurve(reverse(tail(reverse(b.points))))
    b2 = BezierCurve(tail(b.points))
    (oneunit(t) - t) * b1(t) + t * b2(t)
end

@inline function StaticUnivariatePolynomials.derivative(b::BezierCurve{N}) where {N}
    # https://pages.mtu.edu/~shene/COURSES/cs3621/NOTES/spline/Bezier/bezier-der.html
    points = b.points
    n = N - 1
    BezierCurve(ntuple(i -> n * (points[i + 1] - points[i]), n))
end

end # module
