module BezierTest

using PushRecovery.BezierCurves
using Test
using ForwardDiff

@testset "evaluation" begin
    b = BezierCurve(1, 2, 3, 4)
    @test b(0) === 1
    @test b(1) === 4
end

@testset "derivative" begin
    b = BezierCurve(1, 2, 3, 4)
    b′ = derivative(b)
    for t in range(0., 1., length=10)
        @test ForwardDiff.derivative(b, t) ≈ b′(t) atol=1e-10
    end
end

end
