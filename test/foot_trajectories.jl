module FootTrajectoryTest

using Test
using StaticArrays
import QPWalkingControl
using QPWalkingControl: BasicFootTrajectory

@testset "BasicFootTrajectory" begin
    t0 = 0.2
    tf = 1.5
    p0 = SVector(0.0, 0.0, 0.0)
    Δzmid = 0.1
    pf = SVector(0.8, 0.0, 0.05)
    zdf = -0.1
    traj = BasicFootTrajectory(t0, tf, p0, Δzmid, pf, zdf)

    @test traj(t0) ≈ p0 atol=1e-10
    @test traj(tf) ≈ pf atol=1e-10
    @test traj(tf, Val(2))[2][3] ≈ zdf atol=1e-10

    zs = map(t -> traj(t)[3], range(t0, tf, length=1000))
    @test any(z -> z >= max(p0[3], pf[3]) + Δzmid, zs)
    @test all(z -> z >= min(p0[3], pf[3]), zs)

    allocs = @allocated BasicFootTrajectory(t0, tf, p0, Δzmid, pf, zdf)
    @test allocs == 0
end

end
