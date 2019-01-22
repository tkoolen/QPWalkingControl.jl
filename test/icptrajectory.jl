module ICPTrajectoryTest

using PushRecovery
using StaticArrays
using MathOptInterface
using OSQP
using OSQP.MathOptInterfaceOSQP
using Test
using Random
using LinearAlgebra

const MOI = MathOptInterface

import PushRecovery: ICPTrajectoryGenerator, solve!, initial_icp, cop
import PushRecovery: SHRep

using UnicodePlots

function optimizer()
    optimizer = OSQP.Optimizer()
    MOI.set(optimizer, OSQPSettings.Verbose(), false)
    MOI.set(optimizer, OSQPSettings.EpsAbs(), 1e-5)
    MOI.set(optimizer, OSQPSettings.EpsRel(), 1e-8)
    MOI.set(optimizer, OSQPSettings.MaxIter(), 20000)
    # MOI.set(optimizer, OSQPSettings.Polish(), true)
    # MOI.set(optimizer, OSQPSettings.AdaptiveRhoInterval(), 25) # required for deterministic behavior
    optimizer
end

@testset "straight line walking" begin
    g = 9.81
    z = 0.95
    ω = sqrt(g / z)
    num_segments = 4
    generator = ICPTrajectoryGenerator{Float64, 4}(optimizer(), ω, num_segments)

    foot_width = 0.15
    foot_length = 0.3
    step_length = 0.5
    step_width = 0.1
    Δt = 0.75
    foot_polygon = SHRep(
        @SMatrix([
             1.0  0.0;
            -1.0  0.0;
             0.0  1.0;
             0.0 -1.0
        ]),
        SVector(foot_width / 2, foot_width / 2, foot_length / 2, foot_length / 2)
    )
    foot_centers = map(1 : num_segments) do i
        SVector((i - 1) * step_length, (-1)^i * step_width)
    end

    for i = 1 : num_segments
        generator.Δts[i] = Δt
        generator.cop_polyhedra[i] = foot_polygon + foot_centers[i]
        generator.preferred_cops[i] = foot_centers[i]
    end
    generator.initial_icp[] = foot_centers[1] + SVector(0.02, 0.01)
    generator.final_icp[] = foot_centers[end]

    solve!(generator)

    @test initial_icp(generator, 1) ≈ generator.initial_icp[] atol=1e-4

    for i = 1 : num_segments
        C = generator.cop_polyhedra[i]
        @test cop(generator, i) ∈ C
    end

    allocs = @allocated solve!(generator)
    @test allocs == 0

    visualize = true
    if visualize
        cops = cop.(Ref(generator), 1 : num_segments)
        icps = initial_icp.(Ref(generator), 1 : num_segments)
        push!(icps, generator.final_icp[])
        plt = scatterplot(getindex.(cops, 1), getindex.(cops, 2), name = "cop", xlabel= "x", ylabel = "y")
        lineplot!(plt, getindex.(icps, 1), getindex.(icps, 2), name = "icp")
        display(plt)
    end
end

end # module
