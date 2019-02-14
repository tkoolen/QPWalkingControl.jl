module ICPTrajectoryTest

using PushRecovery
using StaticArrays
using MathOptInterface
using OSQP
using OSQP.MathOptInterfaceOSQP
using Test
using Random
using ForwardDiff
using LinearAlgebra
using PlanarConvexHulls
using QPControl.Trajectories

const MOI = MathOptInterface

using PushRecovery: integrate_icp
import PushRecovery: ICPTrajectoryGenerator, push_segment!, solve!, initial_icp, cop, find_segment, initialize!
import PushRecovery: HRep, SHRep

const visualize = false
@static if visualize
    import Plots
    using Plots: plot, plot!, gui
    Plots.gr()
end

@testset "integrate_icp" begin
    ξ0 = 1.0
    p = BezierCurve(1, 2, 3, 4)
    ω = 3.0
    for t in range(0.0, 1.0, length=10)
        ξ = integrate_icp(ξ0, p, ω, t)
        ξd = ForwardDiff.derivative(t -> integrate_icp(ξ0, p, ω, t), t)
        @test ξd ≈ ω * (ξ - p(t))
    end
end

function optimizer()
    optimizer = OSQP.Optimizer()
    MOI.set(optimizer, OSQPSettings.Verbose(), false)
    MOI.set(optimizer, OSQPSettings.EpsAbs(), 1e-6)
    MOI.set(optimizer, OSQPSettings.EpsRel(), 1e-8)
    MOI.set(optimizer, OSQPSettings.MaxIter(), 20000)
    # MOI.set(optimizer, OSQPSettings.Polish(), true)
    MOI.set(optimizer, OSQPSettings.AdaptiveRhoInterval(), 25) # required for deterministic behavior
    optimizer
end

@testset "straight line walking" begin
    visualize = true

    g = 9.81
    z = 0.95
    ω = sqrt(g / z)
    num_segments = 8
    generator = ICPTrajectoryGenerator{Float64, 6}(optimizer(), num_segments, ω)

    foot_width = 0.15
    foot_length = 0.3
    step_length = 0.5
    step_width = 0.2
    foot_points = SVector(
        SVector(-foot_length / 2, -foot_width / 2),
        SVector(-foot_length / 2,  foot_width / 2),
        SVector( foot_length / 2,  foot_width / 2),
        SVector( foot_length / 2, -foot_width / 2)
    )

    foot_locations = [foot_points .+ Ref(SVector(0.0, -step_width)), foot_points .+ Ref(SVector(0.0, step_width))]
    convex_hulls = [DConvexHull{Float64}() for i = 1 : num_segments]
    Δts = fill(NaN, num_segments)
    support_indices = Set([1, 2])
    next_swing_index = 1
    first_step = true
    rng = MersenneTwister(1)
    for i = 1 : num_segments
        convex_hull = convex_hulls[i]
        support_points = collect(Iterators.flatten(foot_locations[i] for i in support_indices))
        jarvis_march!(convex_hull, support_points)
        if length(support_indices) == length(foot_locations)
            Δts[i] = first_step ? 1.0 + 0.2 * rand(rng) : 0.2 + 0.2 * rand(rng)
            delete!(support_indices, next_swing_index)
        else
            Δts[i] = 0.6 + 0.2 * rand(rng)
            for i in eachindex(foot_locations)
                push!(support_indices, i)
            end
            foot_displacement = if first_step
                first_step = false
                SVector(step_length, 0.0)
            else
                SVector(2 * step_length, 0.0)
            end
            foot_locations[next_swing_index] = foot_locations[next_swing_index] .+ Ref(foot_displacement)
            next_swing_index = next_swing_index == 1 ? 2 : 1
        end
    end
    support_polygon = Piecewise(Constant.(convex_hulls), cumsum([0.0; Δts]))

    for num_active_segments in 1 : num_segments
        tf = sum(Δts[1 : num_active_segments]) - 10 * eps()
        empty!(generator)
        p0 = SVector(0.02, 0.01)
        ξ0 = p0
        ξf = centroid(convex_hulls[num_active_segments])
        initialize!(generator, p0, ξ0, ξf)
        for i = 1 : num_active_segments
            push_segment!(generator, Δts[i], convex_hulls[i], centroid(convex_hulls[i]))
        end

        solve!(generator)

        # Initial/final ICP
        @test first(generator(0.0)) ≈ generator.initial_icp[] atol=1e-4
        @test first(generator(tf)) ≈ generator.final_icp[] atol=1e-4

        # CoP in support polygon
        for t in range(0.0, tf; length=100)
            @test cop(generator, t) in support_polygon(t)
        end

        # Continuity
        let t = 0.0
            for i in 1 : num_active_segments
                for t′ in (max(0.0, t - eps()), t, t + eps())
                    ξ, ξd = generator(t′)
                    @test ξ ≈ initial_icp(generator, i) atol=1e-4
                    i′, t0′ = find_segment(generator, t′)
                    if t′ > t
                        @test i′ == i
                        @test t0′ == t
                    elseif t′ < t
                        @test i′ == i - 1
                    end
                end
                t += generator.Δts[i]
            end

            for t′ in (t, t + eps())
                ξ, ξd = generator(t′)
                @test ξ ≈ generator.final_icp[] atol=1e-4
                @test ξd ≈ SVector(0.0, 0.0) atol=1e-4
            end
        end

        allocs = @allocated solve!(generator)
        @test allocs == 0

        @static if visualize
            ts = range(0.0, sum(Δts), length=200)
            xlim = [-0.5, 2.5]
            ylim = [-0.5, 0.5]
            width = 1600
            height = round(Int, width * diff(ylim)[1] / diff(xlim)[1])
            plt = plot(xlabel = "x", ylabel = "y", xlim = xlim, ylim = ylim, size = (width, height))
            cops = cop.(Ref(generator), ts)
            plot!(plt, getindex.(cops, 1), getindex.(cops, 2), label = "cop")
            icps = first.(generator.(ts))
            plot!(plt, getindex.(icps, 1), getindex.(icps, 2), label = "icp")
            for i in 1 : num_active_segments
                hull = convex_hulls[i]
                points = push!(copy(vertices(hull)), first(vertices(hull)))
                plot!(plt, getindex.(points, 1), getindex.(points, 2), label = "hull $i")
            end
            gui(plt)
            sleep(1)
        end
    end
end

end # module
