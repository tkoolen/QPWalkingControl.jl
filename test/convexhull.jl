module ConvexHullTest

using PushRecovery.ConvexHull
using StaticArrays
using MathOptInterface
using OSQP
using OSQP.MathOptInterfaceOSQP
using Test
using Random
using LinearAlgebra

const MOI = MathOptInterface

function optimizer(atol_distance)
    optimizer = OSQP.Optimizer()
    MOI.set(optimizer, OSQPSettings.Verbose(), false)
    MOI.set(optimizer, OSQPSettings.EpsAbs(), atol_distance^2)
    MOI.set(optimizer, OSQPSettings.EpsRel(), 1e-8)
#     MOI.set(optimizer, OSQPSettings.Polish(), true)
#     MOI.set(optimizer, OSQPSettings.AdaptiveRhoInterval(), 25) # required for deterministic behavior
    optimizer
end

@testset "ConvexHull basics" begin
    atol_distance = 1e-3
    problem = ConvexHullProblem{2, 3, Float64}(optimizer(atol_distance))
    set_vertices!(problem, [SVector(0., 0.), SVector(2., 0.), SVector(0., 3.)])

    p = SVector(0.5, 0.5)
    set_point!(problem, p)
    solve!(problem)
    @test is_point_inside(problem; atol=atol_distance)
    @test distance_to_closest_point(problem) ≈ 0 atol=atol_distance
    @test closest_point(problem) ≈ p atol=atol_distance

    p = SVector(-0.5, -0.5)
    set_point!(problem, p)
    solve!(problem)
    @test !is_point_inside(problem; atol=atol_distance)
    @test distance_to_closest_point(problem) ≈ sqrt(0.5) atol=atol_distance
    @test closest_point(problem) ≈ SVector(0.0, 0.0) atol=atol_distance
end

@testset "ConvexHull point inside" begin
    atol_distance = 1e-3
    N = 2
    M = 10
    problem = ConvexHullProblem{N, M, Float64}(optimizer(atol_distance))
    Random.seed!(1)
    for i = 1 : 10000
        vertices = [rand(SVector{N}) for i = 1 : M]
        set_vertices!(problem, vertices)
        weights = rand(M)
        weights ./= sum(weights)
        point = sum(weights .* vertices)
        set_point!(problem, point)
        solve!(problem)
        @test is_point_inside(problem; atol=atol_distance)
    end
end

@testset "ConvexHull point outside" begin
    atol_distance = 1e-3
    N = 2
    M = 10
    problem = ConvexHullProblem{N, M, Float64}(optimizer(atol_distance))
    Random.seed!(1)
    for i = 1 : 10000
        set_vertices!(problem, [(1 - atol_distance) * normalize(randn(SVector{N})) for i = 1 : M])
        point = normalize(randn(SVector{N}))
        set_point!(problem, point)
        solve!(problem)
        @test !is_point_inside(problem; atol=atol_distance)
        closest = closest_point(problem)
        @test norm(point - closest) ≈ distance_to_closest_point(problem) atol=atol_distance

        set_point!(problem, closest)
        solve!(problem)
        @test is_point_inside(problem; atol=atol_distance)
        @test distance_to_closest_point(problem) ≈ 0 atol=atol_distance
    end
end

@testset "ConvexHull allocations" begin
    atol_distance = 1e-3
    N = 2
    M = 10
    problem = ConvexHullProblem{N, M, Float64}(optimizer(atol_distance))
    vertices = [rand(SVector{N}) for i = 1 : M]
    point = rand(SVector{N})

    testfun = function (problem, vertices, point)
        set_vertices!(problem, vertices)
        set_point!(problem, point)
        solve!(problem)
    end

    testfun(problem, vertices, point)
    allocs = @allocated testfun(problem, vertices, point)
    @test allocs == 0
end

end
