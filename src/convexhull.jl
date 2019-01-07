module ConvexHull

export
    ConvexHullProblem,
    set_point!,
    set_vertex!,
    set_vertices!,
    solve!,
    is_point_inside,
    distance_to_closest_point,
    closest_point

using LinearAlgebra
using StaticArrays
using Parametron

import MathOptInterface

const MOI = MathOptInterface

struct ConvexHullProblem{N, M, T, O<:MOI.AbstractOptimizer, L}
    model::Model{T, O}
    point::MVector{N, T}
    vertices::MMatrix{N, M, T, L}
    weights::SVector{M, Variable}
    closest::SVector{N, Variable}

    function ConvexHullProblem{N, M, T}(optimizer::O) where {N, M, T, O<:MOI.AbstractOptimizer}
        model = Model(optimizer)
        L = N * M

        point = Parameter(model, val=zero(MVector{N, T}))
        vertices = Parameter(model, val=zero(MMatrix{N, M, T}))
        weights = SVector(ntuple(i -> Variable(model), M))
        closest = SVector(ntuple(i -> Variable(model), N))

        @constraint model vertices * weights == closest
        @constraint model weights >= zeros(M)
        @constraint model weights <= ones(M)
        @constraint model sum(weights) == 1
        difference = @expression point - closest
        @objective model Minimize difference ⋅ difference

        new{N, M, T, O, L}(model, point.val[], vertices.val[], weights, closest)
    end
end

function set_point!(problem::ConvexHullProblem{N, M, T}, point::AbstractVector) where {N, M, T}
    problem.point .= point
end

Base.@propagate_inbounds function set_vertex!(problem::ConvexHullProblem{N, M, T}, i::Integer, vertex::AbstractVector) where {N, M, T}
    problem.vertices[:, i] = vertex
end

function set_vertices!(problem::ConvexHullProblem{N, M, T}, vertices::AbstractVector) where {N, M, T}
    @boundscheck length(vertices) === M || error()
    @inbounds for i = Base.OneTo(M)
        set_vertex!(problem, i, vertices[i])
    end
end

solve!(problem::ConvexHullProblem) = Parametron.solve!(problem.model)

is_point_inside(problem::ConvexHullProblem; atol=0) = objectivevalue(problem.model) <= atol^2
distance_to_closest_point(problem::ConvexHullProblem) = sqrt(abs(objectivevalue(problem.model)))
closest_point(problem::ConvexHullProblem) = value.(problem.model, problem.closest)

end
