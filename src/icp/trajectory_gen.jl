function integrate_icp(ξ0, p::BezierCurve, ω, θ)
    exp_ω_θ = exp(ω * θ)
    return ξ0 * exp_ω_θ - ω * exp_ω_θ * exponential_integral(p, -ω, θ)
end

const COP_TRAJ_DEGREE = 4

struct ICPTrajectoryGenerator{T, N, O<:MOI.AbstractOptimizer, L}
    ω::T
    model::Model{T, O}
    cop_pieces::Vector{BezierCurve{COP_TRAJ_DEGREE + 1, SVector{2, Variable}}}
    icp_knots::Vector{SVector{2, Variable}}
    Δts::Vector{T}
    cop_polyhedra::Vector{MHRep{N, 2, T, L}}
    preferred_cops::Vector{Vec2{T}}
    initial_cop::Base.RefValue{Vec2{T}}
    initial_icp::Base.RefValue{Vec2{T}}
    final_icp::Base.RefValue{Vec2{T}}
    num_active_segments::Base.RefValue{Int}

    function ICPTrajectoryGenerator{T, N}(optimizer::O, num_segments::Integer, ω::Number) where {T, N, O<:MOI.AbstractOptimizer}
        # TODO:
        # * find the right objective function

        # Variables: (2 * COP_TRAJ_DEGREE + 4) * n + 2:
        # * CoP control points: 2 * (COP_TRAJ_DEGREE + 1) * n
        # * ICP knots: 2 * (n + 1)
        # Equality constraints: 6 * n + 4:
        # * initial CoP: 2
        # * CoP C⁰ continuity: 2 * (n - 1)
        # * CoP C¹ continuity: 2 * (n - 1)
        # * ICP dynamics: 2 * n
        # * initial ICP: 2
        # * final ICP: 2
        # * final CoP == final ICP: 2

        n = num_segments
        L = N * 2
        model = Model(optimizer)

        # Parameters
        Δts = Parameter(model, val=zeros(T, n))
        cop_polyhedra = Parameter(model, val=[zero(MHRep{N, 2, T, L}) for i = 1 : n])
        preferred_cops = Parameter(model, val=[zero(Vec2{T}) for i = 1 : n])
        initial_cop = Parameter(model, val=zero(Vec2{T}))
        initial_icp = Parameter(model, val=zero(Vec2{T}))
        final_icp = Parameter(model, val=zero(Vec2{T}))

        # CoP trajectory
        cop_pieces = [BezierCurve(ntuple(i -> SVector(Variable(model), Variable(model)), Val(COP_TRAJ_DEGREE + 1))) for i = 1 : n]

        # Initial CoP
        @constraint model first(first(cop_pieces).points) == initial_cop

        # CoP trajectory continuity
        for i = 1 : n - 1
            # C⁰
            pᵢ = cop_pieces[i]
            pᵢ₊₁ = cop_pieces[i + 1]
            @constraint model first(pᵢ₊₁.points) == last(pᵢ.points)

            # C¹. Note: clearing denominators.
            pᵢ′ = derivative(pᵢ)
            pᵢ₊₁′ = derivative(pᵢ₊₁)
            # @constraint model first(pᵢ₊₁′.points) * Δts[i] == last(pᵢ′.points) * Δts[i + 1] #  # TODO: allocates. Fix missing optimization for SVector * scalar.
            @constraint model Vector(first(pᵢ₊₁′.points)) * Δts[i] == Vector(last(pᵢ′.points)) * Δts[i + 1]
        end

        # CoP bounds
        for i = 1 : n
            Aᵢ = @expression cop_polyhedra[i].A
            bᵢ = @expression cop_polyhedra[i].b
            pᵢ = cop_pieces[i]
            for point in pᵢ.points
                @constraint model Aᵢ * point <= bᵢ
            end
        end

        # ICP dynamics
        icp_knots = [SVector(Variable(model), Variable(model)) for i = 1 : n + 1]
        for i = 1 : n
            # Use an indicator variable to make icp knots equal if the Δt is zero
            z = @expression Δts[i] != 0
            next_icp = icp_knots[i] * exp(ω) - ω * exp(ω) * exponential_integral(cop_pieces[i] , -ω)
            @constraint model icp_knots[i + 1] == z * Vector(next_icp) + (1 - z) * Vector(icp_knots[i]) # TODO: don't use Vector
        end
        @constraint model first(icp_knots) == initial_icp
        @constraint model last(icp_knots) == final_icp
        @constraint model last(last(cop_pieces).points) == final_icp # implies that final ICP velocity is zero

        # Objective
        objective = Parametron.LazyExpression(identity, zero(QuadraticFunction{Float64}))
        for i = 1 : n
            pᵢ = cop_pieces[i]
            z = @expression Δts[i] != 0 # inactive segments shouldn't incur cost
            for point in pᵢ.points
                e = @expression point - preferred_cops[i]
                objective = @expression objective + z * (e ⋅ e)
            end
        end
        @objective model Minimize objective

        num_active_segments = Ref(0)

        new{T, N, O, L}(
            ω,
            model,
            cop_pieces, icp_knots,
            Δts.val[], cop_polyhedra.val[], preferred_cops.val[],
            initial_cop.val, initial_icp.val, final_icp.val,
            num_active_segments
        )
    end
end

function Base.empty!(generator::ICPTrajectoryGenerator)
    for hrep in generator.cop_polyhedra
        hrep.A .= 0
        hrep.b .= 0
    end
    generator.Δts .= 0
    generator.num_active_segments[] = 0
    generator
end

# TODO: name?
function initialize!(generator::ICPTrajectoryGenerator, initial_cop::AbstractVector, initial_icp::AbstractVector, final_icp::AbstractVector)
    generator.initial_cop[] = initial_cop
    generator.initial_icp[] = initial_icp
    generator.final_icp[] = final_icp
    generator
end

function push_segment!(generator::ICPTrajectoryGenerator,
        Δt::Number, cop_polyhedron::ConvexHull, preferred_cop::SVector = centroid(cop_polyhedron))
    i = generator.num_active_segments[] += 1
    @boundscheck i <= length(generator.Δts) || error()
    @inbounds begin
        generator.Δts[i] = Δt
        hrep = generator.cop_polyhedra[i]
        hrep!(hrep.A, hrep.b, cop_polyhedron)
        generator.preferred_cops[i] = preferred_cop
    end
    generator
end

function solve!(generator::ICPTrajectoryGenerator)
    Parametron.solve!(generator.model)
    checkstatus(generator.model)
end

function checkstatus(model::Parametron.Model)
    ok = terminationstatus(model) == MOI.Success && primalstatus(model) == MOI.FeasiblePoint
    if !ok
        okish = terminationstatus(model) == MOI.AlmostSuccess && primalstatus(model) == MOI.UnknownResultStatus
        if !okish
            """
            Solve failed.
            * termination status: $(terminationstatus(model))
            * primal status: $(primalstatus(model))
            * dual status: $(dualstatus(model))
            """ |> error
        end
    end
end

function initial_icp(generator::ICPTrajectoryGenerator, segment_number::Integer)
    value.(generator.model, generator.icp_knots[segment_number])
end

function cop_piece(generator::ICPTrajectoryGenerator, segment_number::Integer)
    BezierCurve(map(x -> value.(generator.model, x), generator.cop_pieces[segment_number].points))
end

final_time(generator::ICPTrajectoryGenerator) = sum(generator.Δts)

function find_segment(generator::ICPTrajectoryGenerator{T}, t::Number) where T
    i = 1
    t0 = zero(T)
    for Δt in generator.Δts
        if t <= t0 + Δt || i >= generator.num_active_segments[]
            break
        else
            i += 1
            t0 += Δt
        end
    end
    i, t0
end

num_segments(generator::ICPTrajectoryGenerator) = length(generator.Δts)

function (generator::ICPTrajectoryGenerator{T})(t::Number) where T
    t < 0 && throw(ArgumentError("Cannot evaluate at a negative time."))
    if t < sum(generator.Δts)
        # Find segment index and segment start time
        i, t0 = find_segment(generator, t)

        # Segment phase
        θ = (t - t0) / generator.Δts[i]

        # CoP
        pᵢ = cop_piece(generator, i)
        p = pᵢ(θ)

        # ICP
        ω = generator.ω
        ξ0 = initial_icp(generator, i)
        ξ = integrate_icp(ξ0, pᵢ, ω, θ)

        # ICP derivative. Use (12) in Capturability-based Analysis Part 1.
        ξd = ω * (ξ - p)
    else
        # Hold at the final ICP
        ξ = generator.final_icp[]
        ξd = zero(ξ)
    end
    ξ, ξd
end

function cop(generator::ICPTrajectoryGenerator, t::Number)
    ξ, ξd = generator(t)
    ω = generator.ω
    return ξ - ξd / ω
end
