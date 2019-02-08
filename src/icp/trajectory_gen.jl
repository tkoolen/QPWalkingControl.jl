struct ICPTrajectoryGenerator{T, N, O<:MOI.AbstractOptimizer, L}
    model::Model{T, O}
    cops::Vector{SVector{2, Variable}}
    icps::Vector{SVector{2, Variable}}
    Δts::Vector{T}
    ωs::Vector{T}
    cop_polyhedra::Vector{MHRep{N, 2, T, L}}
    preferred_cops::Vector{Vec2{T}}
    initial_icp::Base.RefValue{Vec2{T}}
    final_icp::Base.RefValue{Vec2{T}}
    num_active_segments::Base.RefValue{Int}

    function ICPTrajectoryGenerator{T, N}(optimizer::O, num_segments::Integer) where {T, N, O<:MOI.AbstractOptimizer}
        n = num_segments
        L = N * 2

        model = Model(optimizer)

        # Variables
        cops = [SVector(ntuple(i -> Variable(model), 2)) for i = 1 : n]
        icps = [SVector(ntuple(i -> Variable(model), 2)) for i = 1 : n]

        # Parameters
        Δts = Parameter(model, val=zeros(T, n))
        ωs = Parameter(model, val=zeros(T, n))
        cop_polyhedra = Parameter(model, val=[zero(MHRep{N, 2, T, L}) for i = 1 : n])
        preferred_cops = Parameter(model, val=[zero(Vec2{T}) for i = 1 : n])
        initial_icp = Parameter(model, val=zero(Vec2{T}))
        final_icp = Parameter(model, val=zero(Vec2{T}))

        # CoP bounds
        for i = 1 : n
            A = @expression cop_polyhedra[i].A
            b = @expression cop_polyhedra[i].b
            @constraint model A * cops[i] <= b
        end

        # Initial ICP
        @constraint model icps[1] == initial_icp

        # ICP dynamics
        for i = 1 : n
            w = @expression exp(ωs[i] * Δts[i])
            next_icp = @expression w * icps[i] + (1 - w) * cops[i]
            if i < n
                @constraint model next_icp == icps[i + 1]
            else
                @constraint model next_icp == final_icp
            end
        end

        # Objective
        objective = Parametron.LazyExpression(identity, zero(QuadraticFunction{Float64}))
        for i = 1 : n
            e = @expression cops[i] - preferred_cops[i]
            objective = @expression objective + e ⋅ e
        end
        @objective model Minimize objective

        num_active_segments = Ref(0)

        new{T, N, O, L}(
            model,
            cops, icps,
            Δts.val[], ωs.val[], cop_polyhedra.val[], preferred_cops.val[],
            initial_icp.val, final_icp.val,
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

function push_segment!(generator::ICPTrajectoryGenerator,
        Δt::Number, ω::Number, cop_polyhedron::Union{ConvexHull, HRep}, preferred_cop::SVector)
    i = generator.num_active_segments[] += 1
    @boundscheck i <= length(generator.Δts) || error()
    @inbounds begin
        generator.Δts[i] = Δt
        generator.ωs[i] = ω
        if cop_polyhedron isa HRep
            generator.cop_polyhedra[i] = cop_polyhedron
        else
            hrep = generator.cop_polyhedra[i]
            hrep!(hrep.A, hrep.b, cop_polyhedron)
        end
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
    value.(generator.model, generator.icps[segment_number])
end

function cop(generator::ICPTrajectoryGenerator, segment_number::Integer)
    value.(generator.model, generator.cops[segment_number])
end

final_time(generator::ICPTrajectoryGenerator) = sum(generator.Δts)

function find_segment(generator::ICPTrajectoryGenerator{T}, t::Number) where T
    i = 1
    t0 = zero(T)
    for Δt in generator.Δts
        if t <= t0 + Δt
            break
        else
            i += 1
            t0 += Δt
        end
    end
    i, t0
end

function (generator::ICPTrajectoryGenerator{T})(t::Number) where T
    t < 0 && throw(ArgumentError("Cannot evaluate at a negative time."))

    if t < sum(generator.Δts)
        # Find segment index and segment start time
        i, t0 = find_segment(generator, t)

        # Use (13) in Capturability-based Analysis Part 1 to find ICP ξ
        ξ0 = initial_icp(generator, i)
        p = cop(generator, i)
        Δt = t - t0
        ω = generator.ωs[i]
        w = exp(ω * Δt)
        ξ = w * ξ0 + (1 - w) * p

        # Use (12) in Capturability-based Analysis Part 1 to find ICP time derivative ξd
        ξd = ω * (ξ - p)
    else
        # Hold at the final ICP
        ξ = generator.final_icp[]
        ξd = zero(ξ)
    end

    ξ, ξd
end
