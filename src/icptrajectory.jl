struct ICPTrajectoryGenerator{T, M, O<:MOI.AbstractOptimizer, L}
    model::Model{T, O}
    cops::Vector{SVector{2, Variable}}
    icps::Vector{SVector{2, Variable}}
    Δts::Vector{T}
    ωs::Vector{T}
    cop_polyhedra::Vector{SHRep{M, 2, T, L}}
    preferred_cops::Vector{SVector{2, T}}
    initial_icp::Base.RefValue{SVector{2, T}}
    final_icp::Base.RefValue{SVector{2, T}}

    function ICPTrajectoryGenerator{T, M}(optimizer::O, num_segments::Integer) where {T, M, O<:MOI.AbstractOptimizer}
        n = num_segments
        L = M * 2

        model = Model(optimizer)

        # Variables
        cops = [SVector(ntuple(i -> Variable(model), 2)) for i = 1 : n]
        icps = [SVector(ntuple(i -> Variable(model), 2)) for i = 1 : n]

        # Parameters
        Δts = Parameter(model, val=zeros(T, n))
        ωs = Parameter(model, val=zeros(T, n))
        cop_polyhedra = Parameter(model, val=[zero(SHRep{M, 2, T, L}) for i = 1 : n])
        preferred_cops = Parameter(model, val=[zero(SVector{2, T}) for i = 1 : n])
        initial_icp = Parameter(model, val=zero(SVector{2, T}))
        final_icp = Parameter(model, val=zero(SVector{2, T}))

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
            Δt = @expression Δts[i]
            ω = @expression ωs[i]
            w = @expression exp(ω * Δt)
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

        new{T, M, O, L}(
            model,
            cops, icps,
            Δts.val[], ωs.val[], cop_polyhedra.val[], preferred_cops.val[],
            initial_icp.val, final_icp.val
        )
    end
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
            @show terminationstatus(model)
            @show primalstatus(model)
            @show dualstatus(model)
            println()
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
    # Find segment index and segment start time
    i, t0 = find_segment(generator, t)
    num_segments = length(generator.Δts)

    if i <= num_segments
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
