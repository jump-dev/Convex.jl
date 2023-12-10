mutable struct Problem{T<:Real} <: AbstractExpr
    head::Symbol
    objective::Union{AbstractExpr,Nothing}
    constraints::Array{Constraint}
    status::MOI.TerminationStatusCode
    model::Union{MOI.ModelLike,Nothing}
    function Problem{T}(
        head::Symbol,
        objective::Union{AbstractExpr,Nothing},
        constraints::Array = Constraint[],
    ) where {T<:Real}
        if objective !== nothing && sign(objective) == Convex.ComplexSign()
            error("Objective cannot be a complex expression")
        else
            return new(
                head,
                objective,
                constraints,
                MOI.OPTIMIZE_NOT_CALLED,
                nothing,
            )
        end
    end
end

function Base.getproperty(p::Problem, s::Symbol)
    if s === :optval
        if getfield(p, :status) == MOI.OPTIMIZE_NOT_CALLED
            return nothing
        else
            return objective_value(p)
        end
    elseif s === :size
        return (1, 1)
    else
        return getfield(p, s)
    end
end

dual_status(p::Problem) = MOI.get(p.model, MOI.DualStatus())
primal_status(p::Problem) = MOI.get(p.model, MOI.PrimalStatus())
termination_status(p::Problem) = MOI.get(p.model, MOI.TerminationStatus())
function objective_value(p::Problem)
    # These don't have an objective value, and it would be confusing to return one
    if p.head === :satisfy
        return nothing
    end
    return MOI.get(p.model, MOI.ObjectiveValue())
end

Problem(args...) = Problem{Float64}(args...)

function problem_vexity(p::Problem)
    bad_vex = (ConcaveVexity, NotDcp)
    obj_vex = objective_vexity(p)
    constr_vex = ConstVexity()
    for i in 1:length(p.constraints)
        vex = vexity(p.constraints[i])
        typeof(vex) in bad_vex &&
            @warn "Problem not DCP compliant: constraint $i is not DCP"
        constr_vex += vex
    end
    problem_vex = obj_vex + constr_vex
    # this check is redundant
    typeof(problem_vex) in bad_vex && @warn("Problem not DCP compliant")
    return problem_vex
end

function objective_vexity(p::Problem)
    bad_vex = (ConcaveVexity, NotDcp)
    if p.head == :satisfy
        obj_vex = ConstVexity()
    elseif p.head == :minimize
        obj_vex = vexity(p.objective)
    elseif p.head == :maximize
        obj_vex = -vexity(p.objective)
    else
        error("Unknown type of problem $(p.head)")
    end

    typeof(obj_vex) in bad_vex &&
        @warn "Problem not DCP compliant: objective is not DCP"

    return obj_vex
end

vexity(p::Problem) = objective_vexity(p)

# is this right?
# function problem_vexity(sense, obj_vexity, constr_vexity)
#     if constr_vexity == ConstVexity()
#         return obj_vexity
#     end

#     if constr_vexity == AffineVexity() && obj_vexity == AffineVexity()
#         return sense === :maximize ? ConcaveVexity() : ConvexVexity()
#     end

#     return obj_vexity + constr_vexity
# end

for f in (:sign, :monotonicity, :curvature)
    @eval function $f(p::Problem)
        if p.head === :satisfy
            error("Satisfiability problem cannot be used as subproblem")
        end
        m = $f(p.objective)
        if p.head === :maximize
            return (-).(m)
        else
            return m
        end
    end
end

function new_conic_form!(context::Context, p::Problem)
    for c in p.constraints
        add_constraint!(context, c)
    end
    if p.head !== :satisfy
        return conic_form!(context, p.objective)
    else
        return nothing
    end
end

function Context(problem::Problem{T}, optimizer_factory; kwargs...) where {T}
    optimizer = MOI.instantiate(optimizer_factory)
    return Context(problem, optimizer; kwargs...)
end

function Context(p::Problem{T}, optimizer::MOI.ModelLike) where {T}
    context = Context{T}(optimizer)
    cfp = conic_form!(context, p)

    # Save some memory before the solve;
    # we don't need these anymore,
    # since we don't re-use contexts between solves currently
    empty!(context.conic_form_cache)

    model = context.model

    if p.head == :satisfy
        MOI.set(model, MOI.ObjectiveSense(), MOI.FEASIBILITY_SENSE)
    else
        obj = scalar_fn(cfp)
        MOI.set(model, MOI.ObjectiveFunction{typeof(obj)}(), obj)
        MOI.set(
            model,
            MOI.ObjectiveSense(),
            p.head == :maximize ? MOI.MAX_SENSE : MOI.MIN_SENSE,
        )
    end
    return context
end

function Problem{T}(
    head::Symbol,
    objective::AbstractExpr,
    constraints::Constraint...,
) where {T<:Real}
    return Problem{T}(head, objective, [constraints...])
end

# Allow users to simply type minimize
function minimize(
    objective::AbstractExpr,
    constraints::Constraint...;
    numeric_type = Float64,
)
    return Problem{numeric_type}(:minimize, objective, collect(constraints))
end
function minimize(
    objective::AbstractExpr,
    constraints::Array{<:Constraint} = Constraint[];
    numeric_type = Float64,
)
    return Problem{numeric_type}(:minimize, objective, constraints)
end
function minimize(::Value, constraints::Constraint...; numeric_type = Float64)
    return satisfy(collect(constraints); numeric_type = numeric_type)
end
function minimize(
    ::Value,
    constraints::Array{<:Constraint} = Constraint[];
    numeric_type = Float64,
)
    return satisfy(constraints; numeric_type = numeric_type)
end

# Allow users to simply type maximize
function maximize(
    objective::AbstractExpr,
    constraints::Constraint...;
    numeric_type = Float64,
)
    return Problem{numeric_type}(:maximize, objective, collect(constraints))
end
function maximize(
    objective::AbstractExpr,
    constraints::Array{<:Constraint} = Constraint[];
    numeric_type = Float64,
)
    return Problem{numeric_type}(:maximize, objective, constraints)
end
function maximize(::Value, constraints::Constraint...; numeric_type = Float64)
    return satisfy(collect(constraints); numeric_type = numeric_type)
end
function maximize(
    ::Value,
    constraints::Array{<:Constraint} = Constraint[];
    numeric_type = Float64,
)
    return satisfy(constraints; numeric_type = numeric_type)
end

# Allow users to simply type satisfy (if there is no objective)
function satisfy(constraints::Constraint...; numeric_type = Float64)
    return Problem{numeric_type}(:satisfy, nothing, [constraints...])
end
function satisfy(
    constraints::Array{<:Constraint} = Constraint[];
    numeric_type = Float64,
)
    return Problem{numeric_type}(:satisfy, nothing, constraints)
end
function satisfy(constraint::Constraint; numeric_type = Float64)
    return satisfy([constraint]; numeric_type = numeric_type)
end

# +(constraints, constraints) is defined in constraints.jl
function add_constraints!(p::Problem, constraints::Array{<:Constraint})
    return append!(p.constraints, constraints)
end
function add_constraints!(p::Problem, constraint::Constraint)
    return add_constraints!(p, [constraint])
end

function add_constraint!(p::Problem, constraints::Array{<:Constraint})
    return add_constraints!(p, constraints)
end
function add_constraint!(p::Problem, constraint::Constraint)
    return add_constraints!(p, constraint)
end
