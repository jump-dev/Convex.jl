mutable struct Problem{T<:Real} <: AbstractExpr
    head::Symbol
    objective::Union{AbstractExpr,Nothing}
    constraints::Array{Constraint}
    status::MOI.TerminationStatusCode
    model::Union{MOI.ModelLike,Nothing}
    id_hash::UInt64
    function Problem{T}(
        head::Symbol,
        objective::Union{AbstractExpr,Nothing},
        constraints::Array = Constraint[],
    ) where {T<:Real}
        if objective !== nothing && sign(objective) == Convex.ComplexSign()
            error("Objective cannot be a complex expression")
        else
            p = new(
                head,
                objective,
                constraints,
                MOI.OPTIMIZE_NOT_CALLED,
                nothing,
            )
            p.id_hash = objectid(p) # is this right?
            return p
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
objective_value(p::Problem) = MOI.get(p.model, MOI.ObjectiveValue())

Problem(args...) = Problem{Float64}(args...)

function vexity(p::Problem)
    bad_vex = [ConcaveVexity, NotDcp]

    obj_vex = vexity(p.objective)
    if p.head == :maximize
        obj_vex = -obj_vex
    end
    typeof(obj_vex) in bad_vex &&
        @warn "Problem not DCP compliant: objective is not DCP"

    constr_vex = ConstVexity()
    for i in 1:length(p.constraints)
        vex = vexity(p.constraints[i])
        typeof(vex) in bad_vex &&
            @warn "Problem not DCP compliant: constraint $i is not DCP"
        constr_vex += vex
    end
    # problem_vex = problem_vexity(p.head, obj_vex, constr_vex)
    problem_vex = obj_vex + constr_vex
    # this check is redundant
    # typeof(problem_vex) in bad_vex && warn("Problem not DCP compliant")
    return problem_vex
end

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

sign(p::Problem) = sign(p.objective)
monotonicity(p::Problem) = monotonicity(p.objective)

function template(p::Problem, context::Context)
    for c in p.constraints
        add_constraints_to_context(c, context)
    end
    if p.head !== :satisfy
        return template(p.objective, context)
    else
        return nothing
    end
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
function minimize(
    objective::Value,
    constraints::Constraint...;
    numeric_type = Float64,
)
    return satisfy(collect(constraints); numeric_type = numeric_type)
end
function minimize(
    objective::Value,
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
function maximize(
    objective::Value,
    constraints::Constraint...;
    numeric_type = Float64,
)
    return satisfy(collect(constraints); numeric_type = numeric_type)
end
function maximize(
    objective::Value,
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
