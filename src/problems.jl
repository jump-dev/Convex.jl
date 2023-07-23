mutable struct Problem{T<:Real}
    head::Symbol
    objective::AbstractExpr
    constraints::Array{Constraint}
    status::MOI.TerminationStatusCode
    model::Union{MOI.ModelLike,Nothing}

    function Problem{T}(
        head::Symbol,
        objective::AbstractExpr,
        constraints::Array = Constraint[],
    ) where {T<:Real}
        if sign(objective) == Convex.ComplexSign()
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
    problem_vex = obj_vex + constr_vex
    # this check is redundant
    # typeof(problem_vex) in bad_vex && warn("Problem not DCP compliant")
    return problem_vex
end

function conic_form!(p::Problem, unique_conic_forms::UniqueConicForms)
    objective_var = Variable()
    objective = conic_form!(objective_var, unique_conic_forms)
    conic_form!(p.objective - objective_var == 0, unique_conic_forms)
    for constraint in p.constraints
        conic_form!(constraint, unique_conic_forms)
    end
    return objective, objective_var.id_hash
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
    return minimize(
        convert(AbstractExpr, objective),
        collect(constraints);
        numeric_type = numeric_type,
    )
end
function minimize(
    objective::Value,
    constraints::Array{<:Constraint} = Constraint[];
    numeric_type = Float64,
)
    return minimize(
        convert(AbstractExpr, objective),
        constraints;
        numeric_type = numeric_type,
    )
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
    return maximize(
        convert(AbstractExpr, objective),
        collect(constraints);
        numeric_type = numeric_type,
    )
end
function maximize(
    objective::Value,
    constraints::Array{<:Constraint} = Constraint[];
    numeric_type = Float64,
)
    return maximize(
        convert(AbstractExpr, objective),
        constraints;
        numeric_type = numeric_type,
    )
end

# Allow users to simply type satisfy (if there is no objective)
function satisfy(constraints::Constraint...; numeric_type = Float64)
    return Problem{numeric_type}(:minimize, Constant(0), [constraints...])
end
function satisfy(
    constraints::Array{<:Constraint} = Constraint[];
    numeric_type = Float64,
)
    return Problem{numeric_type}(:minimize, Constant(0), constraints)
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

# caches conic form of x when x is the solution to the optimization problem p
function cache_conic_form!(
    conic_forms::UniqueConicForms,
    x::AbstractExpr,
    p::Problem,
)
    objective = conic_form!(p.objective, conic_forms)
    for c in p.constraints
        conic_form!(c, conic_forms)
    end
    return cache_conic_form!(conic_forms, x, objective)
end

"""
    write_to_file(problem::Problem, filename::String)

Write the current problem to the file at `filename`. Requires solving
the problem at least once using [`solve!`](@ref) to ensure that the
problem is loaded into a MathOptInterface model.

The file format is inferred from the filename extension. Supported file
types depend on the model type.
"""
function write_to_file(problem::Problem{T}, filename::String) where {T}
    isnothing(problem.model) && throw(
        ArgumentError(
            """
            Problem has not been loaded into a MathOptInterface model; 
            call `solve!(problem, optimizer)` before writing problem to file.
            """,
        ),
    )
    dest = MOI.FileFormats.Model(filename = filename)
    src = problem.model
    MOI.copy_to(MOI.Bridges.full_bridge_optimizer(dest, T), src)
    return MOI.write_to_file(dest, filename) # nothing
end
