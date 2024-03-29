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
        end
        return new(
            head,
            objective,
            constraints,
            MOI.OPTIMIZE_NOT_CALLED,
            nothing,
        )
    end
end

function Base.getproperty(p::Problem, s::Symbol)
    if s === :optval
        if getfield(p, :status) == MOI.OPTIMIZE_NOT_CALLED
            return nothing
        else
            return objective_value(p)
        end
    end
    return getfield(p, s)
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
        if typeof(vex) in bad_vex
            @warn "Problem not DCP compliant: constraint $i is not DCP"
        end
        constr_vex += vex
    end
    problem_vex = obj_vex + constr_vex
    # this check is redundant
    if typeof(problem_vex) in bad_vex
        @warn("Problem not DCP compliant")
    end
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
    if typeof(obj_vex) in bad_vex
        @warn "Problem not DCP compliant: objective is not DCP"
    end
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

function monotonicity(p::Problem)
    if p.head === :satisfy
        error("Satisfiability problem cannot be used as subproblem")
    end
    m = monotonicity(p.objective)
    if p.head === :maximize
        return (-).(m)
    end
    return m
end

function curvature(p::Problem)
    if p.head === :satisfy
        error("Satisfiability problem cannot be used as subproblem")
    end
    m = curvature(p.objective)
    if p.head === :maximize
        return (-).(m)
    end
    return m
end

function Base.sign(p::Problem)
    if p.head === :satisfy
        error("Satisfiability problem cannot be used as subproblem")
    end
    m = sign(p.objective)
    if p.head === :maximize
        return (-).(m)
    end
    return m
end

function new_conic_form!(context::Context, p::Problem)
    for c in p.constraints
        add_constraint!(context, c)
    end
    if p.head == :satisfy
        return
    end
    return conic_form!(context, p.objective)
end

function _to_scalar_moi(::Type{T}, x) where {T}
    return _to_scalar_moi(T, only(MOI.Utilities.scalarize(x)))
end

function _to_scalar_moi(::Type{T}, x::SparseTape) where {T}
    return _to_scalar_moi(T, to_vaf(x))
end

function _to_scalar_moi(::Type{T}, x::Number) where {T<:Real}
    return MOI.ScalarAffineFunction{T}(MOI.ScalarAffineTerm{T}[], x)
end

_to_scalar_moi(::Type{T}, f::MOI.AbstractScalarFunction) where {T} = f

function Context(p::Problem{T}, optimizer_factory) where {T}
    context = Context{T}(optimizer_factory)
    cfp = conic_form!(context, p)
    # Save some memory before the solve;
    # we don't need these anymore,
    # since we don't re-use contexts between solves currently
    empty!(context.conic_form_cache)
    if p.head == :satisfy
        MOI.set(context.model, MOI.ObjectiveSense(), MOI.FEASIBILITY_SENSE)
    else
        obj = _to_scalar_moi(T, cfp)
        MOI.set(context.model, MOI.ObjectiveFunction{typeof(obj)}(), obj)
        MOI.set(
            context.model,
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

function minimize(
    objective::Value,
    constraints::Constraint...;
    numeric_type = Float64,
)
    return minimize(objective, collect(constraints); numeric_type)
end

function minimize(
    objective::Value,
    constraints::Array{<:Constraint} = Constraint[];
    numeric_type = Float64,
)
    return minimize(constant(objective), constraints; numeric_type)
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
    return maximize(objective, collect(constraints); numeric_type)
end

function maximize(
    objective::Value,
    constraints::Array{<:Constraint} = Constraint[];
    numeric_type = Float64,
)
    return maximize(constant(objective), constraints; numeric_type)
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

Base.:+(x::Array{<:Constraint}, y::Array{<:Constraint}) = vcat(x, y)

Base.:+(x::Constraint, y::Constraint) = [x, y]

Base.:+(x::Constraint, y::Array{<:Constraint}) = vcat(x, y)

Base.:+(x::Array{<:Constraint}, y::Constraint) = vcat(x, y)

iscomplex(c::Constraint) = iscomplex(c.lhs) || iscomplex(c.rhs)

function add_constraint!(context::Context, c::Constraint)
    if c in keys(context.constr_to_moi_inds)
        return
    end
    return _add_constraint!(context, c)
end

"""
    write_to_file(problem::Problem{Float64}, filename::String)

Write the current problem to the file at `filename`.

Requires solving the problem at least once using [`solve!`](@ref) to ensure that
the problem is loaded into a MathOptInterface model.

The file format is inferred from the filename extension. Supported file
types depend on the model type.

Currently, `Float64` is the only supported coefficient type. This may be
relaxed in future if file formats support other types.
"""
function write_to_file(p::Problem{T}, filename::String) where {T<:Float64}
    if p.model === nothing
        msg = """
        Problem has not been loaded into a MathOptInterface model;
        call `solve!(problem, optimizer)` before writing problem to file.
        """
        throw(ArgumentError(msg))
    end
    dest = MOI.FileFormats.Model(; filename)
    model = MOI.Bridges.full_bridge_optimizer(dest, T)
    MOI.copy_to(model, p.model)
    MOI.write_to_file(dest, filename)
    return
end
