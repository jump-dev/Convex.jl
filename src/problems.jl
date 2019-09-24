using MathOptInterface
const MOI = MathOptInterface
const MOIU = MOI.Utilities
export Problem, minimize, maximize, satisfy, add_constraint!, add_constraints!

mutable struct Problem{T<:Real}
    head::Symbol
    objective::AbstractExpr
    constraints::Array{Constraint}
    status::MOI.TerminationStatusCode
    optval::Union{Real,Nothing}
    model::Union{MOI.ModelLike, Nothing}

    function Problem{T}(head::Symbol, objective::AbstractExpr,
                     constraints::Array=Constraint[]) where {T <: Real}
        if sign(objective)== Convex.ComplexSign()
            error("Objective cannot be a complex expression")
        else
            return new(head, objective, constraints, MOI.OPTIMIZE_NOT_CALLED, nothing, nothing)
        end
    end
end

dual_status(p::Problem) = get(p.model, MOI.DualStatus())
primal_status(p::Problem) = get(p.model, MOI.PrimalStatus())

Problem(args...) = Problem{Float64}(args...)

# If the problem constructed is of the form Ax=b where A is m x n
# returns:
# index: n
# constr_size: m
# var_to_ranges a dictionary mapping from variable id to (start_index, end_index)
# where start_index and end_index are the start and end indexes of the variable in A
function find_variable_ranges(constraints, id_to_variables)
    index = 0
    constr_size = 0
    var_to_ranges = Dict{UInt64, Tuple{Int, Int}}()
    for constraint in constraints
        for i = 1:length(constraint.objs)
            for (id, val) in constraint.objs[i]
                if !haskey(var_to_ranges, id) && id != objectid(:constant)
                    var = id_to_variables[id]
                    if var.sign == ComplexSign()
                        var_to_ranges[id] = (index + 1, index + 2*length(var))
                        index += 2*length(var)
                    else
                        var_to_ranges[id] = (index + 1, index + length(var))
                        index += length(var)
                    end
                end
            end
            constr_size += constraint.sizes[i]
        end
    end
    return index, constr_size, var_to_ranges
end

function vexity(p::Problem)
    bad_vex = [ConcaveVexity, NotDcp]

    obj_vex = vexity(p.objective)
    if p.head == :maximize
        obj_vex = -obj_vex
    end
    typeof(obj_vex) in bad_vex && @warn "Problem not DCP compliant: objective is not DCP"

    constr_vex = ConstVexity()
    for i in 1:length(p.constraints)
        vex = vexity(p.constraints[i])
        typeof(vex) in bad_vex && @warn "Problem not DCP compliant: constraint $i is not DCP"
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

Problem{T}(head::Symbol, objective::AbstractExpr, constraints::Constraint...) where {T<:Real} =
    Problem{T}(head, objective, [constraints...])

# Allow users to simply type minimize
minimize(objective::AbstractExpr, constraints::Constraint...; numeric_type = Float64) =
    Problem{numeric_type}(:minimize, objective, collect(constraints))
minimize(objective::AbstractExpr, constraints::Array{<:Constraint}=Constraint[]; numeric_type = Float64) =
    Problem{numeric_type}(:minimize, objective, constraints)
minimize(objective::Value, constraints::Constraint...; numeric_type = Float64) =
    minimize(convert(AbstractExpr, objective), collect(constraints); numeric_type = numeric_type)
minimize(objective::Value, constraints::Array{<:Constraint}=Constraint[]; numeric_type = Float64) =
    minimize(convert(AbstractExpr, objective), constraints; numeric_type = numeric_type)

# Allow users to simply type maximize
maximize(objective::AbstractExpr, constraints::Constraint...; numeric_type = Float64) =
    Problem{numeric_type}(:maximize, objective, collect(constraints))
maximize(objective::AbstractExpr, constraints::Array{<:Constraint}=Constraint[]; numeric_type = Float64) =
    Problem{numeric_type}(:maximize, objective, constraints)
maximize(objective::Value, constraints::Constraint...; numeric_type = Float64) =
    maximize(convert(AbstractExpr, objective), collect(constraints); numeric_type = numeric_type)
maximize(objective::Value, constraints::Array{<:Constraint}=Constraint[]; numeric_type = Float64) =
    maximize(convert(AbstractExpr, objective), constraints; numeric_type = numeric_type)

# Allow users to simply type satisfy (if there is no objective)
satisfy(constraints::Constraint...; numeric_type = Float64) = Problem{numeric_type}(:minimize, Constant(0), [constraints...])
satisfy(constraints::Array{<:Constraint}=Constraint[]; numeric_type = Float64) =
    Problem{numeric_type}(:minimize, Constant(0), constraints)
satisfy(constraint::Constraint; numeric_type = Float64) = satisfy([constraint]; numeric_type = numeric_type)

# +(constraints, constraints) is defined in constraints.jl
add_constraints!(p::Problem, constraints::Array{<:Constraint}) =
    +(p.constraints, constraints)
add_constraints!(p::Problem, constraint::Constraint) = add_constraints!(p, [constraint])
add_constraint! = add_constraints!

# caches conic form of x when x is the solution to the optimization problem p
function cache_conic_form!(conic_forms::UniqueConicForms, x::AbstractExpr, p::Problem)
    objective = conic_form!(p.objective, conic_forms)
    for c in p.constraints
        conic_form!(c, conic_forms)
    end
    cache_conic_form!(conic_forms, x, objective)
end
