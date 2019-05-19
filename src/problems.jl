import MathProgBase

export Problem, Solution, minimize, maximize, satisfy, add_constraint!, add_constraints!
export Float64OrNothing

const Float64OrNothing = Union{Float64, Nothing}

# TODO: Cleanup
mutable struct Solution{T<:Number}
    primal::Array{T, 1}
    dual::Array{T, 1}
    status::Symbol
    optval::T
    has_dual::Bool
end

Solution(x::Array{T, 1}, status::Symbol, optval::T) where {T} =
    Solution(x, T[], status, optval, false)
Solution(x::Array{T, 1}, y::Array{T, 1}, status::Symbol, optval::T) where {T} =
    Solution(x, y, status, optval, true)

mutable struct Problem
    head::Symbol
    objective::AbstractExpr
    constraints::Array{Constraint}
    status::Symbol
    optval::Float64OrNothing
    model::Union{MathProgBase.AbstractConicModel, Nothing}
    solution::Solution

    function Problem(head::Symbol, objective::AbstractExpr,
                     model::Union{MathProgBase.AbstractConicModel, Nothing},
                     constraints::Array=Constraint[])
        if sign(objective)== Convex.ComplexSign()
            error("objective can not be a complex expression")
        else
            return new(head, objective, constraints, Symbol("not yet solved"), nothing, model)
        end
    end
end

# constructor if model is not specified
function Problem(head::Symbol, objective::AbstractExpr, constraints::Array=Constraint[],
                 solver::Union{MathProgBase.AbstractMathProgSolver, Nothing}=nothing)
    model = solver !== nothing ? MathProgBase.ConicModel(solver) : solver
    Problem(head, objective, model, constraints)
end

# If the problem constructed is of the form Ax=b where A is m x n
# returns:
# index: n
# constr_size: m
# var_to_ranges a dictionary mapping from variable id to (start_index, end_index)
# where start_index and end_index are the start and end indexes of the variable in A
function find_variable_ranges(constraints)
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

function conic_form!(p::Problem, unique_conic_forms::UniqueConicForms=UniqueConicForms())
    objective_var = Variable()
    objective = conic_form!(objective_var, unique_conic_forms)
    conic_form!(p.objective - objective_var == 0, unique_conic_forms)
    for constraint in p.constraints
        conic_form!(constraint, unique_conic_forms)
    end
    return objective, objective_var.id_hash
end

function conic_problem(p::Problem)
    if length(p.objective) != 1
        error("objective must be a scalar")
    end

    # conic problems have the form
    # minimize c'*x
    # st       b - Ax \in cones
    # our job is to take the conic forms of the objective and constraints
    # and convert them into vectors b and c and a matrix A
    # one chunk of rows in b and in A corresponds to each constraint,
    # and one chunk of columns in b and A corresponds to each variable,
    # with the size of the chunk determined by the size of the constraint or of the variable

    # A map to hold unique constraints. Each constraint is keyed by a symbol
    # of which atom generated the constraints, and a integer hash of the child
    # expressions used by the atom
    unique_conic_forms = UniqueConicForms()
    objective, objective_var_id = conic_form!(p, unique_conic_forms)
    constraints = unique_conic_forms.constr_list
    # var_to_ranges maps from variable id to the (start_index, stop_index) pairs of the columns of A corresponding to that variable
    # var_size is the sum of the lengths of all variables in the problem
    # constr_size is the sum of the lengths of all constraints in the problem
    var_size, constr_size, var_to_ranges = find_variable_ranges(constraints)
    c = spzeros(var_size, 1)
    objective_range = var_to_ranges[objective_var_id]
    c[objective_range[1]:objective_range[2]] .= 1

    # slot in all of the coefficients in the conic forms into A and b
    A = spzeros(constr_size, var_size)
    b = spzeros(constr_size, 1)
    cones = Tuple{Symbol, UnitRange{Int}}[]
    constr_index = 0
    for constraint in constraints
        total_constraint_size = 0
        for i = 1:length(constraint.objs)
            sz = constraint.sizes[i]
            for (id, val) in constraint.objs[i]
                if id == objectid(:constant)
                    for l in 1:sz
                        b[constr_index + l] = val[1][l] == 0 ? val[2][l] : val[1][l]
                    end
                    #b[constr_index + sz + 1 : constr_index + 2*sz] = val[2]
                else
                    var_range = var_to_ranges[id]
                    if id_to_variables[id].sign == ComplexSign()
                        A[constr_index + 1 : constr_index + sz, var_range[1] : var_range[1] + length(id_to_variables[id])-1] = -val[1]
                        A[constr_index + 1 : constr_index + sz, var_range[1] + length(id_to_variables[id]) : var_range[2]] = -val[2]
                    else
                        A[constr_index + 1 : constr_index + sz, var_range[1] : var_range[2]] = -val[1]
                    end
                end
            end
            constr_index += sz
            total_constraint_size += sz
        end
        push!(cones, (constraint.cone, constr_index - total_constraint_size + 1 : constr_index))
    end

    # find integral and boolean variables
    vartypes = fill(:Cont, length(c))
    for var_id in keys(var_to_ranges)
        variable = id_to_variables[var_id]
        if :Int in variable.sets
            startidx, endidx = var_to_ranges[var_id]
            for idx in startidx:endidx
                vartypes[idx] = :Int
            end
        end
        if :Bin in variable.sets
            startidx, endidx = var_to_ranges[var_id]
            for idx in startidx:endidx
                vartypes[idx] = :Bin
            end
        end
    end

    if p.head == :maximize
        c = -c
    end

    return c, A, b, cones, var_to_ranges, vartypes, constraints
end

Problem(head::Symbol, objective::AbstractExpr, constraints::Constraint...) =
    Problem(head, objective, [constraints...])

# Allow users to simply type minimize
minimize(objective::AbstractExpr, constraints::Constraint...) =
    Problem(:minimize, objective, collect(constraints))
minimize(objective::AbstractExpr, constraints::Array{<:Constraint}=Constraint[]) =
    Problem(:minimize, objective, constraints)
minimize(objective::Value, constraints::Constraint...) =
    minimize(convert(AbstractExpr, objective), collect(constraints))
minimize(objective::Value, constraints::Array{<:Constraint}=Constraint[]) =
    minimize(convert(AbstractExpr, objective), constraints)

# Allow users to simply type maximize
maximize(objective::AbstractExpr, constraints::Constraint...) =
    Problem(:maximize, objective, collect(constraints))
maximize(objective::AbstractExpr, constraints::Array{<:Constraint}=Constraint[]) =
    Problem(:maximize, objective, constraints)
maximize(objective::Value, constraints::Constraint...) =
    maximize(convert(AbstractExpr, objective), collect(constraints))
maximize(objective::Value, constraints::Array{<:Constraint}=Constraint[]) =
    maximize(convert(AbstractExpr, objective), constraints)

# Allow users to simply type satisfy (if there is no objective)
satisfy(constraints::Constraint...) = Problem(:minimize, Constant(0), [constraints...])
satisfy(constraints::Array{<:Constraint}=Constraint[]) =
    Problem(:minimize, Constant(0), constraints)
satisfy(constraint::Constraint) = satisfy([constraint])

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
