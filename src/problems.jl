import MathProgBase

export Problem, Solution, minimize, maximize, satisfy, add_constraint!, add_constraints!
export Float64OrNothing
export conic_problem

Float64OrNothing = Union(Float64, Nothing)

# TODO: Cleanup
type Solution{T<:Number}
  primal::Array{T, 1}
  dual::Array{T, 1}
  status::Symbol
  optval::T
  has_dual::Bool
end

Solution{T}(x::Array{T, 1}, status::Symbol, optval::T) = Solution(x, T[], status, optval, false)
Solution{T}(x::Array{T, 1}, y::Array{T, 1}, status::Symbol, optval::T) = Solution(x, y, status, optval, true)

type Problem
  head::Symbol
  objective::AbstractExpr
  constraints::Array{Constraint}
  status::Symbol
  optval::Float64OrNothing
  model::MathProgBase.AbstractMathProgModel
  solution::Solution

  function Problem(head::Symbol, objective::AbstractExpr,  
                   model::MathProgBase.AbstractMathProgModel, constraints::Array=Constraint[])
    return new(head, objective, constraints, "not yet solved", nothing, model)
  end
end
# constructor if model is not specified
function Problem(head::Symbol, objective::AbstractExpr, constraints::Array=Constraint[], 
                 solver::MathProgBase.AbstractMathProgSolver = get_default_solver())
  if solver == nothing
    error("The default solver is set to `nothing`
         You must have at least one solver installed.
         You can install a solver such as SCS by running:
         Pkg.add(\"SCS\").
         You will have to restart Julia after that.")
  end
  Problem(head, objective, MathProgBase.model(solver), constraints)
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
  var_to_ranges = Dict{Uint64, @compat Tuple{Int, Int}}()
  for constraint in constraints
    for i = 1:length(constraint.objs)
      for (id, val) in constraint.objs[i]
        if !haskey(var_to_ranges, id) && id != object_id(:constant)
          var = id_to_variables[id]
          var_to_ranges[id] = (index + 1, index + get_vectorized_size(var))
          index += get_vectorized_size(var)
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
  typeof(obj_vex) in bad_vex && warn("Problem not DCP compliant: objective is not DCP")

  constr_vex = ConstVexity()
  for i in 1:length(p.constraints)
    vex = vexity(p.constraints[i])
    typeof(vex) in bad_vex && warn("Problem not DCP compliant: constraint $i is not DCP")
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

function conic_problem(p::Problem)
  if get_vectorized_size(p.objective) != 1
    error("Objective must be a scalar")
  end
  # A map to hold unique constraints. Each constraint is keyed by a symbol
  # of which atom generated the constraints, and a integer hash of the child
  # expressions used by the atom
  unique_conic_forms = UniqueConicForms()
  objective, objective_var_id = conic_form!(p, unique_conic_forms)
  constraints = unique_conic_forms.constr_list
  var_size, constr_size, var_to_ranges = find_variable_ranges(constraints)
  c = spzeros(var_size, 1)
  objective_range = var_to_ranges[objective_var_id]
  c[objective_range[1]:objective_range[2]] = 1

  A = spzeros(constr_size, var_size)
  b = spzeros(constr_size, 1)
  cones = (@compat Tuple{Symbol, UnitRange{Int}})[]
  constr_index = 0
  for constraint in constraints
    total_constraint_size = 0
    for i = 1:length(constraint.objs)
      sz = constraint.sizes[i]
      for (id, val) in constraint.objs[i]
        if id == object_id(:constant)
          b[constr_index + 1 : constr_index + sz] = val
        else
          var_range = var_to_ranges[id]
          A[constr_index + 1 : constr_index + sz, var_range[1] : var_range[2]] = -val
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
  return c, A, b, cones, var_to_ranges, vartypes, constraints
end

Problem(head::Symbol, objective::AbstractExpr, constraints::Constraint...) =
  Problem(head, objective, [constraints...])

# Allow users to simply type minimize
minimize(objective::AbstractExpr, constraints::Constraint...) =
  Problem(:minimize, objective, [constraints...])
minimize{T<:Constraint}(objective::AbstractExpr, constraints::Array{T}=Constraint[]) =
  Problem(:minimize, objective, constraints)
minimize(objective::Value, constraints::Constraint...) =
  minimize(convert(AbstractExpr, objective), [constraints...])
minimize{T<:Constraint}(objective::Value, constraints::Array{T}=Constraint[]) =
  minimize(convert(AbstractExpr, objective), constraints)

# Allow users to simply type maximize
maximize(objective::AbstractExpr, constraints::Constraint...) =
  Problem(:maximize, objective, [constraints...])
maximize{T<:Constraint}(objective::AbstractExpr, constraints::Array{T}=Constraint[]) =
  Problem(:maximize, objective, constraints)
maximize(objective::Value, constraints::Constraint...) =
  maximize(convert(AbstractExpr, objective), [constraints...])
maximize{T<:Constraint}(objective::Value, constraints::Array{T}=Constraint[]) =
  maximize(convert(AbstractExpr, objective), constraints)

# Allow users to simply type satisfy (if there is no objective)
satisfy(constraints::Constraint...) = Problem(:minimize, Constant(0), [constraints...])
satisfy{T<:Constraint}(constraints::Array{T}=Constraint[]) =
  Problem(:minimize, Constant(0), constraints)
satisfy(constraint::Constraint) = satisfy([constraint])

# +(constraints, constraints) is defined in constraints.jl
add_constraints!{T<:Constraint}(p::Problem, constraints::Array{T}) = +(p.constraints, constraints)
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
