export Problem, Solution, minimize, maximize, satisfy, add_constraints
export Float64OrNothing
export dual_conic_problem

Float64OrNothing = Union(Float64, Nothing)

# TODO: Cleanup
type Solution{T<:Number}
  primal::Array{T, 1}
  dual::Array{T, 1}
  slack::Array{T, 1}
  status::Symbol
  optval::T
  has_dual::Bool
end

Solution{T}(x::Array{T, 1}, status::Symbol, optval::T) = Solution(x, T[], T[], status, optval, false)
Solution{T}(x::Array{T, 1}, y::Array{T, 1}, z::Array{T, 1}, status::Symbol, optval::T) = Solution(x, y, z, status, optval, true)

type Problem
  head::Symbol
  objective::AbstractExpr
  constraints::Array{Constraint}
  status::Symbol
  optval::Float64OrNothing
  solution::Solution

  function Problem(head::Symbol, objective::AbstractExpr, constraints::Array=Constraint[])
    return new(head, objective, constraints, "not yet solved", nothing)
  end
end

function find_variable_ranges(constraints)
  index = 0
  constr_size = 0
  var_to_ranges = Dict{Uint64, (Int64, Int64)}()
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

function dual_conic_form(p::Problem, unique_constr)
  objective_var = Variable()
  objective, _ = dual_conic_form(objective_var, unique_constr)
  _, constraints = dual_conic_form(p.objective - objective_var == 0, unique_constr)
  for constraint in p.constraints
    append!(constraints, dual_conic_form(constraint, unique_constr)[2])
  end
  return objective, constraints, objective_var.id
end

function dual_conic_problem(p::Problem)
  unique_constr = Dict{(Symbol, Uint64), (ConicObj, Array{ConicConstr})}()
  objective, constraints, objective_var_id = dual_conic_form(p, unique_constr)
  var_size, constr_size, var_to_ranges = find_variable_ranges(constraints)

  c = spzeros(var_size, 1)
  objective_range = var_to_ranges[objective_var_id]
  c[objective_range[1]:objective_range[2]] = 1

  A = spzeros(constr_size, var_size)
  b = spzeros(constr_size, 1)
  cones = (Symbol, UnitRange{Int64})[]
  constr_index = 0
  for constraint in constraints
    total_constraint_size = 0
    for i = 1:length(constraint.objs)
      sz = constraint.sizes[i]
      for (id, val) in constraint.objs[i]
        if id == object_id(:constant)
          try
            b[constr_index + 1 : constr_index + sz] = val
          catch
            if size(val) == (1,1)
              v = val[1]
              for i=1:sz
                b[constr_index + i] = v
              end
            else
              error("Sizes don't match: $sz vs $(size(val))")
            end
          end
        else
          var_range = var_to_ranges[id]
          try
            A[constr_index + 1 : constr_index + sz, var_range[1] : var_range[2]] = -val
          catch
            if size(val) == (1,1)
              v = -val[1]
              for varidx in var_range[1]:var_range[2]
                for i=1:sz
                  A[constr_index + i, varidx] = v
                end
              end
            else
              error("Sizes don't match: $sz * $var_range vs $(size(val))")
            end
          end
        end
      end
      constr_index += sz
      total_constraint_size += sz
    end
    push!(cones, (constraint.cone, constr_index - total_constraint_size + 1 : constr_index))
  end
  return c, A, b, cones
end

Problem(head::Symbol, objective::AbstractExpr, constraints::Constraint...) =
  Problem(head, objective, [constraints...])

# Allow users to simply type minimize or maximize
minimize(objective::AbstractExpr, constraints::Constraint...) =
  Problem(:minimize, objective, [constraints...])
minimize(objective::AbstractExpr, constraints::Array{Constraint}=Constraint[]) =
  Problem(:minimize, objective, constraints)
minimize(objective::Value, constraints::Constraint...) =
  minimize(convert(AbstractExpr, objective), constraints)
minimize(objective::Value, constraints::Array{Constraint}=Constraint[]) =
  minimize(convert(AbstractExpr, objective), constraints)

maximize(objective::AbstractExpr, constraints::Constraint...) =
  Problem(:maximize, objective, [constraints...])
maximize(objective::AbstractExpr, constraints::Array{Constraint}=Constraint[]) =
  Problem(:maximize, objective, constraints)
maximize(objective::Value, constraints::Constraint...) =
  maximize(convert(AbstractExpr, objective), constraints)
maximize(objective::Value, constraints::Array{Constraint}=Constraint[]) =
  maximize(convert(AbstractExpr, objective), constraints)

satisfy(constraints::Array{Constraint}=Constraint[]) =
  Problem(:minimize, Constant(0), constraints)
satisfy(constraint::Constraint) = satisfy([constraint])

# +(constraints, constraints) is overwritten in constraints.jl
add_constraints(p::Problem, constraints::Array{Constraint}) = +(p.constraints, constraints)
add_constraints(p::Problem, constraint::Constraint) = add_constraints(p, [constraint])
