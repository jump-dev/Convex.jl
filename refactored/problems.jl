export Problem
export Float64OrNothing
export dual_conic_problem

Float64OrNothing = Union(Float64, Nothing)


type Problem
  head::Symbol
  objective::AbstractExpr
  constraints::Array{Constraint}
  status::Symbol
  optval::Float64OrNothing

  function Problem(head::Symbol, objective::AbstractExpr, constraints::Array=Constraint[])
    return new(head, objective, constraints, "not yet solved", nothing)
  end
end


function find_variable_ranges(constraints)
  index = 0
  constr_size = 0
  var_to_ranges = Dict{Uint64, (Int64, Int64)}()
  for constraint in constraints
    constr_size += constraint.size
    for (id, val) in constraint.vars_to_coeffs
      if !haskey(var_to_ranges, id) && id != object_id(:constant)
        var = id_to_variables[id]
        var_to_ranges[id] = (index + 1, index + get_vectorized_size(var))
        index += get_vectorized_size(var)
      end
    end
  end
  return index, constr_size, var_to_ranges
end

function dual_conic_problem(p::Problem)
  objective_var = Variable()
  constraints = dual_conic_form(p.objective - objective_var == 0)[2]
  append!(constraints, dual_conic_form(p.objective)[2])
  for constraint in p.constraints
    append!(constraints, dual_conic_form(constraint)[2])
  end
  var_size, constr_size, var_to_ranges = find_variable_ranges(constraints)

  c = spzeros(var_size, 1)
  objective_range = var_to_ranges[objective_var.id]
  c[objective_range[1]:objective_range[2]] = 1

  A = spzeros(constr_size, var_size)
  b = spzeros(constr_size, 1)
  cones = (Symbol, UnitRange{Int64})[]
  constr_index = 0
  for constraint in constraints
    for (id, val) in constraint.vars_to_coeffs
      if id == object_id(:constant)
        b[constr_index + 1 : constr_index + constraint.size] = val
      else
        var_range = var_to_ranges[id]
        A[constr_index + 1 : constr_index + constraint.size, var_range[1] : var_range[2]] = -val
      end
    end
    push!(cones, (constraint.cone, constr_index + 1 : constr_index + constraint.size))
    constr_index += constraint.size
  end
  return c, A, b, cones
end
