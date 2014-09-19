import ECOS, MathProgBase

export get_var_dict, solve!

function solve!(problem::Problem, m::MathProgBase.AbstractMathProgModel=ECOS.ECOSMathProgModel())
  # maximize -> minimize
  if problem.head == :maximize
    problem.objective = -problem.objective
  end

  ecos_problem, variable_index, eq_constr_index, ineq_constr_index = ECOSConicProblem(problem)
  cp = IneqConicProblem(ecos_problem)
  ECOS.loadineqconicproblem!(m, cp.c, cp.A, cp.b, cp.G, cp.h, cp.cones)
  ECOS.optimize!(m)
  try
    y, z = MathProgBase.getconicdual(m)
    problem.solution = Solution(MathProgBase.getsolution(m), y, z, MathProgBase.status(m), MathProgBase.getobjval(m))
  catch
    problem.solution = Solution(MathProgBase.getsolution(m), MathProgBase.status(m), MathProgBase.getobjval(m))
  end
  # minimize -> maximize
  if (problem.head == :maximize) && (problem.solution.status == :Optimal)
    problem.solution.optval = -problem.solution.optval
  end

  # Populate the problem with the solution
  problem.optval = problem.solution.optval
  problem.status = problem.solution.status

  if problem.status == :Optimal
    populate_variables!(problem, variable_index)
    if problem.solution.dual
      populate_constraints!(problem, eq_constr_index, ineq_constr_index)
    end
  end
end

function canonical_constraints(problem::Problem)
  # need to change objective if problem.head == :maximize?
    canonical_constraints_array = CanonicalConstr[]
    for constraint in problem.constraints
        append!(canonical_constraints_array, constraint.canon_form())
    end
    append!(canonical_constraints_array, problem.objective.canon_form())
    return canonical_constraints_array
end

# Now that the problem has been solved, populate the optimal values of the
# variables back into them
function populate_variables!(problem::Problem, variable_index::Dict{Int64, Int64})
  x = problem.solution.x
  var_dict = get_var_dict(problem.objective, problem.constraints)
  for (id, var) in var_dict
    index = variable_index[id]
    var.value = Base.reshape(x[index : index + get_vectorized_size(var) - 1], var.size)
    if var.size == (1, 1)
      # Make it a scalar
      var.value = var.value[1]
    end
  end
end

function ECOSConicProblem(problem::Problem)
    canonical_constraints_array = canonical_constraints(problem)
    m, n, p, l, ncones, q, G, h, A, b, variable_index, eq_constr_index, ineq_constr_index =
        create_ecos_matrices(canonical_constraints_array, problem.objective)

    # Now, all we need to do is create c
    c = zeros(n, 1)
    objective = problem.objective
    if objective.vexity != :constant
        uid = objective.uid
        c[variable_index[uid] : variable_index[uid] + objective.size[1] - 1] = 1
    end
    return ECOSConicProblem(n=n, m=m, p=p, l=l, ncones=ncones, q=q, G=G, c=c, h=h, A=A, b=b), variable_index, eq_constr_index, ineq_constr_index
end
function populate_constraints!(problem::Problem, eq_constr_index::Dict{Int64, Int64},
    ineq_constr_index::Dict{Int64, Int64})

  y = problem.solution.y
  z = problem.solution.z

  for constraint in problem.constraints
    uid = constraint.canon_uid
    if constraint.head == :(==)
      index = eq_constr_index[uid]
      constraint.dual_value = Base.reshape(y[index : index + get_vectorized_size(constraint.size) - 1], constraint.size)
    else
      index = ineq_constr_index[uid]
      constraint.dual_value = Base.reshape(z[index : index + get_vectorized_size(constraint.size) - 1], constraint.size)
    end
  end
end

# Recursively traverses the AST for the AbstractCvxExpr and finds the variables
# that were defined
# Updates var_dict with the ids of the variables as keys and variables as values
function get_var_dict!(e::AbstractCvxExpr, var_dict::Dict{Int64, Variable})
  if e.head == :variable
    var_dict[e.uid] = e
  elseif e.head == :parameter || e.head == :constant
    return
  else
    for v in e.args
      get_var_dict!(v, var_dict)
    end
  end
end

function get_var_dict(p::Problem)
  return get_var_dict(p.objective, p.constraints)
end

function get_var_dict(objective::AbstractCvxExpr, constraints::Array{CvxConstr})
  var_dict = Dict{Int64, Variable}()

  get_var_dict!(objective, var_dict)
  for c in constraints
    get_var_dict!(c.lhs, var_dict);
    get_var_dict!(c.rhs, var_dict);
  end

  return var_dict
end
