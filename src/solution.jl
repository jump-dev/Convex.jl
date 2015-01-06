import MathProgBase
export solve!

function solve!(problem::Problem,
                s::MathProgBase.AbstractMathProgSolver=get_default_solver();
                warmstart=true)

  m = MathProgBase.model(s)
  # TODO: This is a tiny, temporary hack that will be removed once SCS.jl or SCS
  # starts to take care of symmetry constraints.
  old_solver = get_default_solver()
  set_default_solver(s)
  c, A, b, cones, var_to_ranges, vartypes = conic_problem(problem)
  set_default_solver(old_solver)

  if problem.head == :maximize
    c = -c
  end

  # no conic constraints on variables
  var_cones = fill((:Free, 1:size(A, 2)),1)
  # TODO: Get rid of full once c and b are not sparse
  MathProgBase.loadconicproblem!(m, vec(full(c)), A, vec(full(b)), cones, var_cones)

  # add integer and binary constraints on variables
  if !all(Bool[t==:Cont for t in vartypes])
    try
      MathProgBase.setvartype!(m, vartypes)
    catch
      error("model $(typeof(m)) does not support variables of some of the following types: $(unique(vartypes))")
    end
  end

  # see if we should warmstart (as of 11/19/14, only MILP solvers support this)
  if warmstart && problem.status==:Optimal
    try
      MathProgBase.setwarmstart!(m, problem.solution.primal)
      println("Using warm start from previous solution")
    end
  end

  # optimize problem
  status = MathProgBase.optimize!(m)

  # get the primal (and possibly dual) solution
  try
    dual = MathProgBase.getconicdual(m)
    problem.solution = Solution(MathProgBase.getsolution(m), dual,
                                MathProgBase.status(m), MathProgBase.getobjval(m))
  catch
    problem.solution = Solution(MathProgBase.getsolution(m),
                                MathProgBase.status(m), MathProgBase.getobjval(m))
  end
  populate_variables!(problem, var_to_ranges)

  # minimize -> maximize
  if (problem.head == :maximize)
    problem.solution.optval = -problem.solution.optval
  end

  # Populate the problem with the solution
  problem.optval = problem.solution.optval
  problem.status = problem.solution.status

  if !(problem.status==:Optimal)
    warn("Problem status $(problem.status); solution may be inaccurate.")
  end
end

function populate_variables!(problem::Problem, var_to_ranges::Dict{Uint64, (Int64, Int64)})
  x = problem.solution.primal
  for (id, (start_index, end_index)) in var_to_ranges
    var = id_to_variables[id]
    sz = var.size
    var.value = reshape(x[start_index:end_index], sz[1], sz[2])
    if sz == (1, 1)
      var.value = var.value[1]
    end
  end
end
