import MathProgBase, ECOS, SCS
export solve!

function solve!(problem::Problem, m::MathProgBase.AbstractMathProgModel=ECOS.ECOSMathProgModel())

  c, A, b, cones, var_to_ranges = conic_problem(problem)

  if problem.head == :maximize
    c = -c
  end

  # TODO: Fix once MathProgBase has a loadineqproblem!
  # TODO: Get rid of full once c and b are not sparse
  if typeof(m) == ECOS.ECOSMathProgModel
    ECOS.loadineqconicproblem!(m, full(c), A, full(b), cones)
    ECOS.optimize!(m)
  elseif typeof(m) == SCS.SCSMathProgModel
    SCS.loadineqconicproblem!(m, full(c), A, full(b), cones)
    SCS.optimize!(m)
  else
    error("model type $(typeof(m)) not recognized")
  end

  try
    y, z = MathProgBase.getconicdual(m)
    problem.solution = Solution(MathProgBase.getsolution(m), y, z, MathProgBase.status(m), MathProgBase.getobjval(m))
  catch
    problem.solution = Solution(MathProgBase.getsolution(m), MathProgBase.status(m), MathProgBase.getobjval(m))
    if typeof(m) == ECOS.ECOSMathProgModel
      populate_variables!(problem, var_to_ranges)
    end
  end
  # minimize -> maximize
  if (problem.head == :maximize) && (problem.solution.status == :Optimal)
    problem.solution.optval = -problem.solution.optval
  end

  # Populate the problem with the solution
  problem.optval = problem.solution.optval
  problem.status = problem.solution.status
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
