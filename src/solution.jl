import MathProgBase, SCS, ECOS
export solve!

function solve!(problem::Problem, m::MathProgBase.AbstractMathProgModel=SCS.SCSMathProgModel())

  c, A, b, cones = conic_problem(problem)

  if problem.head == :maximize
    c = -c
  end

  # TODO: Fix once MathProgBase has a loadineqproblem!
  # TODO: Get rid of full once c and b are not sparse
  if typeof(m) == ECOS.ECOSMathProgModel
    ECOS.loadineqconicproblem!(m, full(c), A, full(b), cones)
    ECOS.optimize!(m)
  else
    SCS.loadineqconicproblem!(m, full(c), A, full(b), cones)
    SCS.optimize!(m)
  end

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
end
