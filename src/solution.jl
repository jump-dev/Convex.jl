import MathProgBase
export solve!

SolverOrModel = Union(MathProgBase.AbstractMathProgSolver, MathProgBase.AbstractMathProgModel, Nothing)

function solve!(problem::Problem,
                s::SolverOrModel=get_default_solver();
                warmstart=true, check_vexity=false)

  if isa(s, MathProgBase.AbstractMathProgSolver)
    m = MathProgBase.model(s)
  elseif s != nothing # it is already a model
    warn("deprecated syntax. Use AbstractMathProgSolver instead. eg ECOSSolver() or SCSSolver()")
    m = s
  end

  if s == nothing
    error("The default solver is set to `nothing`
         You must have at least one solver installed.
         You can install a solver such as SCS by running:
         Pkg.add(\"SCS\").
         You will have to restart Julia after that.")
  end

  if check_vexity
    vex = vexity(problem)
  end

  c, A, b, cones, var_to_ranges, vartypes, conic_constraints = conic_problem(problem)

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

  if problem.solution.has_dual
    populate_duals!(conic_constraints, problem.solution.dual)
  end

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

function populate_variables!(problem::Problem, var_to_ranges::Dict{Uint64, (Int, Int)})
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

function populate_duals!{T}(constraints::Array{ConicConstr}, dual::Array{T, 1})
  constr_index = 1
  for constraint in constraints
    # conic_constr_to_constr only has keys for conic constraints with a single objective
    # so this will work
    if haskey(conic_constr_to_constr, constraint)
      sz = constraint.sizes[1]
      c = conic_constr_to_constr[constraint]
      c.dual = reshape(dual[constr_index:constr_index+sz-1], c.size)
      if c.size == (1, 1)
        c.dual = c.dual[1]
      end
      constr_index += sz
    else
      for i = 1:length(constraint.objs)
        sz = constraint.sizes[i]
        constr_index += sz
      end
    end
  end
end
