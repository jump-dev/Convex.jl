using MathProgBase
export can_solve_mip, can_solve_socp, can_solve_sdp, can_solve_exp
export set_default_solver, get_default_solver

function set_default_solver(solver::MathProgBase.MathProgSolverInterface.AbstractMathProgSolver)
  global DEFAULT_SOLVER
  DEFAULT_SOLVER = solver
end

function get_default_solver()
  if DEFAULT_SOLVER == nothing
    warn("No default solver currently set.")
  end
  return DEFAULT_SOLVER
end

# TODO: I have not listed solvers such as CPLEX etc because I have not tested Convex with them
solvers = [("SCS", "SCSSolver"), ("ECOS", "ECOSSolver"), ("Gurobi", "GurobiSolver"), ("Mosek", "MosekSolver"),
          ("GLPKMathProgInterface", "GLPKSolverMIP")]

for (dir, solver) in solvers
  if isdir(Pkg.dir(dir)) && DEFAULT_SOLVER == nothing
    eval(parse("using "*dir))
    eval(parse("set_default_solver("*solver*"())"))
  end
end


if get_default_solver() == nothing
  packages = ""
  for (dir, solver) in solvers
    packages = packages*dir*" | "
  end
  warn("***********************************************************************************************
       You don't have any of
       "*packages*" installed.
       You must have at least one of these solvers. You can install a solver such as SCS by running:
       Pkg.add(\"SCS\")
       You will have to restart Julia after that.
       ***********************************************************************************************")
end

function can_solve_mip(solver)
  name = typeof(solver).name.name
  if name == :GurobiSolver || name == :MosekSolver || name == :GLPKSolverMIP || name == :CPLEXSolver || name == :CbcSolver
    return true
  else
    info("Only GurobiSolver, MosekSolver and GLPKSolverMIP can solve mixed integer programs")
    return false
  end
end

function can_solve_socp(solver)
  name = typeof(solver).name.name
  if name == :ECOSSolver || name == :SCSSolver || name == :SCSMathProgModel || name == :MosekSolver || name == :GurobiSolver
    return true
  else
    info("Only ECOSSolver, SCSSolver, MosekSolver, GurobiSolver and SCSMathProgModel can solve second order cone programs")
    return false
  end
end

function can_solve_exp(solver)
  name = typeof(solver).name.name
  if name == :SCSSolver || name == :SCSMathProgModel #|| name == :MosekSolver
    return true
  else
    info("Only SCSSolver and SCSMathProgModel can solve exponential programs")
    # info("Only SCSSolver, MosekSolver and SCSMathProgModel can solve exponential programs")
    return false
  end
end

function can_solve_sdp(solver)
  name = typeof(solver).name.name
  if name == :SCSSolver || name == :SCSMathProgModel || name == :MosekSolver
    return true
  else
    info("Only SCSSolver, MosekSolver and SCSMathProgModel can solve semidefinite programs")
    return false
  end
end

