import MathProgBase
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
    info("$name cannot solve mixed integer programs. Consider using Gurobi, Mosek, or GLPK.")
    return false
  end
end

function can_solve_socp(solver)
  if :SOC in MathProgBase.supportedcones(solver)
    return true
  else
    name = typeof(solver).name.name
    info("$name cannot solve second order cone programs. Consider using SCS, ECOS, Mosek, or Gurobi.")
    return false
  end
end

function can_solve_exp(solver)
  if :ExpPrimal in MathProgBase.supportedcones(solver)
    return true
  else
    name = typeof(solver).name.name
    info("$name cannot solve exponential programs. Consider using SCS or ECOS.")
    return false
  end
end

function can_solve_sdp(solver)
  if :SDP in MathProgBase.supportedcones(solver)
    return true
  else
    name = typeof(solver).name.name
    info("$name cannot solve semidefinite programs. Consider using SCS or Mosek.")
    return false
  end
end

