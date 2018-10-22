using Pkg
import MathProgBase
export can_solve_mip, can_solve_socp, can_solve_sdp, can_solve_exp
export set_default_solver, get_default_solver

function set_default_solver(solver::MathProgBase.AbstractMathProgSolver)
  global DEFAULT_SOLVER
  DEFAULT_SOLVER = solver
end

function get_default_solver()
  if DEFAULT_SOLVER == nothing
    error("The default solver is set to `nothing`
         You must have at least one solver installed to use Convex.
         You can install a solver such as SCS by running:
         Pkg.add(\"SCS\").
         You will have to restart Julia after that.")
  end
  return DEFAULT_SOLVER
end

# TODO: I have not listed solvers such as CPLEX etc because I have not tested Convex with them
solvers = [("SCS", "SCSSolver"), ("ECOS", "ECOSSolver"), ("Gurobi", "GurobiSolver"), ("Mosek", "MosekSolver"),
          ("GLPKMathProgInterface", "GLPKSolverMIP")]

function isinstalled(pkg)
    for path in Base.DEPOT_PATH
        if isdir(joinpath(path, pkg))
            return true
        elseif isdir(joinpath(path, "packages", pkg))
            return true
        end
    end
    return false
end

for (dir, solver) in solvers
  if isinstalled(dir) && DEFAULT_SOLVER == nothing
    eval(Meta.parse("using "*dir))
    eval(Meta.parse("set_default_solver("*solver*"())"))
  end
end


if get_default_solver() == nothing
  packages = ""
  for (dir, solver) in solvers
  global packages = packages*dir*" | "
  end
  @warn "***********************************************************************************************
       You don't have any of
       "*packages*" installed.
       You must have at least one of these solvers. You can install a solver such as SCS by running:
       Pkg.add(\"SCS\")
       You will have to restart Julia after that.
       ***********************************************************************************************"
end

function can_solve_mip(solver)
  name = typeof(solver).name.name
  if name == :GurobiSolver || name == :MosekSolver || name == :GLPKSolverMIP || name == :CPLEXSolver || name == :CbcSolver
    return true
  else
    @info "$name cannot solve mixed integer programs. Consider using Gurobi, Mosek, or GLPK."
    return false
  end
end

function can_solve_socp(solver)
  if :SOC in MathProgBase.supportedcones(solver)
    return true
  else
    name = typeof(solver).name.name
    @info "$name cannot solve second order cone programs. Consider using SCS, ECOS, Mosek, or Gurobi."
    return false
  end
end

function can_solve_exp(solver)
  if :ExpPrimal in MathProgBase.supportedcones(solver)
    return true
  else
    name = typeof(solver).name.name
    @info "$name cannot solve exponential programs. Consider using SCS or ECOS."
    return false
  end
end

function can_solve_sdp(solver)
  if :SDP in MathProgBase.supportedcones(solver)
    return true
  else
    name = typeof(solver).name.name
    @info "$name cannot solve semidefinite programs. Consider using SCS or Mosek."
    return false
  end
end

