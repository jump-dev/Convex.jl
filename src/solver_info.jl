using MathProgBase
export can_solve_mip, can_solve_sdp
export set_default_solver, get_default_solver

function set_default_solver(solver::MathProgBase.MathProgSolverInterface.AbstractMathProgSolver)
  global DEFAULT_SOLVER
  DEFAULT_SOLVER = solver
end

function get_default_solver()
  return DEFAULT_SOLVER
end

if isdir(Pkg.dir("ECOS"))
  using ECOS
  set_default_solver(ECOSSolver())
end

if isdir(Pkg.dir("SCS")) && DEFAULT_SOLVER == nothing
  using SCS
  set_default_solver(SCSSolver())
end
if isdir(Pkg.dir("Gurobi")) && DEFAULT_SOLVER == nothing
  using Gurobi
  set_default_solver(GurobiSolver())
end
if isdir(Pkg.dir("Mosek")) && DEFAULT_SOLVER == nothing
  using Mosek
  set_default_solver(MosekSolver())
end

if get_default_solver() == nothing
  error("You have any of ECOS.jl, SCS.jl, Mosek.jl or Gurobi.jl installed. Must have at least one of these solvers.")
end

function can_solve_mip(solver)
  name = typeof(solver).name.name
  if name == :GurobiSolver || name == :MosekSolver || name == :GLPKSolverMIP
    return true
  else
    info("Only GurobiSolver, MosekSolver and GLPKSolverMIP can solve mixed integer programs")
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
  if name == :SCSSolver || name == :SCSMathProgModel #|| name == :MosekSolver
    return true
  else
    info("Only SCSSolver and SCSMathProgModel can solve semidefinite programs")
    # info("Only SCSSolver, MosekSolver and SCSMathProgModel can solve semidefinite programs")
    return false
  end
end

