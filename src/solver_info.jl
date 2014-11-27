export can_solve_mip, can_solve_sdp, DEFAULT_SOLVER
export set_default_solver, get_default_solver

global DEFAULT_SOLVER = nothing

if isdir(Pkg.dir("ECOS"))
  using ECOS
  DEFAULT_SOLVER = ECOSSolver
end

if isdir(Pkg.dir("SCS")) && DEFAULT_SOLVER == nothing
  using SCS
  DEFAULT_SOLVER = SCS.SCSMathProgModel
end
if isdir(Pkg.dir("Gurobi")) && DEFAULT_SOLVER == nothing
  using Gurobi
  DEFAULT_SOLVER = GurobiSolver
end
if isdir(Pkg.dir("Mosek")) && DEFAULT_SOLVER == nothing
  using Mosek
  DEFAULT_SOLVER = MosekSolver
end

if DEFAULT_SOLVER == nothing
  error("You have any of ECOS.jl, SCS.jl, Mosek.jl or Gurobi.jl installed. Must have at least one of these solvers.")
end

function set_default_solver(solver)
  global DEFAULT_SOLVER
  DEFAULT_SOLVER = solver
end

function get_default_solver()
  return DEFAULT_SOLVER
end

function can_solve_mip(solver)
  name = solver.name.name
  if name == :GurobiSolver || name == :MosekSolver || name == :GLPKSolverMIP
    return true
  else
    info("Only GurobiSolver, MosekSolver and GLPKSolverMIP can solve mixed integer programs")
    return false
  end
end

function can_solve_exp(solver)
  name = solver.name.name
  if name == :SCSSolver || name == :SCSMathProgModel #|| name == :MosekSolver
    return true
  else
    info("Only SCSSolver and SCSMathProgModel can solve exponential programs")
    # info("Only SCSSolver, MosekSolver and SCSMathProgModel can solve exponential programs")
    return false
  end
end

function can_solve_sdp(solver)
  name = solver.name.name
  if name == :SCSSolver || name == :SCSMathProgModel #|| name == :MosekSolver
    return true
  else
    info("Only SCSSolver and SCSMathProgModel can solve semidefinite programs")
    # info("Only SCSSolver, MosekSolver and SCSMathProgModel can solve semidefinite programs")
    return false
  end
end

