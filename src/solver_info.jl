import MathProgBase
export can_solve_mip, can_solve_socp, can_solve_sdp, can_solve_exp

# NOTE: Convex has been tested with the following solvers:
#   ECOS:   ECOSSolver
#   SCS:    SCSSolver
#   Gurobi: GurobiSolver
#   Mosek:  MosekSolver
#   GLPKMathProgInterface: GLPKSolverMIP
# It has not been tested with other solvers such a CPLEX.

solver_name(solver) = Symbol(typeof(solver).name.module)
function can_solve_mip(solver)
    if solver_name(solver) ∈ (:Gurobi, :MosekTools, :GLPK, :CPLEX, :Cbc)
        return true
    else
        @info "$(solver_name(solver)) cannot solve mixed integer programs. Consider using Gurobi, Mosek, or GLPK."
        return false
    end
end

function can_solve_socp(solver)
    if solver_name(solver) ∈ (:SCS, :ECOS, :MosekTools, :Gurobi)
        return true
    else
        @info "$(solver_name(solver)) cannot solve second order cone programs. Consider using SCS, ECOS, Mosek, or Gurobi."
        return false
    end
end

function can_solve_exp(solver)
    if solver_name(solver) ∈ (:SCS, :ECOS, :MosekTools)
        return true
    else
        @info "$(solver_name(solver)) cannot solve exponential programs. Consider using SCS or ECOS."
        return false
    end
end

function can_solve_sdp(solver)
    if solver_name(solver) ∈ (:SCS, :MosekTools)
        return true
    else
        @info "$(solver_name(solver)) cannot solve semidefinite programs. Consider using SCS or Mosek."
        return false
    end
end
