using Pkg
import MathProgBase
export can_solve_mip, can_solve_socp, can_solve_sdp, can_solve_exp
export isinstalled

function isinstalled(pkg)
    for path in Base.DEPOT_PATH
        if isdir(joinpath(path, pkg)) || isdir(joinpath(path, "packages", pkg))
            return true
        end
    end

    return false
end

# TODO: I have not listed solvers such as CPLEX etc because I have not tested Convex with them
solvers = [("ECOS", "ECOSSolver"), ("SCS", "SCSSolver"), ("Gurobi", "GurobiSolver"), ("Mosek", "MosekSolver"),
           ("GLPKMathProgInterface", "GLPKSolverMIP")]


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
