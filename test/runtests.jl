using Convex
using Test

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

solvers = Any[]

if isinstalled("ECOS")
    using ECOS
    push!(solvers, ECOSSolver(verbose=0))
end

if isinstalled("SCS")
    using SCS
    push!(solvers, SCSSolver(verbose=0, eps=1e-5))
end

if isinstalled("Gurobi")
    using Gurobi
    push!(solvers, GurobiSolver(OutputFlag=0))
end

if isinstalled("Mosek")
    using Mosek
    push!(solvers, MosekSolver(LOG=0))
end

if isinstalled("GLPK") && isinstalled("GLPKMathProgInterface")
    using GLPKMathProgInterface
    push!(solvers, GLPKSolverMIP())
end


for solver in solvers
    println("Running tests with $(solver):")
    set_default_solver(solver)
    println(get_default_solver())
    include("runtests_single_solver.jl")
end

