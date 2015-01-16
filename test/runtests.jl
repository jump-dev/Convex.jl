using Convex
using FactCheck

solvers = Any[]

if isdir(Pkg.dir("ECOS"))
    using ECOS
    push!(solvers, ECOSSolver(verbose=0))
end

if isdir(Pkg.dir("SCS"))
    using SCS
    push!(solvers, SCSSolver(verbose=0))
end

if isdir(Pkg.dir("Gurobi"))
    using Gurobi
    push!(solvers, GurobiSolver())
end

# if isdir(Pkg.dir("Mosek"))
#     using Mosek
#     push!(solvers, MosekSolver())
# end

if isdir(Pkg.dir("GLPK")) && isdir(Pkg.dir("GLPKMathProgInterface"))
    using GLPKMathProgInterface
    push!(solvers, GLPKSolverMIP())
end


for solver in solvers
    println("Running tests with $(solver):")
    set_default_solver(solver)
    println(get_default_solver())
    include("runtests_single_solver.jl")
end

FactCheck.exitstatus()
