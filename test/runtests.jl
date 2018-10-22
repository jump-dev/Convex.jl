using Convex
using Test
using ECOS
using SCS
using GLPKMathProgInterface

# Seed random number stream to improve test reliability
srand(2)

solvers = Any[]

push!(solvers, ECOSSolver(verbose=0))
push!(solvers, GLPKSolverMIP())
push!(solvers, SCSSolver(verbose=0, eps=1e-6))

if isinstalled("Gurobi")
    using Gurobi
    push!(solvers, GurobiSolver(OutputFlag=0))
end

if isinstalled("Mosek")
    using Mosek
    push!(solvers, MosekSolver(LOG=0))
end

for solver in solvers
    println("Running tests with $(solver):")
    set_default_solver(solver)
    println(get_default_solver())
    include("runtests_single_solver.jl")
end
