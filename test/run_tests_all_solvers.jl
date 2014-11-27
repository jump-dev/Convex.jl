using Convex
solvers = Any[]

if isdir(Pkg.dir("ECOS"))
    using ECOS
    push!(solvers, ECOSSolver)
end

if isdir(Pkg.dir("SCS"))
    using SCS
    push!(solvers, SCSSolver)
end

if isdir(Pkg.dir("Gurobi"))
    using Gurobi
    push!(solvers, GurobiSolver)
end

if isdir(Pkg.dir("Mosek"))
    using Mosek
    push!(solvers, MosekSolver)
end

for solver in solvers
    set_default_solver(solver)
    println("Running tests with $(solver):")
    include("run_tests.jl")
end
