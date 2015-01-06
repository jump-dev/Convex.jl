using Convex
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

if isdir(Pkg.dir("Mosek"))
    using Mosek
    push!(solvers, MosekSolver())
end

for solver in solvers
    println("Running tests with $(solver):")
    set_default_solver(solver)
    println(get_default_solver())
    include("run_tests.jl")
end
