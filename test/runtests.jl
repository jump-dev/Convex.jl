using Convex
using Convex: DotMultiplyAtom
using Test
using ECOS
using SCS
using GLPKMathProgInterface
using Random

import LinearAlgebra.eigen
import LinearAlgebra.I
import LinearAlgebra.opnorm
import Random.shuffle
import Statistics.mean

TOL = 1e-3
eye(n) = Matrix(1.0I, n, n)

# Seed random number stream to improve test reliability
Random.seed!(2)

solvers = Any[]

push!(solvers, ECOSSolver(verbose=0))
push!(solvers, GLPKSolverMIP())
push!(solvers, SCSSolver(verbose=0, eps=1e-6))

# If Gurobi is installed, uncomment to test with it:
#using Gurobi
#push!(solvers, GurobiSolver(OutputFlag=0))

# If Mosek is installed, uncomment to test with it:
#using Mosek
#push!(solvers, MosekSolver(LOG=0))

@testset "Convex" begin
    include("test_utilities.jl")
    include("test_constraints.jl")
    include("test_affine.jl")
    include("test_lp.jl")
    include("test_socp.jl")
    include("test_sdp.jl")
    include("test_exp.jl")
    include("test_sdp_and_exp.jl")
    include("test_mip.jl")
end
