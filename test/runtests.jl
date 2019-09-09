using Convex
using Convex.ProblemDepot: run_test
using Test
using SCS, ECOS, GLPKMathProgInterface


# Seed random number stream to improve test reliability
using Random
Random.seed!(2)

@testset "Convex" begin
    include("test_utilities.jl")

    @testset "SCS" begin
        run_test(; exclude=[r"mip"]) do p
            solve!(p, SCSSolver(verbose=0, eps=1e-6))
        end
    end

    @testset "ECOS" begin
        run_test(; exclude=[r"mip", r"sdp"]) do p
            solve!(p, ECOSSolver(verbose=0))
        end
    end

    @testset "GLPK MIP" begin
        run_test(; exclude=[r"socp", r"sdp", r"exp"]) do p
            solve!(p, GLPKSolverMIP())
        end
end
end
