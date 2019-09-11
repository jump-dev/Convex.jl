using Convex
using Convex.ProblemDepot: run_tests
using Test
using SCS, ECOS, GLPK

using MathOptInterface
const MOI = MathOptInterface
const MOIU = MOI.Utilities

# Seed random number stream to improve test reliability
using Random
Random.seed!(2)

@testset "ProblemDepot" begin
    @testset "Problems can run without `solve!`ing if `test==false`" begin
        Convex.ProblemDepot.foreach_problem() do name, func
            @testset "$name" begin
                # We want to check to make sure this does not throw
                func(Val(false), 0.0, 0.0, Float64) do problem
                    model = MOIU.MockOptimizer(MOIU.Model{Float64}())
                    Convex.load_MOI_model!(model, problem)
                end
                @test true
            end
        end
    end
end

@testset "Convex with MOI" begin
    include("test_moi.jl")
end

@testset "Convex" begin
    include("test_utilities.jl")

    @testset "SCS" begin
        run_tests(; exclude=[   r"mip",
                                r"sdp_matrix_frac_atom", # bug: https://github.com/JuliaOpt/SCS.jl/issues/153
                            ]) do p
            solve!(p, SCS.Optimizer(verbose=0, eps=1e-6))
        end
    end

    @testset "ECOS" begin
        run_tests(; exclude=[r"mip", r"sdp"]) do p
            solve!(p, ECOS.Optimizer(verbose=0))
        end
    end

    @testset "GLPK" begin
        run_tests(; exclude=[r"exp", r"sdp", r"socp"]) do p
            solve!(p, GLPK.Optimizer(msg_lev = GLPK.OFF))
        end
    end
end
