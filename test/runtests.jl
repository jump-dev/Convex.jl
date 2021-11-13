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
    @testset "Problems can run without `solve!`ing if `test==false`; T=$T" for T in (Float64, BigFloat)
        Convex.ProblemDepot.foreach_problem() do name, func
            @testset "$name" begin
                # We want to check to make sure this does not throw
                func(Val(false), 0.0, 0.0, T) do problem
                    @test problem isa Convex.Problem{T} # check numeric type
                    model = MOIU.MockOptimizer(MOIU.Model{T}())
                    Convex.load_MOI_model!(model, problem) # make sure it loads without throwing
                end
            end
        end
    end
end

@testset "Definitions" begin
    include("definitions.jl")
end

@testset "Convex" begin
    include("test_utilities.jl")
    include("test_abstract_variable.jl")

    @testset "SCS with warmstarts" begin
        run_tests(; exclude=[r"mip"]) do p
            solve!(p, () -> SCS.Optimizer(verbose=0, eps=1e-6); warmstart = true)
        end
    end

    @testset "SCS" begin
        run_tests(; exclude=[r"mip"]) do p
            solve!(p, () -> SCS.Optimizer(verbose=0, eps=1e-6))
        end
    end

    @testset "ECOS" begin
        run_tests(; exclude=[r"mip", r"sdp"]) do p
            solve!(p, ECOS.Optimizer; silent_solver=true)
        end
    end

    @testset "GLPK" begin
        # Note this is an inclusion, not exclusion;
        # we only test GLPK with MIPs, since those we can't test elsewhere.
        run_tests([r"mip"]) do p
            solve!(p, GLPK.Optimizer; silent_solver=true)
        end
    end
end
