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
        # "sdp_lieb_ando" is currently (14 Nov 2021) failing with SCS
        # on ubuntu in CI (they pass locally on MacOS and in CI with
        # MacOS, and have passed on ubuntu in the past). Disabling
        # them for now; once COSMO or Hypatia is on MOI v0.10, we can
        # try using them, or hope SCS starts solving them again.
        #
        # "sdp_sdp_constraints" is failing with MOI v0.9 on ubuntu.
        run_tests(; exclude=[r"mip", r"sdp_lieb_ando", r"sdp_sdp_constraints"]) do p
            solve!(p, MOI.OptimizerWithAttributes(SCS.Optimizer, "verbose" => 0, "eps_abs" => 1e-6); warmstart = true)
        end
    end

    @testset "SCS" begin
        run_tests(; exclude=[r"mip", r"sdp_lieb_ando", r"sdp_sdp_constraints"]) do p
            solve!(p, MOI.OptimizerWithAttributes(SCS.Optimizer, "verbose" => 0, "eps_abs" => 1e-6))
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
