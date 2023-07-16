using Convex
using Convex.ProblemDepot: run_tests
using Test
using SCS, ECOS, GLPK, Clarabel

using MathOptInterface
const MOI = MathOptInterface
const MOIU = MOI.Utilities
const MOIB = MOI.Bridges

# Seed random number stream to improve test reliability
using Random
Random.seed!(2)

@testset "Convex" begin
    @testset "ProblemDepot" begin
        @testset "Problems can run without `solve!`ing if `test==false`; T=$T" for T in
                                                                                   (
            Float64,
            # BigFloat,
        )
            Convex.ProblemDepot.foreach_problem() do name, func
                @testset "$name" begin
                    # We want to check to make sure this does not throw
                    func(Val(false), 0.0, 0.0, T) do problem
                        @test problem isa Convex.Problem{T} # check numeric type
                        model = MOIU.MockOptimizer(MOIU.Model{T}())

                        # make sure it loads without throwing
                        context = Convex.Context(problem, model)

                        return nothing
                    end
                end
            end
        end
    end

    # @testset "SCS with warmstarts" begin
    #     run_tests(
    #         exclude = [
    #             r"mip",
    #             # TODO(odow): investigate
    #             r"sdp_lieb_ando",
    #             # Tolerance issue with SCS 3.0
    #             r"sdp_Real_Variables_with_complex_equality_constraints",
    #         ],
    #     ) do p
    #         return solve!(
    #             p,
    #             MOI.OptimizerWithAttributes(
    #                 SCS.Optimizer,
    #                 "verbose" => 0,
    #                 "eps_rel" => 1e-6,
    #                 "eps_abs" => 1e-6,
    #             );
    #             warmstart = true,
    #         )
    #     end
    # end

    @testset "Clarabel" begin
        # TODO(odow): investigate sdp_lieb_ando
        run_tests(; exclude = [r"mip", r"sdp_lieb_ando"]) do p
            return solve!(p, Clarabel.Optimizer; silent_solver = true)
        end
    end

    # Slow, and also hits linear algebra issues with BigFloat matrices
    # @testset "Clarabel with BigFloat" begin
    #     # TODO(odow): investigate sdp_lieb_ando
    #     run_tests(; exclude = [r"mip", r"sdp_lieb_ando"], T = BigFloat) do p
    #         return solve!(p, Clarabel.Optimizer{BigFloat}; silent_solver = true)
    #     end
    # end

    # Disabling ECOS due to this segfault:
    # https://github.com/jump-dev/ECOS.jl/issues/144
    # @testset "ECOS" begin
    #     # For `rational_norm` problems:
    #     # >   MathOptInterface.UnsupportedConstraint{MathOptInterface.VectorAffineFunction{Float64}, MathOptInterface.NormCone}: `MathOptInterface.VectorAffineFunction{Float64}`-in-`MathOptInterface.NormCone` constraint is not supported by the model.
    #     # Missing a bridge?
    #     run_tests(; exclude = [r"mip", r"sdp", r"rational_norm"]) do p
    #         return solve!(p, ECOS.Optimizer; silent_solver = true)
    #     end
    # end

    @testset "GLPK" begin
        # Note this is an inclusion, not exclusion;
        # we only test GLPK with MIPs, since those we can't test elsewhere.
        run_tests([r"mip"]) do p
            return solve!(p, GLPK.Optimizer; silent_solver = true)
        end
    end

    include("definitions.jl")
    include("SparseTape.jl")
    include("test_utilities.jl")
    include("test_abstract_variable.jl")
end
