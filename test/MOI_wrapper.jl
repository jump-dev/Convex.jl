module TestMOIWrapper

using Test
import MathOptInterface as MOI
import Convex, ECOS

function test_runtests()
    optimizer = MOI.Utilities.CachingOptimizer(
        MOI.Utilities.UniversalFallback(MOI.Utilities.Model{Float64}()),
        MOI.instantiate(
            () -> Convex.Optimizer(ECOS.Optimizer()),
            with_bridge_type = Float64,
        ),
    )
    config = MOI.Test.Config(; exclude = Any[])
    return MOI.Test.runtests(
        optimizer,
        config;
        exclude = String[
        # GLPK issue
        # glp_simplex: column 1: lb = 1, ub = 0; incorrect bounds
        # glp_intopt: column 1: lb = 1, ub = 0; incorrect bounds
        # test_constraint_ZeroOne_bounds_3: Test Failed at /home/blegat/.julia/packages/MathOptInterface/pgWRA/src/Test/test_constraint.jl:631
        #  Expression: MOI.get(model, MOI.TerminationStatus()) == config.infeasible_status
        #   Evaluated: MathOptInterface.INVALID_MODEL == MathOptInterface.INFEASIBLE
            "test_constraint_ZeroOne_bounds_3",],
    )
end

function runtests()
    for name in names(@__MODULE__; all = true)
        if startswith("$name", "test_")
            @testset "$(name)" begin
                getfield(@__MODULE__, name)()
            end
        end
    end
end

end

TestMOIWrapper.runtests()
