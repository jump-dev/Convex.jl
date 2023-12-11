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
    config = MOI.Test.Config(;
        atol = 1e-4,
        rtol = 1e-4,
        exclude = Any[
            MOI.ConstraintBasisStatus,
            MOI.VariableBasisStatus,
            MOI.ObjectiveBound,
        ],
    )
    return MOI.Test.runtests(
        optimizer,
        config;
        exclude = String[],
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
