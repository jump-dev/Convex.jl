module TestMOIWrapper

using Test
import MathOptInterface as MOI
import Convex, ECOS

function test_runtests()
    T = Float64
    optimizer = MOI.Utilities.CachingOptimizer(
        MOI.Utilities.UniversalFallback(MOI.Utilities.Model{T}()),
        Convex.Optimizer(ECOS.Optimizer),
    )
    MOI.set(optimizer, MOI.Silent(), true)
    MOI.Bridges.remove_bridge(optimizer.optimizer.context.model, MOI.Bridges.Variable.ZerosBridge{T})
    config = MOI.Test.Config(;
        atol = 1e-3,
        rtol = 1e-3,
        exclude = Any[
            MOI.ConstraintBasisStatus,
            MOI.VariableBasisStatus,
            MOI.ObjectiveBound,
        ],
    )
    return MOI.Test.runtests(optimizer, config)
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
