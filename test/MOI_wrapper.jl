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
    MOI.Bridges.remove_bridge(
        optimizer.optimizer.context.model,
        MOI.Bridges.Variable.ZerosBridge{T},
    )
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

function test_scalar_nonlinear_function()
    model = Convex.Optimizer(SCS.Optimizer)
    x = MOI.add_variable(model)
    y = Convex._expr(model, x)
    for (f, g) in Any[
        3.0 => 3.0,
        x => y,
        1.0 * x => 1.0 * y,
        2.0 * x + 3.0 => 2.0 * y + 3.0,
        1.0 * x * x => 1.0 * y * y,
        1.0 * x * x + 2.0 * x + 3.0 => 1.0 * y * y + 2.0 * y + 3.0,
        MOI.ScalarNonlinearFunction(:+, Any[x]) => y,
        MOI.ScalarNonlinearFunction(:+, Any[x, 2]) => +(y, 2),
        MOI.ScalarNonlinearFunction(:+, Any[x, x, x]) => sum([y, y, y]),
        MOI.ScalarNonlinearFunction(:-, Any[x]) => -y,
        MOI.ScalarNonlinearFunction(:-, Any[x, 2]) => -(y, 2),
        MOI.ScalarNonlinearFunction(:*, Any[x, 2]) => *(y, 2),
        MOI.ScalarNonlinearFunction(:/, Any[x, 2]) => y / 2,
        MOI.ScalarNonlinearFunction(:^, Any[x, 2]) => square(y),
        MOI.ScalarNonlinearFunction(:min, Any[x]) => y,
        MOI.ScalarNonlinearFunction(:min, Any[x, x]) => min(y, y),
        MOI.ScalarNonlinearFunction(:min, Any[x, x, x]) => minimum([y, y, y]),
        MOI.ScalarNonlinearFunction(:max, Any[x]) => y,
        MOI.ScalarNonlinearFunction(:max, Any[x, x]) => max(y, y),
        MOI.ScalarNonlinearFunction(:max, Any[x, x, x]) => maximum([y, y, y]),        MOI.ScalarNonlinearFunction(:abs, Any[x]) => abs(y),
        MOI.ScalarNonlinearFunction(:sqrt, Any[x]) => sqrt(y),
        MOI.ScalarNonlinearFunction(:exp, Any[x]) => exp(y),
        MOI.ScalarNonlinearFunction(:log, Any[x]) => log(y),
    ]
        @test sprint(show, Convex._expr(model, f)) == sprint(show, g)
    end
    for f in [
        MOI.ScalarNonlinearFunction(:foo, Any[x]),
        MOI.ScalarNonlinearFunction(:^, Any[2, x]),
        MOI.ScalarNonlinearFunction(:-, Any[x, x, x]),
        MOI.ScalarNonlinearFunction(:/, Any[1, x]),
    ]
        @test_throws MOI.UnsupportedNonlinearOperator Convex._expr(model, f)
    end
    return
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
