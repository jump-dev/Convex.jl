# Copyright (c) 2014: Madeleine Udell and contributors
#
# Use of this source code is governed by a BSD-style license that can be found
# in the LICENSE file or at https://opensource.org/license/bsd-2-clause

module TestMOIWrapper

using Test

import Convex
import ECOS
import JuMP
import LinearAlgebra
import MathOptInterface as MOI
import SCS
import SparseArrays

function runtests()
    for name in names(@__MODULE__; all = true)
        if startswith("$name", "test_")
            @testset "$(name)" begin
                getfield(@__MODULE__, name)()
            end
        end
    end
    return
end

function test_runtests()
    inner = Convex.Optimizer(ECOS.Optimizer)
    model = MOI.Utilities.CachingOptimizer(
        MOI.Utilities.UniversalFallback(MOI.Utilities.Model{Float64}()),
        MOI.Bridges.full_bridge_optimizer(inner, Float64),
    )
    MOI.Bridges.remove_bridge(
        inner.context.model.optimizer,
        MOI.Bridges.Variable.ZerosBridge{Float64},
    )
    MOI.set(model, MOI.Silent(), true)
    MOI.Test.runtests(
        model,
        MOI.Test.Config(;
            atol = 1e-3,
            rtol = 1e-3,
            exclude = Any[
                MOI.ConstraintBasisStatus,
                MOI.VariableBasisStatus,
                MOI.ObjectiveBound,
            ],
        ),
    )
    return
end

function test_issue_564()
    model = MOI.instantiate(
        () -> Convex.Optimizer(ECOS.Optimizer);
        with_bridge_type = Float64,
    )
    MOI.Test.test_add_parameter(model, MOI.Test.Config())
    return
end

function test_constraint_bridging_cost_vector_nonlinear()
    model = Convex.Optimizer(ECOS.Optimizer)
    for S in (
        MOI.Nonnegatives,
        MOI.Nonpositives,
        MOI.Zeros,
        MOI.SecondOrderCone,
        MOI.PositiveSemidefiniteConeTriangle,
        MOI.ExponentialCone,
    )
        @test MOI.supports_constraint(model, MOI.VectorNonlinearFunction, S)
        @test MOI.get(
            model,
            MOI.ConstraintBridgingCost{MOI.VectorNonlinearFunction,S}(),
        ) == 0.0
    end
    return
end

function test_scalar_nonlinear_function()
    model = Convex.Optimizer(ECOS.Optimizer)
    x = MOI.add_variable(model)
    y = Convex._expr(model, x)
    tested_operators = Set{Symbol}()
    for (f, g) in Any[
        3.0=>3.0,
        x=>y,
        1.0*x=>1.0*y,
        2.0*x+3.0=>2.0*y+3.0,
        1.0*x*x=>1.0*y*y,
        1.0*x*x+2.0*x+3.0=>1.0*y*y+2.0*y+3.0,
        MOI.ScalarNonlinearFunction(:+, Any[x])=>y,
        MOI.ScalarNonlinearFunction(:+, Any[x, 2])=>+(y, 2),
        MOI.ScalarNonlinearFunction(:+, Any[x, x, x])=>sum([y, y, y]),
        MOI.ScalarNonlinearFunction(:-, Any[x])=>-y,
        MOI.ScalarNonlinearFunction(:-, Any[x, 2])=>-(y, 2),
        MOI.ScalarNonlinearFunction(:*, Any[x, 2])=>*(y, 2),
        MOI.ScalarNonlinearFunction(:/, Any[x, 2])=>y/2,
        MOI.ScalarNonlinearFunction(:^, Any[x, 2])=>Convex.square(y),
        MOI.ScalarNonlinearFunction(:min, Any[x])=>y,
        MOI.ScalarNonlinearFunction(:min, Any[x, x])=>min(y, y),
        MOI.ScalarNonlinearFunction(:min, Any[x, x, x])=>minimum([y, y, y]),
        MOI.ScalarNonlinearFunction(:max, Any[x])=>y,
        MOI.ScalarNonlinearFunction(:max, Any[x, x])=>max(y, y),
        MOI.ScalarNonlinearFunction(:max, Any[x, x, x])=>maximum([y, y, y]),
        MOI.ScalarNonlinearFunction(:abs, Any[x])=>abs(y),
        MOI.ScalarNonlinearFunction(:sqrt, Any[x])=>sqrt(y),
        MOI.ScalarNonlinearFunction(:exp, Any[x])=>exp(y),
        MOI.ScalarNonlinearFunction(:log, Any[x])=>log(y),
    ]
        @test sprint(show, Convex._expr(model, f)) == sprint(show, g)
        if f isa MOI.ScalarNonlinearFunction
            push!(tested_operators, f.head)
        end
    end
    # Test that we have a test above for every supported operator
    operators = MOI.get(model, MOI.ListOfSupportedNonlinearOperators())
    @test sort(collect(tested_operators)) == sort(operators)
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

function test_not_dcp_constraint()
    inner = Convex.Optimizer(ECOS.Optimizer)
    model = MOI.Utilities.CachingOptimizer(
        MOI.Utilities.UniversalFallback(MOI.Utilities.Model{Float64}()),
        MOI.Bridges.full_bridge_optimizer(inner, Float64),
    )
    x = MOI.add_variable(model)
    f = MOI.ScalarNonlinearFunction(:^, Any[x, 2])
    MOI.add_constraint(model, f, MOI.GreaterThan(1.0))
    F, S = MOI.VectorNonlinearFunction, MOI.Nonnegatives
    @test_throws(
        MOI.AddConstraintNotAllowed{F,S},
        MOI.Utilities.attach_optimizer(model),
    )
    return
end

function test_not_dcp_objective()
    inner = Convex.Optimizer(ECOS.Optimizer)
    model = MOI.Utilities.CachingOptimizer(
        MOI.Utilities.UniversalFallback(MOI.Utilities.Model{Float64}()),
        MOI.Bridges.full_bridge_optimizer(inner, Float64),
    )
    x = MOI.add_variable(model)
    f = MOI.ScalarNonlinearFunction(:^, Any[x, 2])
    MOI.set(model, MOI.ObjectiveSense(), MOI.MAX_SENSE)
    attr = MOI.ObjectiveFunction{typeof(f)}()
    MOI.set(model, attr, f)
    @test_throws(
        MOI.SetAttributeNotAllowed{typeof(attr)},
        MOI.Utilities.attach_optimizer(model),
    )
    return
end

function test_not_dcp_objective_min()
    inner = Convex.Optimizer(ECOS.Optimizer)
    model = MOI.Utilities.CachingOptimizer(
        MOI.Utilities.UniversalFallback(MOI.Utilities.Model{Float64}()),
        MOI.Bridges.full_bridge_optimizer(inner, Float64),
    )
    x = MOI.add_variable(model)
    g = MOI.ScalarNonlinearFunction(:^, Any[x, 2])
    f = MOI.ScalarNonlinearFunction(:-, Any[g])
    MOI.set(model, MOI.ObjectiveSense(), MOI.MIN_SENSE)
    attr = MOI.ObjectiveFunction{typeof(f)}()
    MOI.set(model, attr, f)
    @test_throws(
        MOI.SetAttributeNotAllowed{typeof(attr)},
        MOI.Utilities.attach_optimizer(model),
    )
    return
end

function test_log_det_objective_from_jump()
    n = 2
    model = JuMP.Model(() -> Convex.Optimizer(SCS.Optimizer))
    JuMP.set_silent(model)
    JuMP.@variable(model, X[1:n, 1:n], Symmetric)
    JuMP.@constraint(model, X in JuMP.PSDCone())
    JuMP.@constraint(model, LinearAlgebra.tr(X) <= 4)
    op_det = JuMP.NonlinearOperator(LinearAlgebra.det, :det)
    JuMP.@objective(model, Max, log(op_det(X)))
    attr = MOI.ObjectiveFunction{MOI.ScalarNonlinearFunction}()
    moi_f = MOI.get(model, attr)
    @test moi_f.head == :log
    @test moi_f.args[1] isa MOI.ScalarNonlinearFunction
    @test moi_f.args[1].head == :det
    @test moi_f.args[1].args[1] isa Matrix{MOI.VariableIndex}
    JuMP.optimize!(model)
    @test JuMP.termination_status(model) == MOI.OPTIMAL
    # log det X is maximised on tr(X) ≤ 4, X ⪰ 0 by X = 2I, giving log det = 2 log 2
    @test isapprox(JuMP.objective_value(model), 2 * log(2); atol = 1e-3)
    return
end

end  # TestMOIWrapper

TestMOIWrapper.runtests()
