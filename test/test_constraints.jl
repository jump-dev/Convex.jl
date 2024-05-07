# Copyright (c) 2014: Madeleine Udell and contributors
#
# Use of this source code is governed by a BSD-style license that can be found
# in the LICENSE file or at https://opensource.org/license/bsd-2-clause

module TestConstraints

using Convex
using Test

import MathOptInterface as MOI
import SCS

function runtests()
    for name in names(@__MODULE__; all = true)
        if startswith("$name", "test_")
            @testset "$name" begin
                getfield(@__MODULE__, name)()
            end
        end
    end
    return
end

### constraints/EqualToConstraint

function test_EqualToConstraint()
    @test_throws(
        ErrorException(
            "Cannot create constraint between expressions of size (2, 3) and (3, 2)",
        ),
        Variable(2, 3) == Variable(3, 2),
    )
    @test vexity(square(Variable()) == 1) == Convex.NotDcp()
    @test vexity(-square(Variable()) == 1) == Convex.NotDcp()
    return
end

function test_EqualToConstraint_violated()
    p = satisfy([constant(5) == 0])
    @test_logs (:info,) (:warn,) solve!(p, SCS.Optimizer)
    return
end

function test_EqualToConstraint_dual_minimize()
    x = Variable()
    c = 3 * x == 1
    p = minimize(2.0 * x + 1.0, [c])
    solve!(p, SCS.Optimizer; silent_solver = true)
    @test isapprox(c.dual, 2 / 3; atol = 1e-4)
    return
end

function test_EqualToConstraint_dual_maximize()
    x = Variable()
    c = 3 * x == 1
    p = maximize(2.0 * x + 1.0, [c])
    solve!(p, SCS.Optimizer; silent_solver = true)
    @test isapprox(c.dual, -2 / 3; atol = 1e-4)
    return
end

function test_EqualToConstraint_dual_complex()
    x = ComplexVariable()
    c = 3 * x == 1 - 2im
    p = minimize(real(x) + 2imag(x), [c])
    solve!(p, SCS.Optimizer; silent_solver = true)
    @test isapprox(c.dual, 1 / 3 - 2 / 3im; atol = 1e-4)
    return
end

### constraints/GreaterThanConstraint

function test_GreaterThanConstraint()
    @test_throws(
        ErrorException(
            "Cannot create constraint between expressions of size (2, 3) and (3, 2)",
        ),
        Variable(2, 3) >= Variable(3, 2),
    )
    @test_throws(
        ErrorException(
            "Cannot create constraint between expressions of sign $(Convex.NoSign()) and $(Convex.ComplexSign())",
        ),
        Variable() >= 2 + 3im,
    )
    @test_throws(
        ErrorException(
            "Cannot create constraint between expressions of sign $(Convex.ComplexSign()) and $(Convex.NoSign())",
        ),
        2 + 3im >= Variable(),
    )
    @test vexity(square(Variable()) >= 1) == Convex.NotDcp()
    return
end

function test_GreaterThanConstraint_violated()
    p = satisfy([constant(5) >= 6])
    @test_logs (:info,) (:warn,) solve!(p, SCS.Optimizer)
    return
end

function test_GreaterThanConstraint_dual_minimize()
    x = Variable()
    c = 3 * x >= 1
    p = minimize(2.0 * x + 1.0, [c])
    solve!(p, SCS.Optimizer; silent_solver = true)
    @test isapprox(c.dual, 2 / 3; atol = 1e-4)
    return
end

function test_GreaterThanConstraint_dual_maximize()
    x = Variable()
    c = 1 >= 3 * x
    p = maximize(2.0 * x + 1.0, [c])
    solve!(p, SCS.Optimizer; silent_solver = true)
    @test isapprox(c.dual, 2 / 3; atol = 1e-4)
    return
end

### constraints/LessThanConstraint

function test_LessThanConstraint()
    @test_throws(
        ErrorException(
            "Cannot create constraint between expressions of size (2, 3) and (3, 2)",
        ),
        Variable(2, 3) <= Variable(3, 2),
    )
    @test_throws(
        ErrorException(
            "Cannot create constraint between expressions of sign $(Convex.NoSign()) and $(Convex.ComplexSign())",
        ),
        Variable() <= 2 + 3im,
    )
    @test_throws(
        ErrorException(
            "Cannot create constraint between expressions of sign $(Convex.ComplexSign()) and $(Convex.NoSign())",
        ),
        2 + 3im <= Variable(),
    )
    @test vexity(-square(Variable()) <= 1) == Convex.NotDcp()
    return
end

function test_LessThanConstraint_violated()
    p = satisfy([constant(5) <= 4])
    @test_logs (:info,) (:warn,) solve!(p, SCS.Optimizer)
    return
end

function test_LessThanConstraint_dual_minimize()
    x = Variable()
    c = 1 <= 3 * x
    p = minimize(2.0 * x + 1.0, [c])
    solve!(p, SCS.Optimizer; silent_solver = true)
    @test isapprox(c.dual, -2 / 3; atol = 1e-4)
    return
end

function test_LessThanConstraint_dual_maximize()
    x = Variable()
    c = 3 * x <= 1
    p = maximize(2.0 * x + 1.0, [c])
    solve!(p, SCS.Optimizer; silent_solver = true)
    @test isapprox(c.dual, -2 / 3; atol = 1e-4)
    return
end

### constraints/GenericConstraint{MOI.PositiveSemidefiniteConeSquare}

function test_GenericConstraint_PositiveSemidefiniteConeSquare()
    @test_throws(
        ErrorException("Positive semidefinite expressions must be square"),
        Convex.GenericConstraint{MOI.PositiveSemidefiniteConeSquare}(
            Variable(2, 3),
        ),
    )
    X = Variable(2, 2)
    c = Convex.GenericConstraint{MOI.PositiveSemidefiniteConeSquare}(X)
    p = minimize(tr(X), [c, X >= [1 2; 3 4]])
    solve!(p, SCS.Optimizer; silent_solver = true)
    @test isapprox(X.value, [2.25 3; 3 4]; atol = 1e-3)
    y = (c.dual + c.dual') / 2
    @test isapprox(y[1], 1; atol = 1e-3)
    @test (0 ⪯ X) isa
          Convex.GenericConstraint{MOI.PositiveSemidefiniteConeSquare}
    @test (-X ⪯ 0) isa
          Convex.GenericConstraint{MOI.PositiveSemidefiniteConeSquare}
    @test (-X ⪯ constant(0)) isa
          Convex.GenericConstraint{MOI.PositiveSemidefiniteConeSquare}
    @test (constant(0) ⪯ X) isa
          Convex.GenericConstraint{MOI.PositiveSemidefiniteConeSquare}
    @test_throws(ErrorException("Set PSD not understood"), X in :PSD)
    @test vexity(X ⪯ square(Variable())) == Convex.NotDcp()
    return
end

function test_GenericConstraint_PositiveSemidefiniteConeSquare_violated()
    X = constant([1 2; 3 4])
    p = satisfy([X ⪰ 0])
    @test_logs (:info,) (:warn,) solve!(p, SCS.Optimizer)
    X = constant([1 2; 2 3])
    p = satisfy([X ⪰ 0])
    @test_logs (:info,) (:warn,) solve!(p, SCS.Optimizer)
    return
end

### constraints/GenericConstraint_SecondOrderCone

function test_GenericConstraint_SecondOrderCone()
    x = Variable(3)
    t = Variable()
    c = Convex.GenericConstraint{MOI.SecondOrderCone}(vcat(t, x))
    p = minimize(t, [c, x >= [2, 3, 4]])
    solve!(p, SCS.Optimizer; silent_solver = true)
    @test isapprox(x.value, [2, 3, 4]; atol = 1e-3)
    t_ = sqrt(29)
    @test isapprox(t.value, t_; atol = 1e-3)
    @test isapprox(c.dual, [1, -2 / t_, -3 / t_, -4 / t_]; atol = 1e-3)
    c = Convex.GenericConstraint{MOI.SecondOrderCone}(vcat(square(t), x))
    @test vexity(c) === Convex.NotDcp()
    @test Convex.sprint(Convex.head, c) == "soc"
    return
end

function test_GenericConstraint_SecondOrderCone_set_with_size()
    x = Variable(2, 2)
    S = MOI.SecondOrderCone
    sz = size(x)
    @test_throws(
        ErrorException(
            "Cannot constrain a matrix of size `$sz` to be long to the cone " *
            "`$S`, there should be only one column.",
        ),
        Convex.GenericConstraint{S}(x),
    )
    return
end

### constraints/GenericConstraint_RotatedSecondOrderCone

function test_GenericConstraint_RotatedSecondOrderCone()
    x = Variable(3)
    t = Variable()
    c = Convex.GenericConstraint{MOI.RotatedSecondOrderCone}(vcat(t, 1, x))
    p = minimize(t, [c, x >= [2, 3, 4]])
    solve!(p, SCS.Optimizer; silent_solver = true)
    @test isapprox(x.value, [2, 3, 4]; atol = 1e-3)
    @test isapprox(t.value, 29 / 2; atol = 1e-3)
    @test isapprox(c.dual, [1, 29 / 2, -2, -3, -4]; atol = 1e-3)
    c = Convex.GenericConstraint{MOI.RotatedSecondOrderCone}(vcat(square(t), x))
    @test vexity(c) === Convex.NotDcp()
    @test Convex.sprint(Convex.head, c) == "rsoc"
    return
end

### constraints/GenericConstraint_ExponentialCone

function test_GenericConstraint_ExponentialConeConstraint()
    # y * exp(x / y) <= z  <=>  (x, y, z) in ExpCone
    z = Variable()
    # 1 * exp(1 / 1) <= z
    c = Convex.GenericConstraint{MOI.ExponentialCone}(vcat(1, constant(1), z))
    p = minimize(z, [c])
    solve!(p, SCS.Optimizer; silent_solver = true)
    @test isapprox(z.value, exp(1); atol = 1e-4)
    @test isapprox(c.dual, [-exp(1), 0, 1]; atol = 1e-4)
    z = Variable()
    # 2 * exp(3 / 2) <= z
    c = Convex.GenericConstraint{MOI.ExponentialCone}(vcat(constant(3), 2, z))
    p = minimize(z, [c])
    solve!(p, SCS.Optimizer; silent_solver = true)
    @test isapprox(z.value, 2 * exp(3 / 2); atol = 1e-4)
    @test isapprox(c.dual, [-exp(3 / 2), exp(3 / 2) / 2, 1]; atol = 1e-4)
    c = Convex.GenericConstraint{MOI.ExponentialCone}(vcat(square(z), 1, z))
    @test vexity(c) === Convex.NotDcp()
    @test sprint(Convex.head, c) == "exp"
    return
end

### constraints/GenericConstraint_RelativeEntropyCone

function test_GenericConstraint_RelativeEntropyConeConstraint()
    # u >= sum_i( w_i log (w_i/v_i) )
    u = Variable()
    v = [1, 2]
    w = Variable(2)
    c = Convex.GenericConstraint{MOI.RelativeEntropyCone}(vcat(u, v, w))
    p = minimize(u, [c, w >= 2])
    solve!(p, SCS.Optimizer; silent_solver = true)
    @test isapprox(u.value, 2 * log(2 / 1) + 2 * log(2 / 2); atol = 1e-4)
    @test isapprox(c.dual, [1, 2, 1, -log(2) - 1, -1]; atol = 1e-3)
    c = Convex.GenericConstraint{MOI.RelativeEntropyCone}(vcat(square(u), 1, u))
    @test vexity(c) === Convex.NotDcp()
    @test sprint(Convex.head, c) == "RelativeEntropyCone"
    return
end

end  # TestConstraints

TestConstraints.runtests()
