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
    solve!(p, SCS.Optimizer; silent = true)
    @test isapprox(c.dual, 2 / 3; atol = 1e-4)
    return
end

function test_EqualToConstraint_dual_maximize()
    x = Variable()
    c = 3 * x == 1
    p = maximize(2.0 * x + 1.0, [c])
    solve!(p, SCS.Optimizer; silent = true)
    @test isapprox(c.dual, -2 / 3; atol = 1e-4)
    return
end

function test_EqualToConstraint_dual_complex()
    x = ComplexVariable()
    c = 3 * x == 1 - 2im
    p = minimize(real(x) + 2imag(x), [c])
    solve!(p, SCS.Optimizer; silent = true)
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
    solve!(p, SCS.Optimizer; silent = true)
    @test isapprox(c.dual, 2 / 3; atol = 1e-4)
    return
end

function test_GreaterThanConstraint_dual_maximize()
    x = Variable()
    c = 1 >= 3 * x
    p = maximize(2.0 * x + 1.0, [c])
    solve!(p, SCS.Optimizer; silent = true)
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
    solve!(p, SCS.Optimizer; silent = true)
    @test isapprox(c.dual, -2 / 3; atol = 1e-4)
    return
end

function test_LessThanConstraint_dual_maximize()
    x = Variable()
    c = 3 * x <= 1
    p = maximize(2.0 * x + 1.0, [c])
    solve!(p, SCS.Optimizer; silent = true)
    @test isapprox(c.dual, -2 / 3; atol = 1e-4)
    return
end

### constraints/Constraint{MOI.PositiveSemidefiniteConeSquare}

function test_Constraint_PositiveSemidefiniteConeSquare()
    @test_throws(
        ErrorException("Positive semidefinite expressions must be square"),
        Convex.Constraint{MOI.PositiveSemidefiniteConeSquare}(Variable(2, 3)),
    )
    X = Variable(2, 2)
    c = Convex.Constraint{MOI.PositiveSemidefiniteConeSquare}(X)
    p = minimize(tr(X), [c, X >= [1 2; 3 4]])
    solve!(p, SCS.Optimizer; silent = true)
    @test isapprox(X.value, [2.25 3; 3 4]; atol = 1e-3)
    y = (c.dual + c.dual') / 2
    @test isapprox(y[1], 1; atol = 1e-3)
    @test (0 ⪯ X) isa Convex.Constraint{MOI.PositiveSemidefiniteConeSquare}
    @test (-X ⪯ 0) isa Convex.Constraint{MOI.PositiveSemidefiniteConeSquare}
    @test (-X ⪯ constant(0)) isa
          Convex.Constraint{MOI.PositiveSemidefiniteConeSquare}
    @test (constant(0) ⪯ X) isa
          Convex.Constraint{MOI.PositiveSemidefiniteConeSquare}
    @test_throws(ErrorException("Set PSD not understood"), X in :PSD)
    @test vexity(X ⪯ square(Variable())) == Convex.NotDcp()
    return
end

function test_Constraint_PositiveSemidefiniteConeSquare_violated()
    X = constant([1 2; 3 4])
    p = satisfy([X ⪰ 0])
    @test_logs (:info,) (:warn,) solve!(p, SCS.Optimizer)
    X = constant([1 2; 2 3])
    p = satisfy([X ⪰ 0])
    @test_logs (:info,) (:warn,) solve!(p, SCS.Optimizer)
    return
end

### constraints/Constraint_SecondOrderCone

function test_Constraint_SecondOrderCone()
    x = Variable(3)
    t = Variable()
    c = Convex.Constraint{MOI.SecondOrderCone}(vcat(t, x))
    p = minimize(t, [c, x >= [2, 3, 4]])
    solve!(p, SCS.Optimizer; silent = true)
    @test isapprox(x.value, [2, 3, 4]; atol = 1e-3)
    t_ = sqrt(29)
    @test isapprox(t.value, t_; atol = 1e-3)
    @test isapprox(c.dual, [1, -2 / t_, -3 / t_, -4 / t_]; atol = 1e-3)
    c = Convex.Constraint{MOI.SecondOrderCone}(vcat(square(t), x))
    @test vexity(c) === Convex.NotDcp()
    @test Convex.sprint(Convex.head, c) == "MOI.SecondOrderCone"
    return
end

function test_Constraint_SecondOrderCone_set_with_size()
    x = Variable(2, 2)
    S = MOI.SecondOrderCone
    sz = size(x)
    @test_throws(
        ErrorException(
            "Cannot constrain a matrix of size `$sz` to be long to the cone " *
            "`$S`, there should be only one column.",
        ),
        Convex.Constraint{S}(x),
    )
    return
end

### constraints/Constraint_RotatedSecondOrderCone

function test_Constraint_RotatedSecondOrderCone()
    x = Variable(3)
    t = Variable()
    c = Convex.Constraint{MOI.RotatedSecondOrderCone}(vcat(t, 1, x))
    p = minimize(t, [c, x >= [2, 3, 4]])
    solve!(p, SCS.Optimizer; silent = true)
    @test isapprox(x.value, [2, 3, 4]; atol = 1e-3)
    @test isapprox(t.value, 29 / 2; atol = 1e-3)
    @test isapprox(c.dual, [1, 29 / 2, -2, -3, -4]; atol = 1e-3)
    c = Convex.Constraint{MOI.RotatedSecondOrderCone}(vcat(square(t), x))
    @test vexity(c) === Convex.NotDcp()
    @test Convex.sprint(Convex.head, c) == "MOI.RotatedSecondOrderCone"
    return
end

### constraints/Constraint_ExponentialCone

function test_Constraint_ExponentialConeConstraint()
    # y * exp(x / y) <= z  <=>  (x, y, z) in ExpCone
    z = Variable()
    # 1 * exp(1 / 1) <= z
    c = Convex.Constraint{MOI.ExponentialCone}(vcat(1, constant(1), z))
    p = minimize(z, [c])
    solve!(p, SCS.Optimizer; silent = true)
    @test isapprox(z.value, exp(1); atol = 1e-4)
    @test isapprox(c.dual, [-exp(1), 0, 1]; atol = 1e-4)
    z = Variable()
    # 2 * exp(3 / 2) <= z
    c = Convex.Constraint{MOI.ExponentialCone}(vcat(constant(3), 2, z))
    p = minimize(z, [c])
    solve!(p, SCS.Optimizer; silent = true)
    @test isapprox(z.value, 2 * exp(3 / 2); atol = 1e-4)
    @test isapprox(c.dual, [-exp(3 / 2), exp(3 / 2) / 2, 1]; atol = 1e-4)
    c = Convex.Constraint{MOI.ExponentialCone}(vcat(square(z), 1, z))
    @test vexity(c) === Convex.NotDcp()
    @test sprint(Convex.head, c) == "MOI.ExponentialCone"
    return
end

### constraints/Constraint_RelativeEntropyCone

function test_Constraint_RelativeEntropyConeConstraint()
    # u >= sum_i( w_i log (w_i/v_i) )
    u = Variable()
    v = [1, 2]
    w = Variable(2)
    c = Convex.Constraint{MOI.RelativeEntropyCone}(vcat(u, v, w))
    p = minimize(u, [c, w >= 2])
    solve!(p, SCS.Optimizer; silent = true)
    @test isapprox(u.value, 2 * log(2 / 1) + 2 * log(2 / 2); atol = 1e-4)
    @test isapprox(c.dual, [1, 2, 1, -log(2) - 1, -1]; atol = 1e-3)
    c = Convex.Constraint{MOI.RelativeEntropyCone}(vcat(square(u), 1, u))
    @test vexity(c) === Convex.NotDcp()
    @test sprint(Convex.head, c) == "MOI.RelativeEntropyCone"
    return
end

function test_Constraint_NoVexity()
    x = Variable(2)
    set = MOI.Complements(2)
    @test_throws(
        ErrorException(
            "`Convex.vexity(vex, ::$(typeof(set)))`: is not yet implemented. Please open an issue at https://github.com/jump-dev/Convex.jl",
        ),
        vexity(Convex.Constraint(x, set)),
    )
    return
end

### constraints/Constraint_GeometricMeanEpiConeSquare

function test_GeometricMeanEpiConeSquare()
    T = Variable(2, 2)
    A = Variable(2, 2)
    B = Variable(2, 2)
    C = [1 0; 0 1]
    @test_throws DomainError GeometricMeanEpiConeSquare(-3 // 2, 3)
    @test_throws DomainError GeometricMeanEpiConeSquare(1 // 2, 3)
    @test_throws DomainError GeometricMeanEpiConeSquare(5 // 2, 3)
    set = Convex.GeometricMeanEpiConeSquare(3 // 2, 2)
    @test sprint(Convex.head, set) == "GeometricMeanEpiConeSquare"
    @test MOI.dimension(set) == 12
    for (f, dcp) in (
        (T, A, B) => Convex.ConvexVexity(),
        (T, A, C) => Convex.ConvexVexity(),
        (T, C, B) => Convex.ConvexVexity(),
        (C, A, B) => Convex.ConvexVexity(),
        (T * T, A, B) => Convex.NotDcp(),
        (T .^ 2, A, B) => Convex.NotDcp(),
        (T, A * A, B) => Convex.NotDcp(),
        (T, sqrt(A), B) => Convex.NotDcp(),
        (T, A, B * B) => Convex.NotDcp(),
        (T, A, sqrt(B)) => Convex.NotDcp(),
        (sqrt(T), A, B) => Convex.ConvexVexity(),
        (T + C, A, B) => Convex.ConvexVexity(),
        (T * C, A, B) => Convex.ConvexVexity(),
    )
        c = Convex.Constraint(f, set)
        @test vexity(c) == dcp
        @test sign(c.child) == Convex.NoSign()
    end
    Z = Variable(3, 3)
    @test_throws DimensionMismatch Convex.Constraint((Z, A, B), set)
    @test_throws DimensionMismatch Convex.Constraint((T, Z, B), set)
    @test_throws DimensionMismatch Convex.Constraint((T, A, Z), set)
    return
end

### RelativeEntropyEpiConeSquare

function test_RelativeEntropyEpiConeSquare()
    set = Convex.RelativeEntropyEpiConeSquare(3)
    @test sprint(Convex.head, set) == "RelativeEntropyEpiConeSquare"
    @test MOI.dimension(set) == 3 * 3^2
    e = [1 2; 3 4; 5 6]
    set = Convex.RelativeEntropyEpiConeSquare(3, 3, 2, e)
    @test MOI.dimension(set) == 2^2 + 2 * 3^2
    @test_throws(
        DimensionMismatch("e matrix must have 4 rows"),
        Convex.RelativeEntropyEpiConeSquare(4, 3, 2, e),
    )
    τ = Variable(2, 2)
    X = Variable(3, 3)
    Y = Variable(3, 3)
    C = [1 0 0; 0 1 0; 0 0 1]
    C2 = [1 0; 0 1]
    for (f, dcp) in (
        (τ, X, Y) => Convex.ConvexVexity(),
        (τ, X, C) => Convex.ConvexVexity(),
        (τ, C, Y) => Convex.ConvexVexity(),
        (C2, X, Y) => Convex.ConvexVexity(),
        (τ * τ, X, Y) => Convex.NotDcp(),
        (τ .^ 2, X, Y) => Convex.NotDcp(),
        (τ, X * X, Y) => Convex.NotDcp(),
        (τ, sqrt(X), Y) => Convex.NotDcp(),
        (τ, X, Y * Y) => Convex.NotDcp(),
        (τ, X, sqrt(Y)) => Convex.NotDcp(),
        (sqrt(τ), X, Y) => Convex.ConvexVexity(),
        (τ + C2, X, Y) => Convex.ConvexVexity(),
        (τ * C2, X, Y) => Convex.ConvexVexity(),
    )
        c = Convex.Constraint(f, set)
        @test vexity(c) == dcp
        @test sign(c.child) == Convex.NoSign()
    end
    @test_throws(
        DimensionMismatch,
        Convex.Constraint((Variable(3, 1), X, Y), set)
    )
    @test_throws(DimensionMismatch, Convex.Constraint((Variable(), X, Y), set),)
    @test_throws(
        DimensionMismatch,
        Convex.Constraint((τ, Variable(2, 3), Y), set)
    )
    @test_throws(
        DimensionMismatch,
        Convex.Constraint((τ, X, Variable(2, 3)), set)
    )
    return
end

function test_distance_to_set_matrix()
    x = Variable(2, 2)
    y = Variable()
    fix!(x, [1 0; 0 1])
    # Constraint has a fixed `Matrix` value.
    model = minimize(y, [sum(x; dims = 1) <= 1, y >= 1])
    solve!(model, SCS.Optimizer; silent = true)
    @test ≈(model.optval, 1.0; atol = 1e-6)
    return
end

function test_distance_to_set_undefined()
    t = Variable()
    fix!(t, 2)
    x = Variable(2, 2)
    fix!(x, [1 0; 0 1])
    y = Variable()
    # This constraint is fixed, and `MOI.distance_to_set` is not defined for it,
    # but it should still work without erroring.
    c = Convex.Constraint(vcat(t, vec(x)), MOI.NormSpectralCone(2, 2))
    model = minimize(y, [c, y >= 1])
    solve!(model, SCS.Optimizer; silent = true)
    @test ≈(model.optval, 1.0; atol = 1e-6)
    return
end

end  # TestConstraints

TestConstraints.runtests()
