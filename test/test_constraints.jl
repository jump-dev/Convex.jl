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
            "Cannot create equality constraint between expressions of size (2, 3) and (3, 2)",
        ),
        Variable(2, 3) == Variable(3, 2),
    )
    @test vexity(square(Variable()) == 1) == Convex.NotDcp()
    @test vexity(-square(Variable()) == 1) == Convex.NotDcp()
    return
end

function test_EqualToConstraint_violated()
    p = satisfy([constant(5) == 0])
    @test_logs (:warn,) (:warn,) solve!(p, SCS.Optimizer)
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

### constraints/ExponentialConeConstraint

function test_ExponentialConeConstraint()
    # y * exp(x / y) <= z  <=>  (x, y, z) in ExpCone
    z = Variable()
    # 1 * exp(1 / 1) <= z
    c = Convex.ExponentialConeConstraint(1, constant(1), z)
    p = minimize(z, [c])
    solve!(p, SCS.Optimizer; silent_solver = true)
    @test isapprox(z.value, exp(1); atol = 1e-4)
    @test isapprox(c.dual, [-exp(1), 0, 1]; atol = 1e-4)
    z = Variable()
    # 2 * exp(3 / 2) <= z
    c = Convex.ExponentialConeConstraint(constant(3), 2, z)
    p = minimize(z, [c])
    solve!(p, SCS.Optimizer; silent_solver = true)
    @test isapprox(z.value, 2 * exp(3 / 2); atol = 1e-4)
    @test isapprox(c.dual, [-exp(3 / 2), exp(3 / 2) / 2, 1]; atol = 1e-4)
    return
end

### constraints/GreaterThanConstraint

function test_GreaterThanConstraint()
    @test_throws(
        ErrorException(
            "Cannot create inequality constraint between expressions of size (2, 3) and (3, 2)",
        ),
        Variable(2, 3) >= Variable(3, 2),
    )
    @test_throws(
        ErrorException(
            "Cannot create inequality constraint between expressions of sign $(Convex.NoSign()) and $(Convex.ComplexSign())",
        ),
        Variable() >= 2 + 3im,
    )
    @test_throws(
        ErrorException(
            "Cannot create inequality constraint between expressions of sign $(Convex.ComplexSign()) and $(Convex.NoSign())",
        ),
        2 + 3im >= Variable(),
    )
    @test vexity(square(Variable()) >= 1) == Convex.NotDcp()
    return
end

function test_GreaterThanConstraint_violated()
    p = satisfy([constant(5) >= 6])
    @test_logs (:warn,) (:warn,) solve!(p, SCS.Optimizer)
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
            "Cannot create inequality constraint between expressions of size (2, 3) and (3, 2)",
        ),
        Variable(2, 3) <= Variable(3, 2),
    )
    @test_throws(
        ErrorException(
            "Cannot create inequality constraint between expressions of sign $(Convex.NoSign()) and $(Convex.ComplexSign())",
        ),
        Variable() <= 2 + 3im,
    )
    @test_throws(
        ErrorException(
            "Cannot create inequality constraint between expressions of sign $(Convex.ComplexSign()) and $(Convex.NoSign())",
        ),
        2 + 3im <= Variable(),
    )
    @test vexity(-square(Variable()) <= 1) == Convex.NotDcp()
    return
end

function test_LessThanConstraint_violated()
    p = satisfy([constant(5) <= 4])
    @test_logs (:warn,) (:warn,) solve!(p, SCS.Optimizer)
    return
end

function test_LessThanConstraint_dual_minimize()
    x = Variable()
    c = 1 <= 3 * x
    p = minimize(2.0 * x + 1.0, [c])
    solve!(p, SCS.Optimizer; silent_solver = true)
    # This dual is currently positive despite a <= constraint because Convex
    # implements rhs - lhs in R_+ instead of lhs - rhs in R_-.
    @test isapprox(c.dual, 2 / 3; atol = 1e-4)
    return
end

function test_LessThanConstraint_dual_maximize()
    x = Variable()
    c = 3 * x <= 1
    p = maximize(2.0 * x + 1.0, [c])
    solve!(p, SCS.Optimizer; silent_solver = true)
    # This dual is currently positive despite a <= constraint because Convex
    # implements rhs - lhs in R_+ instead of lhs - rhs in R_-.
    @test isapprox(c.dual, 2 / 3; atol = 1e-4)
    return
end

end  # TestConstraints

TestConstraints.runtests()
