module TestAtoms

using Convex
using Test

import MathOptInterface as MOI

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

function _test_structural_identical(a::MOI.ModelLike, b::MOI.ModelLike)
    # Test that the variables are the same. We make the strong assumption that
    # the variables are added in the same order to both models.
    a_x = MOI.get(a, MOI.ListOfVariableIndices())
    b_x = MOI.get(b, MOI.ListOfVariableIndices())
    attr = MOI.NumberOfVariables()
    Test.@test MOI.get(a, attr) == MOI.get(b, attr)
    Test.@test length(a_x) == length(b_x)
    # A dictionary that maps things from `b`-space to `a`-space.
    x_map = Dict(bx => a_x[i] for (i, bx) in enumerate(b_x))
    # To check that the constraints, we need to first cache all of the
    # constraints in `a`.
    constraints = Dict{Any,Any}()
    for (F, S) in MOI.get(a, MOI.ListOfConstraintTypesPresent())
        Test.@test MOI.supports_constraint(a, F, S)
        constraints[(F, S)] =
            map(MOI.get(a, MOI.ListOfConstraintIndices{F,S}())) do ci
                return (
                    MOI.get(a, MOI.ConstraintFunction(), ci),
                    MOI.get(a, MOI.ConstraintSet(), ci),
                )
            end
    end
    # Now compare the constraints in `b` with the cache in `constraints`.
    b_constraint_types = MOI.get(b, MOI.ListOfConstraintTypesPresent())
    # There may be constraint types reported in `a` that are not in `b`, but
    # have zero constraints in `a`.
    for (F, S) in keys(constraints)
        attr = MOI.NumberOfConstraints{F,S}()
        Test.@test (F, S) in b_constraint_types || MOI.get(a, attr) == 0
    end
    for (F, S) in b_constraint_types
        # Check that the same number of constraints are present
        attr = MOI.NumberOfConstraints{F,S}()
        if !haskey(constraints, (F, S))
            # Constraint is reported in `b`, but not in `a`. Check that there
            # are no actual constraints in `b`.
            Test.@test MOI.get(b, attr) == 0
            continue
        else
            Test.@test MOI.get(a, attr) == MOI.get(b, attr)
        end
        # Check that supports_constraint is implemented
        Test.@test MOI.supports_constraint(b, F, S)
        # Check that each function in `b` matches a function in `a`
        for ci in MOI.get(b, MOI.ListOfConstraintIndices{F,S}())
            f_b = MOI.get(b, MOI.ConstraintFunction(), ci)
            f_b = MOI.Utilities.map_indices(x_map, f_b)
            s_b = MOI.get(b, MOI.ConstraintSet(), ci)
            # We don't care about the order that constraints are added, only
            # that one matches.
            Test.@test any(constraints[(F, S)]) do (f, s)
                return s_b == s && isapprox(f, f_b) && typeof(f) == typeof(f_b)
            end
        end
    end
    # Test model attributes are set, like ObjectiveSense and ObjectiveFunction.
    a_attrs = MOI.get(a, MOI.ListOfModelAttributesSet())
    b_attrs = MOI.get(b, MOI.ListOfModelAttributesSet())
    Test.@test length(a_attrs) == length(b_attrs)
    for attr in b_attrs
        Test.@test attr in a_attrs
        if attr == MOI.ObjectiveSense()
            # map_indices isn't defined for `OptimizationSense`
            Test.@test MOI.get(a, attr) == MOI.get(b, attr)
        else
            attr_b = MOI.Utilities.map_indices(x_map, MOI.get(b, attr))
            Test.@test isapprox(MOI.get(a, attr), attr_b)
        end
    end
    return
end

function _test_model_equal(f, target_string::String)
    context = Convex.Context{Float64}(MOI.Utilities.Model{Float64})
    f(context)
    target = MOI.Utilities.Model{Float64}()
    MOI.Utilities.loadfromstring!(target, target_string)
    _test_structural_identical(context.model, target)
    return
end

function test_RationalNormAtom()
    _test_model_equal(
        """
        variables: x1, x2, t
        [-1.0+1.0*x1, -2.0+1.0*x2] in Nonnegatives(2)
        [1.0*t, 1.0*x1, 1.0*x2] in NormCone(1.5, 3)
        """,
    ) do context
        x = Variable(2)
        Convex.add_constraint!(context, x >= [1, 2])
        t = Convex.conic_form!(context, rationalnorm(x, 3 // 2))
    end
    _test_model_equal(
        """
        variables: x1, x2, x3, x4, t
        [-1.0+1.0*x1, -3.0+1.0*x2, -2.0+1.0*x3, -4.0+1.0*x4] in Nonnegatives(4)
        [1.0*t, 1.0*x1, 1.0*x2, 1.0*x3, 1.0*x4] in NormCone(2.0, 5)
        """,
    ) do context
        x = Variable(2, 2)
        Convex.add_constraint!(context, x >= [1 2; 3 4])
        t = Convex.conic_form!(context, rationalnorm(x, 2 // 1))
    end
    return
end

function test_RationalNormAtom_complex_matrix()
    x = Variable(2, 2)
    @test_throws(
        ErrorException("[RationalNormAtom] not defined for complex matrices"),
        rationalnorm(im * x, 3 // 2),
    )
    return
end

function test_RationalNormAtom_complex_vector()
    x = Variable(2)
    y = Variable(2)
    x.value = [1.0, -2.0]
    y.value = [-3.0, 4.0]
    atom = rationalnorm(x + im * y, 3 // 2)
    z = abs.([1 - 3im, -2 + 4im])
    @test evaluate(atom) ≈ sum(abs.(z) .^ (3 // 2))^(2 // 3)
    return
end

function test_RationalNormAtom_matrix()
    x = Variable(2, 2)
    atom = rationalnorm(x, 3 // 2)
    x.value = [1.0 2.0; 3.0 4.0]
    @test evaluate(atom) ≈ sum(abs.(x.value) .^ (3 // 2))^(2 // 3)
    return
end

function test_RationalNormAtom_less_than_1()
    x = Variable(3)
    k = 1 // 2
    @test_throws(
        ErrorException(
            "[RationalNormAtom] p-norms not defined for p < 1. Got $k",
        ),
        rationalnorm(x, k),
    )
    return
end

end  # TestAtoms

TestAtoms.runtests()
