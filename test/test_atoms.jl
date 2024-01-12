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

function _to_moi(x::MOI.AbstractVectorFunction)
    if MOI.output_dimension(x) > 1
        return x
    end
    return only(MOI.Utilities.scalarize(x))
end

_to_moi(x::Convex.SparseTape) = _to_moi(Convex.to_vaf(x))

_to_moi(v::MOI.AbstractScalarFunction) = v

"""
    _test_atom(f, target_string::String; value_type = Float64)
"""
function _test_atom(build_fn, target_string::String; value_type = Float64)
    context = Convex.Context{value_type}(MOI.Utilities.Model{value_type})
    atom = build_fn(context)
    # All atoms must be an AbstractExpr
    @test atom isa Convex.AbstractExpr
    # All atoms must be mutable
    @test ismutable(atom)
    @test sprint(Convex.head, atom) isa String
    @test Base.sign(atom) isa Convex.Sign
    N = length(atom.children)
    @test Convex.monotonicity(atom) isa NTuple{N,<:Convex.Monotonicity}
    @test Convex.curvature(atom) isa Convex.Vexity
    t = Convex.conic_form!(context, atom)
    MOI.set(context.model, MOI.ObjectiveSense(), MOI.MIN_SENSE)
    obj = _to_moi(t)
    MOI.set(context.model, MOI.ObjectiveFunction{typeof(obj)}(), obj)
    target = MOI.Utilities.Model{value_type}()
    MOI.Utilities.loadfromstring!(target, target_string)
    # Use the same names for each model
    for (x, y) in zip(
        MOI.get(context.model, MOI.ListOfVariableIndices()),
        MOI.get(target, MOI.ListOfVariableIndices()),
    )
        name = MOI.get(target, MOI.VariableName(), y)
        MOI.set(context.model, MOI.VariableName(), x, name)
    end
    context_string = sprint(print, context.model)
    target_string = sprint(print, target)
    if context_string != target_string
        @info "Target model\n$target_string"
        @info "context.model\n$context_string"
    end
    @test context_string == target_string
    return
end

### affine/AdditionAtom

function test_AdditionAtom()
    target = """
    variables: x
    minobjective: 1.0 + 1.0 * x
    """
    _test_atom(target) do context
        x = Variable()
        return x + 1
    end
    target = """
    variables: x
    minobjective: 2.0 + 1.0 * x
    """
    _test_atom(target) do context
        x = Variable()
        return 2.0 + x
    end
    target = """
    variables: x
    minobjective: 2.0 + 3.0 * x
    """
    _test_atom(target) do context
        x = Variable()
        return 2.0 + 3.0 * x
    end
    target = """
    variables: x
    minobjective: 0.0 + 5.0 * x
    """
    _test_atom(target) do context
        x = Variable()
        return 2.0 * x + 3.0 * x
    end
    target = """
    variables: x1, x2
    minobjective: [2.0 + x1, 2.0 + x2]
    """
    _test_atom(target) do context
        x = Variable(2)
        return x + 2
    end
    target = """
    variables: x1, x2
    minobjective: [2.0 + x1, 2.0 + x2]
    """
    _test_atom(target) do context
        x = Variable(2)
        return 2 + x
    end
    target = """
    variables: x
    minobjective: 2.0 + 3.0 * x
    """
    _test_atom(target) do context
        x = Variable()
        a = x + x
        b = 2 + x
        return a + b
    end
    target = """
    variables: x1, x2, y
    minobjective: [1.0 * x1 + 1.0 * y, 1.0 * x2 + 1.0 * y]
    """
    _test_atom(target) do context
        x = Variable(2)
        y = Variable()
        return x + y
    end
    # TODO(odow): fix me
    # target = """
    # variables: x1, x2, y
    # minobjective: [1.0 * x1 + 1.0 * y, 1.0 * x2 + 1.0 * y]
    # """
    # _test_atom(target) do context
    #     x = Variable(2)
    #     y = Variable()
    #     return y + x
    # end
    return
end

function test_AdditionAtom_negate()
    target = """
    variables: x
    minobjective: -1.0 + 1.0 * x
    """
    _test_atom(target) do context
        x = Variable()
        return x - 1
    end
    target = """
    variables: x
    minobjective: 2.0 + -1.0 * x
    """
    _test_atom(target) do context
        x = Variable()
        return 2.0 - x
    end
    target = """
    variables: x
    minobjective: 2.0 + -3.0 * x
    """
    _test_atom(target) do context
        x = Variable()
        return 2.0 - 3.0 * x
    end
    target = """
    variables: x
    minobjective: 0.0 + -1.0 * x
    """
    _test_atom(target) do context
        x = Variable()
        return 2.0 * x - 3.0 * x
    end
    return
end

function test_AdditionAtom_errors()
    x = Variable(2, 2)
    y = Variable(2, 3)
    @test_throws(
        ErrorException(
            "[AdditionAtom] cannot add expressions of sizes $(x.size) and $(y.size)",
        ),
        x + y,
    )
    return
end

### second_order_cone/RationalNormAtom

function test_RationalNormAtom()
    target = """
    variables: x1, x2, t
    minobjective: 1.0 * t + 0.0
    [-1.0+1.0*x1, -2.0+1.0*x2] in Nonnegatives(2)
    [1.0*t, 1.0*x1, 1.0*x2] in NormCone(1.5, 3)
    """
    _test_atom(target) do context
        x = Variable(2)
        Convex.add_constraint!(context, x >= [1, 2])
        return rationalnorm(x, 3 // 2)
    end
    _test_atom(
        """
        variables: x1, x2, x3, x4, t
        minobjective: 1.0 * t + 0.0
        [-1.0+1.0*x1, -3.0+1.0*x2, -2.0+1.0*x3, -4.0+1.0*x4] in Nonnegatives(4)
        [1.0*t, 1.0*x1, 1.0*x2, 1.0*x3, 1.0*x4] in NormCone(2.0, 5)
        """,
    ) do context
        x = Variable(2, 2)
        Convex.add_constraint!(context, x >= [1 2; 3 4])
        return rationalnorm(x, 2 // 1)
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
