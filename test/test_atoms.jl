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
    @test Convex.monotonicity(atom) isa NTuple{N,Convex.Monotonicity}
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
    target = """
    variables: y, x1, x2
    minobjective: [1.0 * y + 1.0 * x1, 1.0 * y + 1.0 * x2]
    """
    _test_atom(target) do context
        x = Variable(2)
        y = Variable()
        return y + x
    end
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

### affine/ConjugateAtom

function test_ConjugateAtom()
    x = Variable()
    y = constant(2.0)
    @test isequal(conj(x), x)
    @test isequal(conj(y), y)
    z = conj(im * x)
    @test z isa Convex.ConjugateAtom
    @test sprint(Convex.head, z) == "conj"
    @test Base.sign(z) == Convex.ComplexSign()
    @test Convex.monotonicity(z) == (Convex.Nondecreasing(),)
    @test Convex.curvature(z) == Convex.ConstVexity()
    @test Convex.evaluate(conj(Convex.ComplexConstant(y, y))) == 2.0 - 2.0im
    x.value = [-2.0]
    @test evaluate(conj(x)) ≈ -2.0
    @test evaluate(conj(im * x)) ≈ 0 + 2.0im
    return
end

### affine/DiagAtom

function test_DiagAtom()
    target = """
    variables: x11, x21, x21, x22
    minobjective: [1.0 * x11, 1.0 * x22]
    """
    _test_atom(target) do context
        return diag(Variable(2, 2))
    end
    target = """
    variables: x11, x21, x12
    minobjective: 0.0 + 1.0 * x12
    """
    _test_atom(target) do context
        return diag(Variable(2, 2), 1)
    end
    target = """
    variables: x11, x21
    minobjective: 1.0 * x21
    """
    _test_atom(target) do context
        return diag(Variable(2, 2), -1)
    end
    @test_throws(
        ErrorException(
            "[DiagAtom] bounds error in calling diag. Got 3 but it must be in -2..2",
        ),
        diag(Variable(2, 2), 3),
    )
    @test_throws(
        ErrorException(
            "[DiagAtom] bounds error in calling diag. Got -5 but it must be in -3..3",
        ),
        diag(Variable(3, 4), -5),
    )
    return
end

### affine/DiagMatrixAtom

function test_DiagMatrixAtom()
    target = """
    variables: x, y
    minobjective: [1.0 * x, 0.0, 0.0, 1.0 * y]

    """
    _test_atom(target) do context
        return diagm(Variable(2))
    end
    _test_atom(target) do context
        return diagm(0 => Variable(2))
    end
    _test_atom(target) do context
        return Diagonal(Variable(2))
    end
    _test_atom(target) do context
        return Diagonal(Variable(1, 2))
    end
    _test_atom(target) do context
        return Diagonal(Variable(2, 1))
    end
    @test_throws(
        ArgumentError(
            "[DiagMatrixAtom] only vectors are allowed for `diagm(x)` and `Diagonal(x). Did you mean to use `diag(x, 0)`?",
        ),
        diagm(Variable(2, 2)),
    )
    @test_throws(
        ArgumentError(
            "[DiagMatrixAtom] only the main diagonal is supported. Got `d=1`",
        ),
        diagm(1 => Variable(2)),
    )
    return
end

### affine/HcatAtom

function test_HcatAtom()
    target = """
    variables: x
    minobjective: [1.0 * x, 1.0 * x]
    """
    _test_atom(target) do context
        x = Variable()
        return hcat(x, x)
    end
    _test_atom(target) do context
        x = Variable()
        return vcat(x, x)
    end
    target = """
    variables: x1, x2
    minobjective: [1.0 * x1, 1.0 * x2, 2.0]
    """
    _test_atom(target) do context
        x = Variable(2)
        y = constant(2)
        return vcat(x, y)
    end
    _test_atom(target) do context
        x = Variable(2)
        return vcat(x, 2)
    end
    _test_atom(target) do context
        x = Variable(1, 2)
        return hcat(x, 2)
    end
    @test_throws(
        DimensionMismatch(
            "[HcatAtom] cannot stack expressions of incompatible size. Got 1 expected 2.",
        ),
        hcat(Variable(2), constant(2)),
    )
    @test_throws(
        DimensionMismatch(
            "[HcatAtom] cannot stack expressions of incompatible size. Got 2 expected 1.",
        ),
        vcat(Variable(2, 1), Variable(1, 2)),
    )
    return
end

### affine/ImaginaryAtom

function test_ImaginaryAtom()
    target = """
    variables: x
    minobjective: 1.0 * x
    """
    _test_atom(target) do context
        return imag(im * Variable())
    end
    target = """
    variables: x
    minobjective: 1.0 * x + 3.0
    """
    _test_atom(target) do context
        y = constant(2 + 3im)
        return Variable() + imag(y)
    end
    target = """
    variables: x
    minobjective: 1.0 * x
    """
    _test_atom(target) do context
        y = constant(2)
        return Variable() + imag(y)
    end
    return
end

### affine/IndexAtom

function test_IndexAtom()
    target = """
    variables: x1, x2
    minobjective: [1.0 * x1, 1.0 * x2]
    """
    _test_atom(target) do context
        return Variable(2)[:, 1]
    end
    _test_atom(target) do context
        return Variable(2)[:, :]
    end
    target = """
    variables: x1, x2, x3
    minobjective: [1.0 * x1, 1.0 * x3]
    """
    _test_atom(target) do context
        return Variable(3)[[1, 3]]
    end
    target = """
    variables: x1, x2, x3
    minobjective: [1.0 * x1, 1.0 * x3]
    """
    _test_atom(target) do context
        return Variable(2, 2)[1, :]
    end
    target = """
    variables: x1, x2, x3, x4
    minobjective: [1.0 * x2, 1.0 * x4]
    """
    _test_atom(target) do context
        return Variable(2, 2)[2, :]
    end
    _test_atom(target) do context
        return Variable(2, 2)[2:2, :]
    end
    target = """
    variables: x1, x2
    minobjective: [1.0 * x1, 1.0 * x2]
    """
    _test_atom(target) do context
        return Variable(2, 2)[:, 1]
    end
    target = """
    variables: x1, x2, x3, x4
    minobjective: [1.0 * x3, 1.0 * x4]
    """
    _test_atom(target) do context
        return Variable(2, 2)[:, 2]
    end
    y = [true, false, true]
    x = Variable(3)
    @test string(x[y]) == string([x[1], x[3]])
    return
end

### affine/MultiplyAtom

function test_MultiplyAtom()
    target = """
    variables: x
    minobjective: 2.0 * x
    """
    _test_atom(target) do context
        return 2 * Variable()
    end
    _test_atom(target) do context
        return Variable() * 2
    end
    target = """
    variables: x, y
    minobjective: [2.0 * x, 2.0 * y]
    """
    _test_atom(target) do context
        return 2 * Variable(2)
    end
    _test_atom(target) do context
        return Variable(2) * 2
    end
    target = """
    variables: x11, x21, x12, x22
    minobjective: [2.0 * x11 + 3.0 * x12, 2.0 * x21 + 3.0 * x22]
    """
    _test_atom(target) do context
        return Variable(2, 2) * [2, 3]
    end
    target = """
    variables: x11, x21, x12, x22
    minobjective: [2.0 * x11 + 3.0 * x21, 2.0 * x12 + 3.0 * x22]
    """
    _test_atom(target) do context
        return [2, 3]' * Variable(2, 2)
    end
    target = """
    variables: t, x
    minobjective: 1.0 * t
    [1.0 + 1.0 * t, 1.0 + -1.0 * t, 2.0 * x] in SecondOrderCone(3)
    """
    _test_atom(target) do context
        x = Variable()
        return x * x
    end
    target = """
    variables: x
    minobjective: 0.5 * x
    """
    _test_atom(target) do context
        return Variable() / 2
    end
    target = """
    variables: x1, x2
    minobjective: [0.25 * x1, 0.25 * x2]
    """
    _test_atom(target) do context
        return Variable(2) / 4
    end
    _test_atom(target) do context
        return 0.25 .* Variable(2)
    end
    _test_atom(target) do context
        return Variable(2) .* 0.25
    end
    _test_atom(target) do context
        return Variable(2) ./ 4
    end
    target = """
    variables: x1, x2
    minobjective: [0.5 * x1, 0.25 * x2]
    """
    _test_atom(target) do context
        return Variable(2) ./ [2, 4]
    end
    @test_throws(
        ErrorException(
            "[MultiplyAtom] cannot multiply two expressions of sizes (2, 2) and (3, 3)",
        ),
        Variable(2, 2) * Variable(3, 3)
    )
    @test_throws(
        ErrorException(
            "[MultiplyAtom] multiplication of two non-constant expressions is not DCP compliant",
        ),
        _test_atom(_ -> Variable(2)' * Variable(2), ""),
    )
    @test_throws(
        ErrorException(
            "[MultiplyAtom] multiplication of two non-constant expressions is not DCP compliant",
        ),
        _test_atom(_ -> Variable(2) .* Variable(2), ""),
    )
    @test_throws(
        ErrorException(
            "[MultiplyAtom] multiplication of two non-constant expressions is not DCP compliant",
        ),
        _test_atom(_ -> Variable() * Variable(), ""),
    )
    return
end

### affine/NegateAtom

function test_NegateAtom()
    target = """
    variables: x
    minobjective: -1.0 + -1.0 * x
    """
    _test_atom(target) do context
        return -(1 + Variable())
    end
    target = """
    variables: x
    minobjective: 1.0 * x + -1.0
    """
    _test_atom(target) do context
        return Variable() + -constant(1.0)
    end
    return
end

### affine/RealAtom

function test_RealAtom()
    target = """
    variables: x
    minobjective: 1.0 * x
    """
    _test_atom(target) do context
        return real(Variable())
    end
    _test_atom(target) do context
        return real(Variable() + im * Variable())
    end
    target = """
    variables: x
    minobjective: 1.0 * x + 2.0
    """
    _test_atom(target) do context
        y = constant(2 + 3im)
        return Variable() + real(y)
    end
    target = """
    variables: x
    minobjective: 1.0 * x + 2.0
    """
    _test_atom(target) do context
        y = constant(2)
        return Variable() + real(y)
    end
    return
end

### affine/ReshapeAtom

function test_ReshapeAtom()
    target = """
    variables: x1, x2, x3, x4
    minobjective: [1.0 * x1, 1.0 * x2, 1.0 * x3, 1.0 * x4]
    """
    _test_atom(target) do context
        return reshape(Variable(4), 2, 2)
    end
    target = """
    variables: x1, x2, x3, x4
    minobjective: [1.0 * x1 + 2.0 * x3, 1.0 * x2 + 2.0 * x4]
    [-1.0 + x1, -2.0 + x2, -3.0 + x3, -4.0 + x4] in Nonnegatives(4)
    """
    _test_atom(target) do context
        x = Variable(4)
        add_constraint!(context, x - [1, 2, 3, 4] >= 0)
        return reshape(x, 2, 2) * [1, 2]
    end
    @test_throws(
        ErrorException(
            "[ReshapeAtom] cannot reshape expression of size (4, 1) to (2, 3)",
        ),
        reshape(Variable(4), 2, 3),
    )
    x = Variable(4)
    x.value = [1, 2, 3, 4]
    atom = reshape(x, 2, 2)
    @test evaluate(atom) == [1 3; 2 4]
    x = Variable()
    x.value = 2
    atom = reshape(x, 1, 1)
    @test evaluate(atom) == 2
    return
end

### affine/SumAtom

function test_SumAtom()
    target = """
    variables: x1, x2
    minobjective: 1.0 * x1 + 1.0 * x2
    """
    _test_atom(target) do context
        return sum(Variable(2))
    end
    target = """
    variables: x1, x2, x3, x4
    minobjective: 1.0 * x1 + 1.0 * x2 + 1.0 * x3 + 1.0 * x4
    """
    _test_atom(target) do context
        return sum(Variable(2, 2))
    end
    target = """
    variables: x1, x2, x3, x4
    minobjective: [1.0 * x1 + 1.0 * x2,  1.0 * x3 + 1.0 * x4]
    """
    _test_atom(target) do context
        return sum(Variable(2, 2); dims = 1)
    end
    target = """
    variables: x1, x2, x3, x4
    minobjective: [1.0 * x1 + 1.0 * x3,  1.0 * x2 + 1.0 * x4]
    """
    _test_atom(target) do context
        return sum(Variable(2, 2); dims = 2)
    end
    @test_throws(
        ErrorException("[SumAtom] sum not implemented for `dims=3`"),
        sum(Variable(2, 2); dims = 3),
    )
    return
end

### exp_+_sdp_cone/LogDetAtom

function test_LogDetAtom()
    target = """
    variables: t, x11, x12, x21, x22
    minobjective: 1.0 * t + 0.0
    [1.0*t, 1.0, 1.0*x11, 1.0 *x12, 1.0*x21, 1.0*x22] in LogDetConeSquare(2)
    """
    _test_atom(target) do context
        return logdet(Variable(2, 2))
    end
    return
end

### exp_cone/EntropyAtom

function test_EntropyAtom()
    target = """
    variables: t1, t2, x1, x2
    minobjective: 1.0 * t1 + 1.0 * t2
    [1.0 * t1, 1.0 * x1, 1.0] in ExponentialCone()
    [1.0 * t2, 1.0 * x2, 1.0] in ExponentialCone()
    """
    _test_atom(target) do context
        return entropy(Variable(2))
    end
    @test_throws(
        ErrorException(
            "[EntropyAtom] the argument should be real but it's instead complex",
        ),
        entropy(im * Variable(2)),
    )
    x = Variable(2)
    atom = entropy(x)
    x.value = [1.0 2.0]
    @test evaluate(atom) ≈ -sum(xi * log(xi) for xi in x.value)
    return
end

### exp_cone/ExpAtom

function test_ExpAtom()
    target = """
    variables: x, z
    minobjective: 1.0 * z + 0.0
    [1.0 * x, 1.0, 1.0 * z] in ExponentialCone()
    """
    _test_atom(target) do context
        return exp(Variable())
    end
    target = """
    variables: x1, x2, z1, z2
    minobjective: [1.0 * z1, 1.0 * z2]
    [1.0 * x1, 1.0, 1.0 * z1] in ExponentialCone()
    [1.0 * x2, 1.0, 1.0 * z2] in ExponentialCone()
    """
    _test_atom(target) do context
        return exp(Variable(2))
    end
    @test_throws(
        ErrorException(
            "[ExpAtom] the argument should be real but it's instead complex",
        ),
        exp(im * Variable()),
    )
    x = Variable(2)
    atom = exp(x)
    x.value = [1.0, 2.0]
    @test evaluate(atom) ≈ exp.([1.0, 2.0])
    return
end

### exp_cone/LogAtom

function test_LogAtom()
    target = """
    variables: x, z
    minobjective: 1.0 * x + 0.0
    [1.0 * x, 1.0, 1.0 * z] in ExponentialCone()
    """
    _test_atom(target) do context
        return log(Variable())
    end
    target = """
    variables: x1, x2, z1, z2
    minobjective: [1.0 * x1, 1.0 * x2]
    [1.0 * x1, 1.0, 1.0 * z1] in ExponentialCone()
    [1.0 * x2, 1.0, 1.0 * z2] in ExponentialCone()
    """
    _test_atom(target) do context
        return log(Variable(2))
    end
    @test_throws(
        ErrorException(
            "[LogAtom] the argument should be real but it's instead complex",
        ),
        log(im * Variable()),
    )
    x = Variable(2)
    atom = log(x)
    x.value = [1.0, 2.0]
    @test evaluate(atom) ≈ log.([1.0, 2.0])
    return
end

### exp_cone/LogSumExp

function test_LogSumExpAtom()
    target = """
    variables: x1, x2, t, z1, z2
    minobjective: 1.0 * t
    [1.0 * x1 + -1.0 * t, 1.0, 1.0 * z1] in ExponentialCone()
    [1.0 * x2 + -1.0 * t, 1.0, 1.0 * z2] in ExponentialCone()
    [1.0 + -1.0*z1 + -1.0*z2] in Nonnegatives(1)
    """
    _test_atom(target) do context
        return logsumexp(Variable(2))
    end
    target = """
    variables: x1, t, z1, z2
    minobjective: 1.0 * t
    [1.0 * x1 + -1.0 * t, 1.0, 1.0 * z1] in ExponentialCone()
    [-1.0 * t, 1.0, 1.0 * z2] in ExponentialCone()
    [1.0 + -1.0*z1 + -1.0*z2] in Nonnegatives(1)
    """
    _test_atom(target) do context
        return logisticloss(Variable())
    end
    target = """
    variables: x1, x1_, t, z1, z2, t_, z1_, z2_
    minobjective: 1.0 * t + 1.0 * t_
    [1.0 * x1 + -1.0 * t, 1.0, 1.0 * z1] in ExponentialCone()
    [-1.0 * t, 1.0, 1.0 * z2] in ExponentialCone()
    [1.0 * x1_ + -1.0 * t_, 1.0, 1.0 * z1_] in ExponentialCone()
    [-1.0 * t_, 1.0, 1.0 * z2_] in ExponentialCone()
    [1.0 + -1.0*z1 + -1.0*z2] in Nonnegatives(1)
    [1.0 + -1.0*z1_ + -1.0*z2_] in Nonnegatives(1)
    """
    _test_atom(target) do context
        return logisticloss(Variable(2))
    end
    @test_throws(
        ErrorException(
            "[LogSumExpAtom] the argument should be real but it's instead complex",
        ),
        logsumexp(im * Variable()),
    )
    x = Variable(2)
    atom = logsumexp(x)
    x.value = [1.0 1_000.0]
    @test evaluate(atom) ≈ 1_000.0
    return
end

### exp_cone/RelativeEntropyAtom

function test_RelativeEntropyAtom()
    target = """
    variables: w1, w2, v1, v2, u
    minobjective: 1.0 * u + 0.0
    [1.0*u, 1.0*v1, 1.0*v2, 1.0*w1, 1.0*w2] in RelativeEntropyCone(5)
    """
    _test_atom(target) do context
        x = Variable(2)
        y = Variable(2)
        return relative_entropy(x, y)
    end
    target = """
    variables: w1, w2, v1, v2, u
    minobjective: -1.0 * u + 1.0
    [1.0*u, 1.0*v1, 1.0*v2, 1.0*w1, 1.0*w2] in RelativeEntropyCone(5)
    """
    _test_atom(target) do context
        x = Variable(2)
        y = Variable(2)
        return 1.0 + log_perspective(x, y)
    end
    x, y = Variable(2), im * Variable(2)
    @test_throws(
        ErrorException(
            "[RelativeEntropyAtom] both the arguments should be real but these are instead $(sign(x)) and $(sign(y))",
        ),
        relative_entropy(x, y),
    )
    @test_throws(
        ErrorException(
            "[RelativeEntropyAtom] both the arguments should be real but these are instead $(sign(y)) and $(sign(x))",
        ),
        relative_entropy(y, x),
    )
    x = Variable(2)
    y = Variable(2)
    atom = relative_entropy(x, y)
    x.value = [1.0, 1_000.0]
    y.value = [1.0, NaN]
    @test evaluate(atom) ≈ Inf
    x.value = [0.0, 1.0]
    y.value = [1.0, 2.0]
    @test evaluate(atom) ≈ log(0.5)
    x.value = [5.0, 1.0]
    y.value = [3.0, 2.0]
    @test evaluate(atom) ≈ 5 * log(5 / 3) + log(0.5)
    return
end

### lp_cone/AbsAtom

function test_AbsAtom()
    target = """
    variables: t, x
    minobjective: 1.0 * t
    [1.0 * t + -1.0 * x] in Nonnegatives(1)
    [1.0 * t + 1.0 * x] in Nonnegatives(1)
    """
    _test_atom(target) do context
        return abs(Variable())
    end
    target = """
    variables: w, t, x
    minobjective: 1.0 * w
    [1.0 * t + -1.0 * x] in Nonnegatives(1)
    [1.0 * t + 1.0 * x] in Nonnegatives(1)
    [1.0 + 1.0*w, 1.0 + -1.0*w, 2.0*t] in SecondOrderCone(3)
    """
    _test_atom(target) do context
        return abs2(Variable())
    end
    target = """
    variables: t1, t2, x1, x2, w1, w2
    minobjective: [1.0 * t1, 1.0*t2]
    [1.0 * t1 + -1.0 * x1] in Nonnegatives(1)
    [1.0 * t2 + -1.0 * w2] in Nonnegatives(1)
    [1.0 * x1, 1.0 * x2, 2.0] in SecondOrderCone(3)
    [1.0 * w2, 1.0 * w1, 2.0] in SecondOrderCone(3)
    """
    _test_atom(target) do context
        return abs(Variable(2) + 2im)
    end
    return
end

### lp_cone/DotSortAtom

function test_DotSortAtom()
    target = """
    variables: v1, v2, u1, u2, x1, x2
    minobjective: 1.0 * u1 + u2 + v1 + v2
    [v1+u1, v1+u2, v2+u1+-1.0*x1, v2+u2+-1.0*x2] in Nonnegatives(4)
    """
    _test_atom(target) do context
        return dotsort(Variable(2), [0, 1])
    end
    target = """
    variables: v1, v2, u1, u2, x1, x2
    minobjective: 1.0 * u1 + u2 + v1 + v2
    [v1+u1+2.0*x1, v1+u2+2.0*x2, v2+u1+-3.0*x1, v2+u2+-3.0*x2] in Nonnegatives(4)
    """
    _test_atom(target) do context
        return dotsort(Variable(2), [-2, 3])
    end
    target = """
    variables: v1, v2, u1, u2, x1, x2
    minobjective: 1.0 * u1 + u2 + v1 + v2
    [v1+u1+2.0*x1, v1+u2+2.0*x2, v2+u1+3.0*x1, v2+u2+3.0*x2] in Nonnegatives(4)
    """
    _test_atom(target) do context
        return dotsort(Variable(2), [-2, -3])
    end
    x = im * Variable(2)
    @test_throws(
        ErrorException(
            "[DotSortAtom] argument should be real instead it is $(sign(x))",
        ),
        dotsort(x, [0, 1]),
    )
    @test_throws(
        ErrorException("[DotSortAtom] x and w must be the same size"),
        dotsort(Variable(2), [0, 1, 2]),
    )
    x = Variable(2, 2)
    atom = dotsort(x, [2 0; 0 1])
    x.value = [1.5 2.5; 3.0 2.0]
    @test evaluate(atom) ≈ 8.5
    x = Variable(4)
    atom = dotsort(x, [1 2; 0 0])
    x.value = [1.5, 2.5, 3.0, 2.0]
    @test evaluate(atom) ≈ 8.5
    x = Variable(2)
    atom = dotsort(square(x), [1, 2])
    @test curvature(atom) == Convex.ConvexVexity()
    @test sign(atom) == Convex.Positive()
    @test monotonicity(atom) == (Convex.Nondecreasing(),)
    atom = dotsort(square(x), [-1, -2])
    @test sign(atom) == Convex.Negative()
    @test monotonicity(atom) == (Convex.Nonincreasing(),)
    atom = dotsort(square(x), [-1, 2])
    @test sign(atom) == Convex.NoSign()
    @test monotonicity(atom) == (Convex.NoMonotonicity(),)
    return
end

### lp_cone/MaxAtom

function test_MaxAtom()
    target = """
    variables: t, x, y
    minobjective: 1.0 * t
    [1.0 * t + -1.0 * x] in Nonnegatives(1)
    [1.0 * t + -1.0 * y] in Nonnegatives(1)
    """
    _test_atom(target) do context
        return max(Variable(), Variable())
    end
    target = """
    variables: t, x, y
    minobjective: 1.0 * t
    [1.0 * t + -1.0 * x] in Nonnegatives(1)
    [1.0 * t + -1.0] in Nonnegatives(1)
    """
    _test_atom(target) do context
        return max(Variable(), 1)
    end
    _test_atom(target) do context
        return max(Variable(), constant(1))
    end
    target = """
    variables: t, x, y
    minobjective: 1.0 * t
    [1.0 * t + -1.0] in Nonnegatives(1)
    [1.0 * t + -1.0 * x] in Nonnegatives(1)
    """
    _test_atom(target) do context
        return max(1, Variable())
    end
    _test_atom(target) do context
        return max(constant(1), Variable())
    end
    target = """
    variables: t1, t2, x1, x2
    minobjective: [1.0 * t1, 1.0 * t2]
    [1.0 * t1 + -1.0 * x1, 1.0 * t2 + -1.0 * x2] in Nonnegatives(2)
    [1.0 * t1, 1.0 * t2] in Nonnegatives(2)
    """
    _test_atom(target) do context
        return pos(Variable(2))
    end
    target = """
    variables: t, x
    minobjective: 1.0 * t
    [-1.0 + 1.0 * t + 1.0 * x] in Nonnegatives(1)
    [1.0 * t] in Nonnegatives(1)
    """
    _test_atom(target) do context
        return hinge_loss(Variable())
    end
    x = Variable()
    @test sign(max(square(x), 1)) == Convex.Positive()
    @test sign(max(-square(x), -1)) == Convex.Negative()
    @test sign(max(x, -1)) == Convex.NoSign()
    @test sign(max(x, 1)) == Convex.Positive()
    return
end

### lp_cone/MaximumAtom

function test_MaximumAtom()
    target = """
    variables: t, x, y
    minobjective: 1.0 * t
    [1.0 * t + -1.0 * x, 1.0 * t + -1.0 * y] in Nonnegatives(2)
    """
    _test_atom(target) do context
        return maximum(Variable(2))
    end
    x = im * Variable(2)
    @test_throws(
        ErrorException(
            "[MaximumAtom] argument should be real instead it is $(sign(x))",
        ),
        maximum(x),
    )
    return
end

### lp_cone/MinAtom

function test_MinAtom()
    target = """
    variables: t, x, y
    minobjective: 1.0 * t
    [-1.0 * t + 1.0 * x] in Nonnegatives(1)
    [-1.0 * t + 1.0 * y] in Nonnegatives(1)
    """
    _test_atom(target) do context
        return min(Variable(), Variable())
    end
    target = """
    variables: t, x, y
    minobjective: 1.0 * t
    [-1.0 * t + 1.0 * x] in Nonnegatives(1)
    [-1.0 * t + 1.0] in Nonnegatives(1)
    """
    _test_atom(target) do context
        return min(Variable(), 1)
    end
    _test_atom(target) do context
        return min(Variable(), constant(1))
    end
    target = """
    variables: t, x, y
    minobjective: 1.0 * t
    [-1.0 * t + 1.0] in Nonnegatives(1)
    [-1.0 * t + 1.0 * x] in Nonnegatives(1)
    """
    _test_atom(target) do context
        return min(1, Variable())
    end
    _test_atom(target) do context
        return min(constant(1), Variable())
    end
    target = """
    variables: t1, t2, x1, x2
    minobjective: [1.0 * t1, 1.0 * t2]
    [1.0 * t1 + 1.0 * x1, 1.0 * t2 + 1.0 * x2] in Nonnegatives(2)
    [1.0 * t1, 1.0 * t2] in Nonnegatives(2)
    """
    _test_atom(target) do context
        return neg(Variable(2))
    end
    x = Variable()
    @test sign(min(square(x), 1)) == Convex.Positive()
    @test sign(min(-square(x), -1)) == Convex.Negative()
    @test sign(min(x, -1)) == Convex.Negative()
    @test sign(min(x, 1)) == Convex.NoSign()
    x, y = Variable(2), Variable(3)
    @test_throws(
        ErrorException(
            "[MinAtom] got different sizes for x as $(x.size) and y as $(y.size)",
        ),
        min(x, y),
    )
    x = im * x
    @test_throws(
        ErrorException(
            "[MinAtom] both the arguments should be real instead they are $(sign(x)) and $(sign(y))",
        ),
        min(x, y),
    )
    return
end

### lp_cone/MinimumAtom

function test_MimimumAtom()
    target = """
    variables: x, y, t
    minobjective: 1.0 * t
    [1.0 * x + -1.0 * t, 1.0 * y + -1.0 * t] in Nonnegatives(2)
    """
    _test_atom(target) do context
        return minimum(Variable(2))
    end
    x = im * Variable(2)
    @test_throws(
        ErrorException(
            "[MinimumAtom] argument should be real instead it is $(sign(x))",
        ),
        minimum(x),
    )
    return
end

### lp_cone/SumLargestAtom

function test_SumLargestAtom()
    target = """
    variables: t1, t2, t3, q, x1, x2, x3
    minobjective: 1.0 * t1 + 1.0 * t2 + 1.0 * t3 + 2.0 * q
    [1.0*t1 + 1.0*q + -1.0*x1, 1.0*t2 + 1.0*q + -1.0*x2, 1.0*t3 + 1.0*q + -1.0*x3] in Nonnegatives(3)
    [1.0 * t1, 1.0 * t2, 1.0 * t3] in Nonnegatives(3)
    """
    _test_atom(target) do context
        return sumlargest(Variable(3), 2)
    end
    target = """
    variables: t1, t2, t3, q, x1, x2, x3
    minobjective: 1.0 + -1.0 * t1 + -1.0 * t2 + -1.0 * t3 + -2.0 * q
    [1.0*t1 + 1.0*q + 1.0*x1, 1.0*t2 + 1.0*q + 1.0*x2, 1.0*t3 + 1.0*q + 1.0*x3] in Nonnegatives(3)
    [1.0 * t1, 1.0 * t2, 1.0 * t3] in Nonnegatives(3)
    """
    _test_atom(target) do context
        return 1 + sumsmallest(Variable(3), 2)
    end
    @test string(sumlargest(Variable(3), 0)) == "0"
    x = im * Variable(3)
    @test_throws(
        ErrorException(
            "[SumLargestAtom] argument should be real instead it is $(sign(x))",
        ),
        sumlargest(x, 2),
    )
    @test_throws(
        ErrorException(
            "[SumLargestAtom] sumlargest and sumsmallest only support positive values of k",
        ),
        sumlargest(Variable(3), -2),
    )
    @test_throws(
        ErrorException(
            "[SumLargestAtom] k cannot be larger than the number of entries in x",
        ),
        sumlargest(Variable(3), 4),
    )
    return
end

### sdp_cone/EigMaxAtom

function test_EigMaxAtom()
    target = """
    variables: t, x11, x21, x12, x22
    minobjective: 1.0 * t
    [1.0 * t + -1.0 * x11, -1.0*x21, -1.0*x12, 1.0 * t + -1.0 * x22] in PositiveSemidefiniteConeSquare(2)
    """
    _test_atom(target) do context
        return eigmax(Variable(2, 2))
    end
    @test_throws(
        ErrorException(
            "[EigMaxAtom] eigmax can only be applied to a square matrix.",
        ),
        eigmax(Variable(2, 3)),
    )
    return
end

### sdp_cone/EigMinAtom

function test_EigMinAtom()
    target = """
    variables: x11, x21, x12, x22, t
    minobjective: 1.0 * t
    [-1.0 * t + 1.0 * x11, 1.0*x21, 1.0*x12, -1.0 * t + 1.0 * x22] in PositiveSemidefiniteConeSquare(2)
    """
    _test_atom(target) do context
        return eigmin(Variable(2, 2))
    end
    @test_throws(
        ErrorException(
            "[EigMinAtom] eigmin can only be applied to a square matrix.",
        ),
        eigmin(Variable(2, 3)),
    )
    return
end

### sdp_cone/MatrixFracAtom

function test_MatrixFracAtom()
    target = """
    variables: x1, x2, t, P11, P21, P12, P22
    minobjective: 1.0 * t
    [1.0 * x1, 1.0 * x2] in Nonnegatives(2)
    [1.0*t, 1.0*x1, 1.0*x2, 1.0*x1, 1.0*P11, 1.0*P21, 1.0*x2, 1.0*P12, 1.0*P22] in PositiveSemidefiniteConeSquare(3)
    """
    _test_atom(target) do context
        x = Variable(2)
        add_constraint!(context, x >= 0)
        return matrixfrac(x, Variable(2, 2))
    end
    target = """
    variables: t, P11, P21, P12, P22
    minobjective: 1.0 * t
    [1.0*t, 1.0, 2.0, 1.0, 1.0*P11, 1.0*P21, 2.0, 1.0*P12, 1.0*P22] in PositiveSemidefiniteConeSquare(3)
    """
    _test_atom(target) do context
        return matrixfrac([1, 2], Variable(2, 2))
    end
    target = """
    variables: t, x1, x2
    minobjective: 1.0 * t
    [1.0*t, 1.0*x1, 1.0*x2, 1.0*x1, 2.0, 0.0, 1.0*x2, 0.0, 3.0] in PositiveSemidefiniteConeSquare(3)
    """
    _test_atom(target) do context
        return matrixfrac(Variable(2), [2 0; 0 3])
    end
    @test_throws(
        ErrorException(
            "[MatrixFracAtom] first argument of matrixfrac must be a vector",
        ),
        matrixfrac(Variable(2, 2), Variable(2, 2)),
    )
    @test_throws(
        ErrorException(
            "[MatrixFracAtom] second argument of matrixfrac must be square",
        ),
        matrixfrac(Variable(2), Variable(2, 3)),
    )
    @test_throws(
        ErrorException(
            "[MatrixFracAtom] sizes must agree for arguments of matrixfrac",
        ),
        matrixfrac(Variable(2), Variable(3, 3)),
    )
    return
end

### sdp_cone/NuclearNormAtom

# TODO

### sdp_cone/OperatorNormAtom

# TODO

### sdp_cone/QuantumEntropyAtom

# TODO

### sdp_cone/QuantumRelativeEntropyAtom

# TODO

### sdp_cone/SumLargestEigsAtom

# TODO

### sdp_cone/TraceLogmAtom

# TODO

### sdp_cone/TraceMpowerAtom

# TODO

### second_order_cone/EuclideanNormAtom

function test_EuclideanNormAtom()
    target = """
    variables: t, x, y
    minobjective: 1.0 * t
    [1.0 * t, 1.0 * x, 1.0 * y] in SecondOrderCone(3)
    """
    _test_atom(target) do context
        return norm2(Variable(2))
    end
    target = """
    variables: t, x, y
    minobjective: 1.0 * t
    [1.0 * t, 1.0 * x, 2.0 * x] in SecondOrderCone(3)
    """
    _test_atom(target) do context
        x = Variable()
        return norm2((1 + 2im) * x)
    end
    return
end

### second_order_cone/GeoMeanAtom

function test_GeoMeanAtom()
    target = """
    variables: t, x, y
    minobjective: 1.0 * t
    [1.0 * t, 1.0 * x, 1.0 * y] in GeometricMeanCone(3)
    """
    _test_atom(target) do context
        return geomean(Variable(), Variable())
    end
    target = """
    variables: t1, t2, x1, x2, y1, y2
    minobjective: [1.0 * t1, 1.0 * t2]
    [1.0 * t1, 1.0 * x1, 1.0 * y1] in GeometricMeanCone(3)
    [1.0 * t2, 1.0 * x2, 1.0 * y2] in GeometricMeanCone(3)
    """
    _test_atom(target) do context
        return geomean(Variable(2), Variable(2))
    end
    @test_throws(
        ErrorException(
            "[GeoMeanAtom] geomean must take arguments of the same size",
        ),
        geomean(Variable(2), Variable(3)),
    )
    @test_throws(
        ErrorException("[GeoMeanAtom] the arguments must be real, not complex"),
        geomean(Variable(), 2im),
    )
    x = Variable(2)
    x.value = [2.0, 3.0]
    y = Variable(2)
    y.value = [4.0, 5.0]
    atom = geomean(x, y)
    @test evaluate(atom) ≈ [sqrt(8), sqrt(15)]
    return
end

### second_order_cone/HuberAtom

function test_HuberAtom()
    target = """
    variables: c, s, n, t, n_abs
    minobjective: 1.0 * t + 4.0 * n_abs
    [1.0 * c + -1.0 * s + -1.0 * n] in Zeros(1)
    [-1.0 * n + 1.0 * n_abs] in Nonnegatives(1)
    [1.0 * n + 1.0 * n_abs] in Nonnegatives(1)
    [1.0 + 1.0 * t, 1.0 + -1.0 * t, 2.0 * s] in SecondOrderCone(3)
    """
    _test_atom(target) do context
        return huber(Variable(), 2.0)
    end
    @test_throws(
        ErrorException(
            "[HuberAtom] parameter must by a positive scalar. Got `M=-2.0`",
        ),
        huber(Variable(2), -2.0),
    )
    @test_throws(
        ErrorException("[HuberAtom] argument must be real"),
        huber(im * Variable(2)),
    )
    x = Variable(2)
    x.value = [2.0, 3.0]
    atom = huber(x, 2.0)
    @test evaluate(atom) ≈ [4.0, 8.0]
    atom = huber(x)
    @test evaluate(atom) ≈ [3.0, 5.0]
    return
end

### second_order_cone/QolElemAtom

function test_QolElemAtom()
    target = """
    variables: y1, y2, t1, t2, x1, x2
    minobjective: [1.0 * t1, 1.0 * t2]
    [1.0 * y1, 1.0 * y2] in Nonnegatives(2)
    [1.0 * y1 + 1.0 * t1, 1.0 * y1 + -1.0 * t1, 2.0 * x1] in SecondOrderCone(3)
    [1.0 * y2 + 1.0 * t2, 1.0 * y2 + -1.0 * t2, 2.0 * x2] in SecondOrderCone(3)
    """
    _test_atom(target) do context
        return qol_elementwise(Variable(2), Variable(2))
    end
    target = """
    variables: t1, t2, x1, x2
    minobjective: [1.0 * t1, 1.0 * t2]
    [1.0 + 1.0 * t1, 1.0 + -1.0 * t1, 2.0 * x1] in SecondOrderCone(3)
    [1.0 + 1.0 * t2, 1.0 + -1.0 * t2, 2.0 * x2] in SecondOrderCone(3)
    """
    _test_atom(target) do context
        return qol_elementwise(Variable(2), constant([1, 1]))
    end
    _test_atom(target) do context
        return square(Variable(2))
    end
    _test_atom(target) do context
        return Variable(2) .^ 2
    end
    _test_atom(target) do context
        a = 2
        return Variable(2) .^ a
    end
    target = """
    variables: y1, y2, t1, t2
    minobjective: [1.0 * t1, 1.0 * t2]
    [1.0 * y1, 1.0 * y2] in Nonnegatives(2)
    [1.0 * y1 + 1.0 * t1, 1.0 * y1 + -1.0 * t1, 2.0] in SecondOrderCone(3)
    [1.0 * y2 + 1.0 * t2, 1.0 * y2 + -1.0 * t2, 2.0] in SecondOrderCone(3)
    """
    _test_atom(target) do context
        return invpos(Variable(2))
    end
    _test_atom(target) do context
        return 1 ./ Variable(2)
    end
    target = """
    variables: y, t
    minobjective: 3.0 * t
    [1.0 * y] in Nonnegatives(1)
    [1.0 * y + 1.0 * t, 1.0 * y + -1.0 * t, 2.0] in SecondOrderCone(3)
    """
    _test_atom(target) do context
        return 3 / Variable()
    end
    y = Variable(2)
    @test_throws(
        ErrorException("cannot divide by a variable of size $(size(y))"),
        1 / y,
    )

    @test_throws(
        ErrorException(
            "square of a complex number is not DCP. Did you mean square_modulus?",
        ),
        square(im * Variable()),
    )
    @test_throws(
        ErrorException(
            "raising variables to powers other than 2 is not implemented",
        ),
        Variable() .^ 3,
    )
    @test_throws(
        ErrorException(
            "[QolElemAtom] elementwise quad over lin must take two arguments of the same size",
        ),
        qol_elementwise(Variable(2), Variable())
    )
    x = Variable(2)
    x.value = [2.0, 3.0]
    y = Variable(2)
    y.value = [4.0, 5.0]
    atom = qol_elementwise(x, y)
    @test evaluate(atom) ≈ [1.0, 9.0 / 5.0]
    return
end

### second_order_cone/QuadOverLinAtom

function test_QuadOverLinAtom()
    target = """
    variables: y, t, x1, x2
    minobjective: 1.0 * t
    [1.0 * y] in Nonnegatives(1)
    [1.0 * y + 1.0 * t, 1.0 * y + -1.0 * t, 2.0 * x1, 2.0 * x2] in SecondOrderCone(4)
    """
    _test_atom(target) do context
        return quadoverlin(Variable(2), Variable())
    end
    @test_throws(
        ErrorException(
            "[QuadOverLinAtom] quadoverlin arguments must be a vector and a scalar",
        ),
        quadoverlin(Variable(2), Variable(2))
    )
    x = Variable(2)
    x.value = [2.0, 3.0]
    atom = quadoverlin(x, constant(2.0))
    @test evaluate(atom) ≈ 13 / 2
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
