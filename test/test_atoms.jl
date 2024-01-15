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

# TODO

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

# TODO

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

# TODO

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
