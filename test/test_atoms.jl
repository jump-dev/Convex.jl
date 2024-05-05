# Copyright (c) 2014: Madeleine Udell and contributors
#
# Use of this source code is governed by a BSD-style license that can be found
# in the LICENSE file or at https://opensource.org/license/bsd-2-clause

module TestAtoms

using Convex
using Test

# Do not use `using LinearAlgebra` to check symbols are re-exported by Convex.
import LinearAlgebra
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
    _test_atom(build_fn::Function, target_string::String; value_type = Float64)

The same arguments and behavior as `_test_reformulation`, but in addition, tests
a number of properties that all atoms must satisfy.

 1. the atom is a subtype of `AbstractExpr`
 2. the atom is mutable
 3. `Convex.head` is implemented and prints a string
 4. `Base.sign` is implemented and returns a `Convex.Sign` object
 5. `Covnex.monotonicity` is implemented and returns a tuple of
    `Convex.Monotonicity` objects, with one element for each child
 6. `Convex.curvature` is implemented and returns a `Convex.Vexity` object

 ## Example

 ```julia
 target = \"\"\"
 variables: x
 minobjective: 1.0 + 1.0 * x
 \"\"\"
 _test_atom(target) do context
     return Variable() + 1
 end
 ```
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
    # Now, we will test if our atom's `conic_form` and `evaluate`
    # agree for constants (and that they are well-defined for constants).
    # To do so, we first flip our atom to have non-concave vexity:
    if Convex.vexity(atom) == Convex.ConcaveVexity()
        atom = -atom
    end
    _test_reformulation(build_fn, target_string; value_type)
    _test_constant_atom(build_fn; value_type)
    return
end

# Recursively _make_constant! an `AbstractExpr` to replace variables with
# constants
function _to_constant(expr::Convex.AbstractVariable)
    if Convex.iscomplex(expr)
        val = rand(Float64, size(expr)) + im * rand(Float64, size(expr))
    else
        val = rand(Float64, size(expr))
    end
    if size(expr, 1) == size(expr, 2)
        val = val' * val
    end
    return constant(val)
end

_to_constant(x::Convex.Value) = x
_to_constant(x::Union{Convex.Constant,Convex.ComplexConstant}) = x

function _to_constant(e::Convex.AbstractExpr)
    e.children = map(_to_constant, e.children)
    return e
end

function _test_constant_atom(build_fn; value_type)
    context = Convex.Context{value_type}(MOI.Utilities.Model{value_type})
    atom = _to_constant(build_fn(context))
    if Convex.iscomplex(atom) || any(Convex.iscomplex, atom.children)
        return
    end
    form = Convex.conic_form!(context, atom)
    if !(form isa AbstractVector)
        return  # The reformulation is still in terms of MOI variables.
    end
    answer = evaluate(atom)
    if answer isa Number
        @test only(form) ≈ answer
    else
        @test form ≈ vec(answer)
    end
    return
end

"""
    _test_reformulation(
        build_fn::Function,
        target_string::String;
        value_type = Float64,
    )

Test the reformulation constructed by `build_fn(::Convex.Context)` produces an
optimization problem given by `target_string` when the result returned by
`build_fn` is set as the objective function.

## Arguments

 * `build_fn`: a function called with a new `context::Convex.Context` object
   that returns an expression for the objective and may optionally modify
   `context` as well.
 * `target_string`: the string representation of the model as needed by
   `MOI.Utilities.loadfromstring!`.
 * `value_type`: the numeric value type of the model.

## Example

```julia
target = \"\"\"
variables: x
minobjective: 1.0 + 1.0 * x
\"\"\"
_test_reformulation(target) do context
    return Variable() + 1
end
```
"""
function _test_reformulation(
    build_fn,
    target_string::String;
    value_type = Float64,
)
    context = Convex.Context{value_type}(MOI.Utilities.Model{value_type})
    atom = build_fn(context)
    t = Convex.conic_form!(context, atom)
    MOI.set(context.model, MOI.ObjectiveSense(), MOI.MIN_SENSE)
    obj = _to_moi(t)
    @test prod(size(atom)) == MOI.output_dimension(obj)
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
    target = """
    variables: x1, x2
    minobjective: [1.0 * x1, 1.0 * x2, 2.0]
    """
    _test_atom(target) do context
        x = Variable(1, 2)
        y = constant(2)
        return hcat(x, y)
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
    [t, 0.5, x] in RotatedSecondOrderCone(3)
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

### affine/VcatAtom

function test_VcatAtom()
    target = """
    variables: x
    minobjective: [1.0 * x, 1.0 * x]
    """
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
    target = """
    variables: x1, x2
    minobjective: [1.0 * x1, 2.0, 1.0 * x2, 3.0]
    """
    _test_atom(target) do context
        x = Variable(1, 2)
        y = constant([2 3])
        return vcat(x, y)
    end
    target = """
    variables: x1, x2, x3
    minobjective: [2.0, 1.0 * x1, 2.0, 3.0, 1.0 * x2, 3.0, 4.0, 1.0 * x3, 4.0]
    """
    _test_atom(target) do context
        x = Variable(1, 3)
        y = constant([2 3 4])
        return vcat(y, x, y)
    end
    target = """
    variables: x1, x2, x3, x4
    minobjective: [x1, x2, 2.0, x3, x4, 3.0]
    """
    _test_atom(target) do context
        x = Variable(2, 2)
        y = constant([2 3])
        return vcat(x, y)
    end
    @test_throws(
        DimensionMismatch(
            "[VcatAtom] cannot stack expressions of incompatible size. Got 2 expected 1.",
        ),
        vcat(Variable(2, 1), Variable(1, 2)),
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

function test_EntropyAtom_elementwise()
    target = """
    variables: t1, t2, x1, x2
    minobjective: [1.0 * t1, + 1.0 * t2]
    [1.0 * t1, 1.0 * x1, 1.0] in ExponentialCone()
    [1.0 * t2, 1.0 * x2, 1.0] in ExponentialCone()
    """
    _test_atom(target) do context
        return entropy_elementwise(Variable(2))
    end
    x = Variable(2)
    atom = entropy_elementwise(x)
    x.value = [1.0, 2.0]
    @test evaluate(atom) ≈ -x.value .* log.(x.value)
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
    variables: u, v1, v2, w1, w2
    minobjective: 1.0 * u + 0.0
    [1.0*u, 1.0*v1, 1.0*v2, 1.0*w1, 1.0*w2] in RelativeEntropyCone(5)
    """
    _test_atom(target) do context
        x = Variable(2)
        y = Variable(2)
        return relative_entropy(x, y)
    end
    target = """
    variables: u, v1, v2, w1, w2
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
    [w, 0.5, t] in RotatedSecondOrderCone(3)
    """
    _test_atom(target) do context
        return abs2(Variable())
    end
    target = """
    variables: t1, t2, x1, x2
    minobjective: [1.0 * t1, 1.0*t2]
    [1.0 * t1, 1.0 * x1, 2.0] in SecondOrderCone(3)
    [1.0 * t2, 1.0 * x2, 2.0] in SecondOrderCone(3)
    """
    _test_atom(target) do context
        return abs(Variable(2) + 2im)
    end
    target = """
    variables: y1, y2
    minobjective: [2 + y1, 5 + y2]
    """
    _test_atom(target) do context
        x = ComplexVariable(2)
        y = Variable(2)
        fix!(x, [2, 3 - 4im])
        return y + abs(x)
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
    _test_atom(target) do context
        return dotsort([0, 1], Variable(2))
    end
    target = """
    variables: v1, v2, u1, u2, x1, x2
    minobjective: 1.0 * u1 + u2 + v1 + v2
    [v1+u1+2.0*x1, v1+u2+2.0*x2, v2+u1+-3.0*x1, v2+u2+-3.0*x2] in Nonnegatives(4)
    """
    _test_atom(target) do context
        return dotsort(Variable(2), [-2, 3])
    end
    _test_atom(target) do context
        return dotsort([-2, 3], Variable(2))
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
    variables: t1, t2, x, y, z
    minobjective: [1.0 * t1, 1.0 * t2]
    [1.0 * t1 + -1.0 * x, t2 + -1.0 * y] in Nonnegatives(2)
    [1.0 * t1 + -1.0 * z, t2 + -1.0 * z] in Nonnegatives(2)
    """
    _test_atom(target) do context
        return max(Variable(2), Variable())
    end
    target = """
    variables: t1, t2, x, y, z
    minobjective: [1.0 * t1, 1.0 * t2]
    [1.0 * t1 + -1.0 * x, t2 + -1.0 * x] in Nonnegatives(2)
    [1.0 * t1 + -1.0 * y, t2 + -1.0 * z] in Nonnegatives(2)
    """
    _test_atom(target) do context
        return max(Variable(), Variable(2))
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
    x, y = Variable(2), Variable(3)
    @test_throws(
        ErrorException(
            "[MaxAtom] got different sizes for x as $(x.size) and y as $(y.size)",
        ),
        max(x, y),
    )
    x = im * x
    @test_throws(
        ErrorException(
            "[MaxAtom] both the arguments should be real instead they are $(sign(x)) and $(sign(y))",
        ),
        max(x, y),
    )
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
    variables: t1, t2, x, y, z
    minobjective: [1.0 * t1, 1.0 * t2]
    [-1.0 * t1 + x, -1.0 * t2 + y] in Nonnegatives(2)
    [-1.0 * t1 + z, -1.0 * t2 + z] in Nonnegatives(2)
    """
    _test_atom(target) do context
        return min(Variable(2), Variable())
    end
    target = """
    variables: t1, t2, x, y, z
    minobjective: [1.0 * t1, 1.0 * t2]
    [-1.0 * t1 + x, -1.0 * t2 + x] in Nonnegatives(2)
    [-1.0 * t1 + y, -1.0 * t2 + z] in Nonnegatives(2)
    """
    _test_atom(target) do context
        return min(Variable(), Variable(2))
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

function test_TraceLogmAtom()
    @test_throws(
        DimensionMismatch("X and C must be the same size"),
        trace_logm(Variable(2, 3), [1 0; 0 1]),
    )
    @test_throws(
        DimensionMismatch("X and C must be square"),
        trace_logm(Variable(2, 3), [1 2 3; 4 5 6]),
    )
    C = [1 2; 5 6]
    @test_throws(
        DomainError(C, "C must be Hermitian"),
        trace_logm(Variable(2, 2), C),
    )
    C = [-1 0; 0 2]
    @test_throws(
        DomainError(C, "C must be positive semidefinite"),
        trace_logm(Variable(2, 2), C),
    )

    C = [1 0; 0 1]
    X = Variable(2, 2)
    atom = trace_logm(X, C)
    @test curvature(atom) == Convex.ConcaveVexity()
    X.value = [2 1; 1 2]
    @test evaluate(atom) ≈ LinearAlgebra.tr(C * log(X.value))
    # TODO(odow): add a test for the reformulation
    return
end

### sdp_cone/TraceMpowerAtom

function test_TraceMpowerAtom()
    @test_throws(
        DimensionMismatch("A and C must be the same size"),
        trace_mpower(Variable(2, 3), 1 // 2, [1 0; 0 1]),
    )
    @test_throws(
        DimensionMismatch("A and C must be square"),
        trace_mpower(Variable(2, 3), 1 // 2, [1 2 3; 4 5 6]),
    )
    C = [1 2; 5 6]
    @test_throws(
        DomainError(C, "C must be Hermitian"),
        trace_mpower(Variable(2, 2), 1 // 2, C),
    )
    C = [-1 0; 0 2]
    @test_throws(
        DomainError(C, "C must be positive semidefinite"),
        trace_mpower(Variable(2, 2), 1 // 2, C),
    )
    @test_throws(
        DomainError(-2 // 1, "t must be in the range [-1, 2]"),
        trace_mpower(Variable(2, 2), -2 // 1, [1 0; 0 1]),
    )
    t = -1 // 2
    C = [1 0; 0 1]
    A = Variable(2, 2)
    atom = trace_mpower(A, t, C)
    @test curvature(trace_mpower(A, 0 // 1, C)) == Convex.ConcaveVexity()
    @test curvature(trace_mpower(A, 1 // 2, C)) == Convex.ConcaveVexity()
    @test curvature(trace_mpower(A, 1 // 1, C)) == Convex.ConcaveVexity()
    @test curvature(trace_mpower(A, -1 // 1, C)) == Convex.ConvexVexity()
    @test curvature(trace_mpower(A, -1 // 2, C)) == Convex.ConvexVexity()
    @test curvature(trace_mpower(A, 3 // 2, C)) == Convex.ConvexVexity()
    A.value = [2 1; 1 2]
    @test evaluate(atom) ≈ LinearAlgebra.tr(C * A.value^t)
    # TODO(odow): add a test for the reformulation
    return
end

### sdp_cone/RootDetAtom

function test_RootDetAtom()
    target = """
    variables: t, x11, x12, x21, x22
    minobjective: 1.0 * t + 0.0
    [1.0*t, 1.0*x11, 1.0 *x12, 1.0*x21, 1.0*x22] in RootDetConeSquare(2)
    """
    _test_atom(target) do context
        return rootdet(Variable(2, 2))
    end
    x = Variable(2, 2)
    x.value = [2.0 -1.5; -1.5 3.0]
    atom = rootdet(x)
    @test evaluate(atom) ≈ LinearAlgebra.det(x.value)^(1 / 2)
    return
end

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
    [t, 0.5, s] in RotatedSecondOrderCone(3)
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
    variables: t1, t2, y1, y2, x1, x2
    minobjective: [1.0 * t1, 1.0 * t2]
    [t1, 0.5 * y1, x1] in RotatedSecondOrderCone(3)
    [t2, 0.5 * y2, x2] in RotatedSecondOrderCone(3)
    """
    _test_atom(target) do context
        return qol_elementwise(Variable(2), Variable(2))
    end
    target = """
    variables: t1, t2, x1, x2
    minobjective: [1.0 * t1, 1.0 * t2]
    [t1, 0.5, x1] in RotatedSecondOrderCone(3)
    [t2, 0.5, x2] in RotatedSecondOrderCone(3)
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
    variables: t1, t2, y1, y2
    minobjective: [1.0 * t1, 1.0 * t2]
    [t1, 0.5 * y1, 1.0] in RotatedSecondOrderCone(3)
    [t2, 0.5 * y2, 1.0] in RotatedSecondOrderCone(3)
    """
    _test_atom(target) do context
        return invpos(Variable(2))
    end
    _test_atom(target) do context
        return 1 ./ Variable(2)
    end
    target = """
    variables: t, y
    minobjective: 3.0 * t
    [t, 0.5 * y, 1.0] in RotatedSecondOrderCone(3)
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
    variables: t, y, x1, x2
    minobjective: 1.0 * t
    [t, 0.5 * y, x1, x2] in RotatedSecondOrderCone(4)
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

### reformulations/conv

function test_conv()
    target = """
    variables: x1, x2
    minobjective: [1.0 * x1, 2.0 * x1 + 1.0 * x2, 2.0 * x2]
    """
    _test_reformulation(target) do context
        return conv(Variable(2), [1, 2])
    end
    _test_reformulation(target) do context
        return conv([1, 2], Variable(2))
    end
    target = """
    variables: x1, x2, x3
    minobjective: [1.0*x1, -3.0*x1+1.0*x2, 2.0*x1+-3.0*x2+1.0*x3, 2.0*x2+-3.0*x3, 2.0*x3]
    """
    _test_reformulation(target) do context
        return conv(Variable(3), [1, -3, 2])
    end
    @test_throws(
        ErrorException("convolution only supported between two vectors"),
        conv([1 2; 3 4], Variable(2, 2)),
    )
    @test_throws(
        ErrorException("convolution only supported between two vectors"),
        conv([1, 2], Variable(2, 2)),
    )
    return
end

### reformulations/dot

function test_dot()
    target = """
    variables: x1, x2
    minobjective: 1.0 * x1 + 2.0 * x2
    """
    _test_reformulation(target) do context
        return dot(Variable(2), [1, 2])
    end
    _test_reformulation(target) do context
        return dot([1, 2], Variable(2))
    end
    _test_reformulation(target) do context
        return dot(constant([1, 2]), Variable(2))
    end
    _test_reformulation(target) do context
        return dot(Variable(2), constant([1, 2]))
    end
    return
end

### reformulations/inner_product

function test_inner_product()
    target = """
    variables: x
    minobjective: 1.0 * x
    """
    _test_reformulation(target) do context
        return inner_product(Variable(1), [1;;])
    end
    _test_reformulation(target) do context
        return inner_product([1;;], Variable(1))
    end
    _test_reformulation(target) do context
        return inner_product(constant([1;;]), Variable(1))
    end
    target = """
    variables: x11, x21, x12, x22
    minobjective: 1.0 * x11 + 3.0 * x21 + 2.0 * x12 + 4.0 * x22
    """
    _test_reformulation(target) do context
        return inner_product(Variable(2, 2), [1 2; 3 4])
    end
    _test_reformulation(target) do context
        return inner_product([1 2; 3 4], Variable(2, 2))
    end
    @test_throws(
        ErrorException(
            "arguments must be square matrices of the same dimension",
        ),
        inner_product([1;;], Variable(2)),
    )
    @test_throws(
        ErrorException(
            "arguments must be square matrices of the same dimension",
        ),
        inner_product([1 2; 3 4], Variable(3, 3)),
    )
    return
end

### reformulations/norm

function test_norm()
    target = """
    variables: t1, t2, x1, x2
    minobjective: 1.0 * t1 + t2
    [1.0*t1+-1.0*x1, 1.0*t2+-1.0*x2] in Nonnegatives(2)
    [1.0*t1+1.0*x1, 1.0*t2+1.0*x2] in Nonnegatives(2)
    """
    _test_reformulation(target) do context
        return norm(Variable(2), 1)
    end
    target = """
    variables: t1, t2, x1, x2
    minobjective: 1.0 * t1 + t2
    [1.0*t1+-1.0*x1, 1.0*t2+-1.0*x2] in Nonnegatives(2)
    [1.0*t1+1.0*x1, 1.0*t2+1.0*x2] in Nonnegatives(2)
    """
    _test_reformulation(target) do context
        return norm(Variable(1, 2), 1)
    end
    target = """
    variables: t, x1, x2
    minobjective: 1.0 * t
    [1.0*t, 1.0*x1, 1.0*x2] in SecondOrderCone(3)
    """
    _test_reformulation(target) do context
        return norm(Variable(2), 2)
    end
    target = """
    variables: t, t1, t2, x1, x2
    minobjective: 1.0 * t
    [1.0*t1+-1.0*x1, 1.0*t2+-1.0*x2] in Nonnegatives(2)
    [1.0*t1+1.0*x1, 1.0*t2+1.0*x2] in Nonnegatives(2)
    [1.0*t + -1.0t1, 1.0*t + -1.0*t2] in Nonnegatives(2)
    """
    _test_reformulation(target) do context
        return norm(Variable(2), Inf)
    end
    target = """
    variables: x1, x2, t
    minobjective: 1.0 * t
    [1.0*t, 1.0*x1, 1.0*x2] in NormCone(1.5, 3)
    """
    _test_atom(target) do context
        return norm(Variable(2), 1.5)
    end
    @test_throws(
        ErrorException("vector p-norms not defined for p < 1"),
        norm(Variable(2), 0.5),
    )
    return
end

### reformulations/partialtranspose

function test_partialtranspose()
    target = """
    variables: x1, x2, x3, x4
    minobjective: [1.0*x1, 1.0*x2, 1.0*x3, 1.0*x4]
    [1.0+x1, 2.0+x2, 3.0+x3, 4.0+x4] in Nonnegatives(4)
    """
    _test_reformulation(target) do context
        x = Variable(2, 2)
        add_constraint!(context, [1 3; 2 4] + x >= 0)
        return partialtranspose(x, 1, [1, 2])
    end
    target = """
    variables: x1, x2, x3, x4
    minobjective: [1.0*x1, 1.0*x3, 1.0*x2, 1.0*x4]
    [1.0+x1, 2.0+x2, 3.0+x3, 4.0+x4] in Nonnegatives(4)
    """
    _test_reformulation(target) do context
        x = Variable(2, 2)
        add_constraint!(context, [1 3; 2 4] + x >= 0)
        return partialtranspose(x, 1, [2, 1])
    end
    x = [1 2 3; 4 5 6; 7 8 9]
    @test partialtranspose(x, 1, [1, 3]) == x
    @test partialtranspose(x, 1, [3, 1]) == x'
    @test partialtranspose(x, 2, [3, 1]) == x
    @test partialtranspose(x, 2, [1, 3]) == x'
    for x in (Variable(2, 3), [1 2; 3 4; 5 6])
        @test_throws(
            ArgumentError("Only square matrices are supported"),
            partialtranspose(x, 1, [2, 3]),
        )
    end
    for x in (Variable(2, 2), [1 2; 3 4])
        @test_throws(
            ArgumentError("Invalid system, should between 1 and 2; got 0"),
            partialtranspose(x, 0, [2, 2]),
        )
        @test_throws(
            ArgumentError("Invalid system, should between 1 and 2; got 3"),
            partialtranspose(x, 3, [2, 2]),
        )
        @test_throws(
            ArgumentError(
                "Dimension of system doesn't correspond to dimension of subsystems",
            ),
            partialtranspose(x, 1, [2, 2]),
        )
    end
    return
end

### reformulations/quadform

function test_quadform()
    target = """
    variables: x11, x21, x12, x22
    minobjective: 1.0 * x11 + 2.0 * x21 + 2.0 * x12 + 4.0 * x22
    """
    _test_reformulation(target) do context
        return quadform([1, 2], Variable(2, 2))
    end
    _test_reformulation(target) do context
        return quadform([1, 2], Variable(2, 2); assume_psd = true)
    end
    _test_reformulation(target) do context
        return quadform(constant([1, 2]), Variable(2, 2))
    end
    target = """
    variables: u, t, x1, x2
    minobjective: 1.0 * u
    [t, 3.999999999999999*x1+1.9999999999999998*x2, 1.9999999999999998*x1+3.999999999999999*x2] in SecondOrderCone(3)
    [u, 0.5, t] in RotatedSecondOrderCone(3)
    """
    _test_reformulation(target) do context
        return quadform(Variable(2), [20.0 16.0; 16.0 20.0])
    end
    _test_reformulation(target) do context
        return quadform(Variable(2), [20.0 16.0; 16.0 20.0]; assume_psd = true)
    end
    _test_reformulation(target) do context
        return quadform(Variable(2), constant([20.0 16.0; 16.0 20.0]))
    end
    target = """
    variables: u, t, x1, x2
    minobjective: 1.0 + -1.0 * u
    [t, 3.999999999999999*x1+1.9999999999999998*x2, 1.9999999999999998*x1+3.999999999999999*x2] in SecondOrderCone(3)
    [u, 0.5, t] in RotatedSecondOrderCone(3)
    """
    _test_reformulation(target) do context
        return 1 + quadform(Variable(2), -[20 16; 16 20])
    end
    _test_reformulation(target) do context
        return 1 + quadform(Variable(2), constant(-[20 16; 16 20]))
    end
    H = Variable(2, 2)
    fix!(H, [1 0; 0 1])
    @test_throws(
        ErrorException(
            "Convex.jl v0.13.5 introduced the ability to use `fix!`ed variables " *
            "in `quadform`. However, this did not consider the case that the " *
            "value of `fix!`ed variables is changed between solves. Due to the " *
            "risk that this may silently produce incorrect solutions, this " *
            "behavior has been removed. Use `evaluate(H)` to obtain the value of " *
            "a fixed variable. If the value changes between solves, rebuild the " *
            "problem for the change to take effect.",
        ),
        quadform(Variable(2), H),
    )
    @test_throws(
        ErrorException("quadform only takes square matrices"),
        quadform(Variable(2), [1 2; 3 4; 5 6]),
    )
    @test_throws(
        ErrorException("quadform only defined for Hermitian matrices"),
        quadform(Variable(2), [1 0; -2 1]),
    )
    return
end

### reformulations/tr

function test_tr()
    target = """
    variables: x1, x2, x3, x4
    minobjective: 1.0 * x1 + 1.0 * x4
    """
    _test_reformulation(target) do context
        return tr(Variable(2, 2))
    end
    return
end

### reformulations/transpose

function test_transpose()
    target = """
    variables: x1, x2, x3, x4
    minobjective: [1.0 * x1, 1.0 * x3, 1.0 * x2, 1.0 * x4]
    [1.0 + x1, 3.0 + x2, 2.0 + x3, 4.0 + x4] in Nonnegatives(4)
    """
    _test_reformulation(target) do context
        x = Variable(2, 2)
        add_constraint!(context, x + [1 2; 3 4] >= 0)
        return x'
    end
    _test_reformulation(target) do context
        x = Variable(2, 2)
        add_constraint!(context, x + [1 2; 3 4] >= 0)
        return LinearAlgebra.transpose(x)
    end
    _test_reformulation(target) do context
        x = Variable(2, 2)
        add_constraint!(context, x + [1 2; 3 4] >= 0)
        return LinearAlgebra.adjoint(x)
    end
    for a in (2, 2.0, [1, 2], [1 2; 3 4], 2 + 3im, 2 - 3im)
        b = constant(a)
        @test evaluate(LinearAlgebra.transpose(b)) ==
              LinearAlgebra.transpose(evaluate(b))
        @test evaluate(LinearAlgebra.adjoint(b)) ==
              LinearAlgebra.adjoint(evaluate(b))
    end
    return
end

end  # TestAtoms

TestAtoms.runtests()
