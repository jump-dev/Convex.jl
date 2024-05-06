# Copyright (c) 2014: Madeleine Udell and contributors
#
# Use of this source code is governed by a BSD-style license that can be found
# in the LICENSE file or at https://opensource.org/license/bsd-2-clause

mutable struct MultiplyAtom <: AbstractExpr
    children::Tuple{AbstractExpr,AbstractExpr}
    size::Tuple{Int,Int}

    function MultiplyAtom(x::AbstractExpr, y::AbstractExpr)
        sz = if x.size == (1, 1)
            y.size
        elseif y.size == (1, 1)
            x.size
        elseif x.size[2] == y.size[1]
            (x.size[1], y.size[2])
        else
            error(
                "[MultiplyAtom] cannot multiply two expressions of sizes $(x.size) and $(y.size)",
            )
        end
        return new((x, y), sz)
    end
end

head(io::IO, ::MultiplyAtom) = print(io, "*")

Base.sign(x::MultiplyAtom) = sign(x.children[1]) * sign(x.children[2])

function monotonicity(x::MultiplyAtom)
    return (
        sign(x.children[2]) * Nondecreasing(),
        sign(x.children[1]) * Nondecreasing(),
    )
end

# Multiplication has an indefinite hessian, so if neither children are constants,
# the curvature of the atom will violate DCP.
function curvature(x::MultiplyAtom)
    if vexity(x.children[1]) != ConstVexity() &&
       vexity(x.children[2]) != ConstVexity()
        return NotDcp()
    end
    return ConstVexity()
end

evaluate(x::MultiplyAtom) = evaluate(x.children[1]) * evaluate(x.children[2])

function _complex_convert(::Type{T}, x) where {T}
    if iscomplex(x)
        return real_convert(Complex{T}, x)
    end
    return real_convert(T, x)
end

real_convert(::Type{T}, x::Number) where {T} = T(x)

real_convert(::Type{T}, x::AbstractMatrix) where {T} = create_sparse(T, x)

real_convert(::Type{T}, x::SPARSE_MATRIX{T}) where {T} = x

real_convert(::Type{T}, x::SPARSE_VECTOR{T}) where {T} = x

real_convert(::Type{T}, x::AbstractVector) where {T} = SPARSE_VECTOR{T}(x)

function new_conic_form!(context::Context{T}, x::MultiplyAtom) where {T}
    if curvature(x) == NotDcp()
        error(
            "[MultiplyAtom] multiplication of two non-constant expressions is not DCP compliant",
        )
    end
    if x.children[1].size == (1, 1) || x.children[2].size == (1, 1)
        # scalar multiplication
        if vexity(x.children[1]) == ConstVexity()
            lhs, expr = x.children
        else
            @assert vexity(x.children[2]) == ConstVexity()
            expr, lhs = x.children
        end
        lhs = if lhs.size == (1, 1)
            # Make sure all 1x1 sized objects are interpreted as scalars, since
            # [1] * [1, 2, 3] is illegal in julia, but 1 * [1, 2, 3] is ok
            only(evaluate(lhs))
        else
            reshape(evaluate(lhs), length(lhs), 1)
        end
        lhs = _complex_convert(T, lhs)
        rhs = conic_form!(context, expr)
        return operate(add_operation, T, sign(x), lhs, rhs)
    elseif vexity(x.children[1]) == ConstVexity()
        # left matrix multiplication
        lhs = _complex_convert(T, evaluate(x.children[1]))
        return operate(
            add_operation,
            T,
            sign(x),
            kron(spidentity(T, x.size[2]), lhs),
            conic_form!(context, x.children[2]),
        )
    else
        @assert vexity(x.children[2]) == ConstVexity()
        # right matrix multiplication
        rhs = _complex_convert(T, evaluate(x.children[2]))
        return operate(
            add_operation,
            T,
            sign(x),
            kron(transpose(rhs), spidentity(T, x.size[1])),
            conic_form!(context, x.children[1]),
        )
    end
end

function Base.:*(x::AbstractExpr, y::AbstractExpr)
    if isequal(x, y) && x.size == (1, 1)
        return square(x)
    end
    return MultiplyAtom(x, y)
end

Base.:*(x::Value, y::AbstractExpr) = MultiplyAtom(constant(x), y)

Base.:*(x::AbstractExpr, y::Value) = MultiplyAtom(x, constant(y))

Base.:/(x::AbstractExpr, y::Value) = MultiplyAtom(x, constant(1 ./ y))

function _dot_multiply(x, y)
    if size(x) == (1, 1) || size(y) == (1, 1)
        return x * y
    end
    if vexity(x) != ConstVexity()
        if vexity(y) != ConstVexity()
            error(
                "[MultiplyAtom] multiplication of two non-constant expressions is not DCP compliant",
            )
        end
        x, y = y, x
    end
    # promote the size of the coefficient matrix, so e.g., 3 .* x works
    # regardless of the size of x
    coeff = evaluate(x) .* ones(size(y))
    # Promote the size of the variable. We've previously ensured neither x nor y
    # is 1x1 and that the sizes are compatible, so if the sizes aren't equal the
    # smaller one is size 1.
    if size(y, 1) < size(coeff, 1)
        y = ones(size(coeff, 1)) * y
    elseif size(y, 2) < size(coeff, 2)
        y = y * ones(1, size(coeff, 1))
    end
    ret = SparseArrays.sparse(LinearAlgebra.Diagonal(vec(coeff))) * vec(y)
    return reshape(ret, size(y, 1), size(y, 2))
end

# if neither is a constant it's not DCP, but might be nice to support anyway for
# eg MultiConvex
function Base.Broadcast.broadcasted(
    ::typeof(*),
    x::AbstractExpr,
    y::AbstractExpr,
)
    if isequal(x, y)
        return square(x)
    end
    return _dot_multiply(x, y)
end

function Base.Broadcast.broadcasted(::typeof(*), x::Value, y::AbstractExpr)
    return _dot_multiply(constant(x), y)
end

function Base.Broadcast.broadcasted(::typeof(*), x::AbstractExpr, y::Value)
    return _dot_multiply(constant(y), x)
end

function Base.Broadcast.broadcasted(::typeof(/), x::AbstractExpr, y::Value)
    return _dot_multiply(constant(1 ./ y), x)
end

# x ./ y and x / y for x constant, y variable is defined in
# second_order_cone/qol_elemwise.jl
