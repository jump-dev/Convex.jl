# Copyright (c) 2014: Madeleine Udell and contributors
#
# Use of this source code is governed by a BSD-style license that can be found
# in the LICENSE file or at https://opensource.org/license/bsd-2-clause

mutable struct DotMultiplyAtom <: AbstractExpr
    children::Tuple{AbstractExpr,AbstractExpr}
    size::Tuple{Int,Int}

    function DotMultiplyAtom(x::AbstractExpr, y::AbstractExpr)
        sz = if x.size[1] == y.size[1] && (x.size[2] == 1 || y.size[2] == 1)
            (x.size[1], max(x.size[2], y.size[2]))
        elseif x.size[2] == y.size[2] && (x.size[1] == 1 || y.size[1] == 1)
            (max(x.size[1], y.size[1]), y.size[2])
        elseif x.size == y.size
            x.size
        else
            error(
                "[DotMultiplyAtom] cannot multiply two expressions of sizes $(x.size) and $(y.size)",
            )
        end
        return new((x, y), sz)
    end
end

head(io::IO, ::DotMultiplyAtom) = print(io, ".*")

Base.sign(x::DotMultiplyAtom) = sign(x.children[1]) * sign(x.children[2])

function monotonicity(x::DotMultiplyAtom)
    return (
        sign(x.children[2]) * Nondecreasing(),
        sign(x.children[1]) * Nondecreasing(),
    )
end

function curvature(x::DotMultiplyAtom)
    lhs, rhs = x.children
    if vexity(lhs) != ConstVexity() && vexity(rhs) != ConstVexity()
        return NotDcp()
    end
    return ConstVexity()
end

function evaluate(x::DotMultiplyAtom)
    return reshape(evaluate(x.children[1]) .* evaluate(x.children[2]), size(x))
end

function new_conic_form!(context::Context{T}, x::DotMultiplyAtom) where {T}
    lhs, rhs = x.children
    if vexity(lhs) != ConstVexity()
        if vexity(rhs) != ConstVexity()
            error(
                "[DotMultiplyAtom] multiplication of two non-constant expressions is not DCP compliant",
            )
        end
        # Switch arguments so that the left-hand side is constant
        lhs, rhs = rhs, lhs
    end
    # Start by assuming that the constant lhs matrix is the smaller object that
    # will be broadcast over the larger RHS object. Let Julia automatically
    # resize it by .* by `ones`.
    coef = evaluate(lhs) .* ones(T, size(rhs))
    if size(coef) != size(rhs)
        # If coef is not the same size as rhs, then we must be broadcasting the
        # smaller rhs object over the larger coef. In this case, rhs must be a
        # row or column vector.
        if size(rhs, 1) == 1
            # rhs is a row vector. Stretch it out to have the same number of
            # rows as coef.
            rhs = ones(T, size(coef, 1)) * rhs
        else
            @assert size(rhs, 2) == 1
            # rhs is a col vector. Stretch it out to have the same number of
            # columns as coef.
            rhs = rhs * ones(T, 1, size(coef, 2))
        end
    end
    # For sanity, check that these are the same size.
    @assert size(coef) == size(rhs)
    # Represent the array x .* y as D(x) * y
    ret = SparseArrays.sparse(LinearAlgebra.Diagonal(vec(coef))) * vec(rhs)
    return conic_form!(context, reshape(ret, size(rhs, 1), size(rhs, 2)))
end

function Base.Broadcast.broadcasted(
    ::typeof(*),
    x::AbstractExpr,
    y::AbstractExpr,
)
    if isequal(x, y)
        return square(x)
    elseif x.size == (1, 1) || y.size == (1, 1)
        return x * y
    end
    return DotMultiplyAtom(x, y)
end

function Base.Broadcast.broadcasted(::typeof(*), x::Value, y::AbstractExpr)
    return constant(x) .* y
end

function Base.Broadcast.broadcasted(::typeof(*), x::AbstractExpr, y::Value)
    return x .* constant(y)
end

function Base.Broadcast.broadcasted(::typeof(/), x::AbstractExpr, y::Value)
    return x .* constant(1 ./ y)
end
