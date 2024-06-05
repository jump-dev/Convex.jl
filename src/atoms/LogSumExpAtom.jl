# Copyright (c) 2014: Madeleine Udell and contributors
#
# Use of this source code is governed by a BSD-style license that can be found
# in the LICENSE file or at https://opensource.org/license/bsd-2-clause

"""
    LogSumExpAtom(x::AbstractExpr, dims::Union{Colon,Int} = :)

Represents the expression `log.(sum(exp.(x); dims))`.
"""
mutable struct LogSumExpAtom <: AbstractExpr
    children::Tuple{AbstractExpr}
    size::Tuple{Int,Int}
    dims::Union{Colon,Int}

    function LogSumExpAtom(x::AbstractExpr, dims::Union{Colon,Int} = :)
        @assert dims == Colon() || 1 <= dims <= 2
        if sign(x) == ComplexSign()
            error(
                "[LogSumExpAtom] the argument should be real but it's instead complex",
            )
        end
        m = dims == 2 ? size(x, 1) : 1
        n = dims == 1 ? size(x, 2) : 1
        return new((x,), (m, n), dims)
    end
end

head(io::IO, ::LogSumExpAtom) = print(io, "logsumexp")

Base.sign(::LogSumExpAtom) = NoSign()

monotonicity(::LogSumExpAtom) = (Nondecreasing(),)

curvature(::LogSumExpAtom) = ConvexVexity()

function evaluate(x::LogSumExpAtom)
    _x = evaluate(x.children[1])
    max_x = maximum(_x)
    return max_x .+ log.(sum(exp.(_x .- max_x); x.dims))
end

function new_conic_form!(context::Context, e::LogSumExpAtom)
    # log(sum(exp(x))) <= t  <=>  sum(exp(x)) <= exp(t) <=> sum(exp(x - t)) <= 1
    x = only(e.children)
    t = Variable(size(e))
    y = if e.dims == 1  # t is a row-vector
        ones(size(x, 1), 1) * t
    elseif e.dims == 2  # t is a col-vector
        t * ones(1, size(x, 2))
    else
        t * ones(size(x))
    end
    add_constraint!(context, 1 >= sum(exp(x - y); dims = e.dims))
    return conic_form!(context, t)
end
