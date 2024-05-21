# Copyright (c) 2014: Madeleine Udell and contributors
#
# Use of this source code is governed by a BSD-style license that can be found
# in the LICENSE file or at https://opensource.org/license/bsd-2-clause

# Returns:
# [logsumexp(v) for v in eachcol(x)]
# where `logsumexp(v) = log(sum(exp, v))` except hopefully more numerically stable
mutable struct ColwiseLogSumExpAtom <: AbstractExpr
    children::Tuple{AbstractExpr}
    size::Tuple{Int,Int}

    function ColwiseLogSumExpAtom(x::AbstractExpr)
        if sign(x) == ComplexSign()
            error(
                "[ColwiseLogSumExpAtom] the argument should be real but it's instead complex",
            )
        end
        ncols = size(x, 2)
        return new((x,), (1, ncols))
    end
end

head(io::IO, ::ColwiseLogSumExpAtom) = print(io, "logsumexp")

Base.sign(::ColwiseLogSumExpAtom) = NoSign()

monotonicity(::ColwiseLogSumExpAtom) = (Nondecreasing(),)

curvature(::ColwiseLogSumExpAtom) = ConvexVexity()

function evaluate(x::ColwiseLogSumExpAtom)
    _x = evaluate(x.children[1])
    return map(eachcol(_x)) do col
        max_x = maximum(col)
        return max_x + log(sum(exp.(col .- max_x)))
    end
end

# We call `vec` first, to treat this as one long vector to `logsumexp`
# (otherwise, `ColwiseLogSumExpAtom` computes `logsumexp` over each column)
logsumexp(x::AbstractExpr) = ColwiseLogSumExpAtom(vec(x))[1]

function new_conic_form!(context::Context, e::ColwiseLogSumExpAtom)
    x = e.children[1]
    ncols = size(x, 2)
    # log(sum(exp(x))) <= t  <=>  sum(exp(x)) <= exp(t) <=> sum(exp(x - t)) <= 1
    t = Variable(ncols)
    z = sum(exp(x - transpose(t) .* ones(size(x))); dims = 1)
    add_constraint!(context, 1 >= z)
    return conic_form!(context, t)
end

function logisticloss(e::AbstractExpr)
    if length(e) == 1
        return logsumexp([e; 0])
    end
    n = length(e)
    e = transpose(vec(e))
    return sum(ColwiseLogSumExpAtom(vcat(e, zeros(1, n))))
end
