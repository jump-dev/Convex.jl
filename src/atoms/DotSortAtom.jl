# Copyright (c) 2014: Madeleine Udell and contributors
#
# Use of this source code is governed by a BSD-style license that can be found
# in the LICENSE file or at https://opensource.org/license/bsd-2-clause

mutable struct DotSortAtom <: AbstractExpr
    children::Tuple{AbstractExpr}
    size::Tuple{Int,Int}
    w::Value

    function DotSortAtom(x::AbstractExpr, w::Value)
        if sign(x) == ComplexSign()
            error(
                "[DotSortAtom] argument should be real instead it is $(sign(x))",
            )
        elseif length(w) != length(x)
            error("[DotSortAtom] x and w must be the same size")
        end
        return new((x,), (1, 1), reshape(w, length(x)))
    end
end

head(io::IO, ::DotSortAtom) = print(io, "dotsort")

function Base.sign(x::DotSortAtom)
    if all(>=(0), x.w)
        return sign(x.children[1])
    elseif all(<=(0), x.w)
        return -sign(x.children[1])
    end
    return NoSign()
end

function monotonicity(x::DotSortAtom)
    if all(>=(0), x.w)
        return (Nondecreasing(),)
    elseif all(<=(0), x.w)
        return (Nonincreasing(),)
    end
    return (NoMonotonicity(),)
end

curvature(::DotSortAtom) = ConvexVexity()

function evaluate(x::DotSortAtom)
    y = only(x.children)
    return sum(a * b for (a, b) in zip(sort(vec(evaluate(y))), sort(x.w)))
end

function new_conic_form!(context::Context{T}, x::DotSortAtom) where {T}
    y = only(x.children)
    if size(y, 1) > 1 && size(y, 2) > 1
        y = vec(y)
    end
    μ, ν, e = Variable(size(y)), Variable(size(y)), ones(T, size(y))
    add_constraint!(context, e * ν' + μ * e' >= y * x.w')
    return conic_form!(context, sum(μ) + sum(ν))
end
