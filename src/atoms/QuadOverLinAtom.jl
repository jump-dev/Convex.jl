# Copyright (c) 2014: Madeleine Udell and contributors
#
# Use of this source code is governed by a BSD-style license that can be found
# in the LICENSE file or at https://opensource.org/license/bsd-2-clause

mutable struct QuadOverLinAtom <: AbstractExpr
    children::Tuple{AbstractExpr,AbstractExpr}
    size::Tuple{Int,Int}

    function QuadOverLinAtom(x::AbstractExpr, y::AbstractExpr)
        if x.size[2] != 1 || y.size != (1, 1)
            error(
                "[QuadOverLinAtom] quadoverlin arguments must be a vector and a scalar",
            )
        end
        return new((x, y), (1, 1))
    end
end

head(io::IO, ::QuadOverLinAtom) = print(io, "qol")

Base.sign(::QuadOverLinAtom) = Positive()

function monotonicity(q::QuadOverLinAtom)
    return (sign(q.children[1]) * Nondecreasing(), Nonincreasing())
end

curvature(::QuadOverLinAtom) = ConvexVexity()

function evaluate(q::QuadOverLinAtom)
    x = evaluate(q.children[1])
    return x' * x / evaluate(q.children[2])
end

function new_conic_form!(context::Context{T}, q::QuadOverLinAtom) where {T}
    t = Variable()
    x, y = q.children
    f = vcat(t, (1 / T(2)) * y, x)
    add_constraint!(context, Constraint{MOI.RotatedSecondOrderCone}(f))
    return conic_form!(context, t)
end

quadoverlin(x::AbstractExpr, y::AbstractExpr) = QuadOverLinAtom(x, y)

function sumsquares(x::AbstractExpr)
    if sign(x) == ComplexSign()
        return sumsquares(vcat(real(x), imag(x)))
    elseif size(x, 2) != 1
        return QuadOverLinAtom(reshape(x, length(x), 1), constant(1))
    end
    return QuadOverLinAtom(x, constant(1))
end
