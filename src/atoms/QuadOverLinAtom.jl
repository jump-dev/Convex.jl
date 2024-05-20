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
        if iscomplex(y)
            error(
                "[QuadOverLinAtom] the second argument to quadoverlin must be real, not complex.",
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
    # `real` is only necessary to fix the type; `x'*x` will always be real-valued.
    return real(output(x' * x)) / evaluate(q.children[2])
end

function new_conic_form!(context::Context{T}, q::QuadOverLinAtom) where {T}
    t = Variable()
    x, y = q.children
    if iscomplex(x)
        # ||x||₂² = ∑ᵢ |xᵢ|^2 = ∑ᵢ [re(xᵢ)² + im(xᵢ)²]
        #         = ||re(x)||₂² + ||im(x)||₂²
        #         = || vcat(re(x), im(y)) ||₂²
        f = vcat(t, (1 / T(2)) * y, vcat(real(x), imag(x)))
    else
        f = vcat(t, (1 / T(2)) * y, x)
    end
    add_constraint!(context, Constraint{MOI.RotatedSecondOrderCone}(f))
    return conic_form!(context, t)
end

quadoverlin(x::AbstractExpr, y::AbstractExpr) = QuadOverLinAtom(x, y)
