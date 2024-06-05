# Copyright (c) 2014: Madeleine Udell and contributors
#
# Use of this source code is governed by a BSD-style license that can be found
# in the LICENSE file or at https://opensource.org/license/bsd-2-clause

mutable struct NegateAtom <: AbstractExpr
    children::Tuple{AbstractExpr}
    size::Tuple{Int,Int}

    NegateAtom(x::AbstractExpr) = new((x,), x.size)
end

Base.sign(x::NegateAtom) = -sign(x.children[1])

monotonicity(::NegateAtom) = (Nonincreasing(),)

curvature(::NegateAtom) = ConstVexity()

evaluate(x::NegateAtom) = -evaluate(x.children[1])

function new_conic_form!(context::Context{T}, A::NegateAtom) where {T}
    subobj = conic_form!(context, only(AbstractTrees.children(A)))
    if subobj isa Value
        return -subobj
    end
    return operate(-, T, sign(A), subobj)
end
