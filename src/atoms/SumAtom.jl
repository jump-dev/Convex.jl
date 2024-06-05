# Copyright (c) 2014: Madeleine Udell and contributors
#
# Use of this source code is governed by a BSD-style license that can be found
# in the LICENSE file or at https://opensource.org/license/bsd-2-clause

mutable struct SumAtom <: AbstractExpr
    children::Tuple{AbstractExpr}
    size::Tuple{Int,Int}

    SumAtom(x::AbstractExpr) = new((x,), (1, 1))
end

head(io::IO, ::SumAtom) = print(io, "sum")

Base.sign(x::SumAtom) = sign(x.children[1])

monotonicity(::SumAtom) = (Nondecreasing(),)

curvature(::SumAtom) = ConstVexity()

evaluate(x::SumAtom) = sum(evaluate(x.children[1]))

function new_conic_form!(context::Context{T}, A::SumAtom) where {T}
    subobj = conic_form!(context, only(AbstractTrees.children(A)))
    return operate(sum, T, sign(A), subobj)
end
