# Copyright (c) 2014: Madeleine Udell and contributors
#
# Use of this source code is governed by a BSD-style license that can be found
# in the LICENSE file or at https://opensource.org/license/bsd-2-clause

mutable struct ConjugateAtom <: AbstractExpr
    children::Tuple{AbstractExpr}
    size::Tuple{Int,Int}

    ConjugateAtom(x::AbstractExpr) = new((x,), x.size)
end

head(io::IO, ::ConjugateAtom) = print(io, "conj")

Base.sign(x::ConjugateAtom) = sign(x.children[1])

monotonicity(::ConjugateAtom) = (Nondecreasing(),)

curvature(::ConjugateAtom) = ConstVexity()

evaluate(x::ConjugateAtom) = conj(evaluate(x.children[1]))

function new_conic_form!(context::Context{T}, x::ConjugateAtom) where {T}
    objective = conic_form!(context, only(AbstractTrees.children(x)))
    return operate(conj, T, sign(x), objective)
end
