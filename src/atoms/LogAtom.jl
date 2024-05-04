# Copyright (c) 2014: Madeleine Udell and contributors
#
# Use of this source code is governed by a BSD-style license that can be found
# in the LICENSE file or at https://opensource.org/license/bsd-2-clause

mutable struct LogAtom <: AbstractExpr
    children::Tuple{AbstractExpr}
    size::Tuple{Int,Int}

    function LogAtom(x::AbstractExpr)
        if sign(x) == ComplexSign()
            error(
                "[LogAtom] the argument should be real but it's instead complex",
            )
        end
        return new((x,), x.size)
    end
end

head(io::IO, ::LogAtom) = print(io, "log")

Base.sign(::LogAtom) = NoSign()

monotonicity(::LogAtom) = (Nondecreasing(),)

curvature(::LogAtom) = ConcaveVexity()

evaluate(x::LogAtom) = log.(evaluate(x.children[1]))

Base.log(x::AbstractExpr) = LogAtom(x)

function new_conic_form!(context::Context, e::LogAtom)
    # log(x) \geq t  <=> (t,1,x) \in ExpCone
    x = e.children[1]

    # to choose the permutation, we want the elements of the constraint to be
    # (t, 1, x)
    # but with the identity permutation, the default is
    # (x, 1, t)
    # So (3, 2, 1) permutes it to the correct order.
    return vectorized_exp_cone_triples!(context, x, (3, 2, 1))
end
