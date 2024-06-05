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

function new_conic_form!(context::Context, e::LogAtom)
    # log(x) >= t  <=> (t, 1, x) in ExponentialCone()
    x = e.children[1]
    # To choose the permutation, we want the elements of the constraint to be:
    #   (t, 1, x) in ExponentialCone()
    # The default is:
    #   (x, 1, t) in ExponentialCone()
    # so [3, 2, 1] permutes it to the correct order.
    return _add_vectorized_exp_cone(context, x, [3, 2, 1])
end
