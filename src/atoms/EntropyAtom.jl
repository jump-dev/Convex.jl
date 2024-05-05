# Copyright (c) 2014: Madeleine Udell and contributors
#
# Use of this source code is governed by a BSD-style license that can be found
# in the LICENSE file or at https://opensource.org/license/bsd-2-clause

mutable struct EntropyAtom <: AbstractExpr
    children::Tuple{AbstractExpr}
    size::Tuple{Int,Int}

    function EntropyAtom(x::AbstractExpr)
        if sign(x) == ComplexSign()
            error(
                "[EntropyAtom] the argument should be real but it's instead complex",
            )
        end
        # TODO check positivity or enforce it
        return new((x,), size(x))
    end
end

head(io::IO, ::EntropyAtom) = print(io, "entropy")

Base.sign(::EntropyAtom) = NoSign()

monotonicity(::EntropyAtom) = (NoMonotonicity(),)

curvature(::EntropyAtom) = ConcaveVexity()

function evaluate(x::EntropyAtom)
    c = evaluate(x.children[1])
    return -c .* log.(c)
end

entropy(x::AbstractExpr) = sum(EntropyAtom(x))

entropy_elementwise(x::AbstractExpr) = EntropyAtom(x)

function new_conic_form!(context::Context, e::EntropyAtom)
    # -x log(x) >= t  <=>  x exp(t/x) <= 1 <=> (t, x, 1) in ExponentialCone()
    x = e.children[1]
    # To choose the permutation, we want the elements of the constraint to be
    #   (t, x, 1) in ExponentialCone()
    # but with the identity permutation, the default is:
    #   (x, 1, t) in ExponentialCone()
    # so [3, 1, 2] permutes it to the correct order.
    return _add_vectorized_exp_cone(context, x, [3, 1, 2])
end
