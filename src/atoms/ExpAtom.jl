# Copyright (c) 2014: Madeleine Udell and contributors
#
# Use of this source code is governed by a BSD-style license that can be found
# in the LICENSE file or at https://opensource.org/license/bsd-2-clause

mutable struct ExpAtom <: AbstractExpr
    children::Tuple{AbstractExpr}
    size::Tuple{Int,Int}

    function ExpAtom(x::AbstractExpr)
        if sign(x) == ComplexSign()
            error(
                "[ExpAtom] the argument should be real but it's instead complex",
            )
        end
        return new((x,), x.size)
    end
end

head(io::IO, ::ExpAtom) = print(io, "exp")

Base.sign(::ExpAtom) = Positive()

monotonicity(::ExpAtom) = (Nondecreasing(),)

curvature(::ExpAtom) = ConvexVexity()

evaluate(x::ExpAtom) = exp.(evaluate(x.children[1]))

Base.exp(x::AbstractExpr) = ExpAtom(x)

function new_conic_form!(context::Context{T}, e::ExpAtom) where {T}
    # exp(x) \leq t  <=>  (x,1,t) \in ExpCone
    x = e.children[1]
    x_tape = conic_form!(context, x)
    if x_tape isa SPARSE_VECTOR
        return exp.(x_tape)
    end
    return vectorized_exp_cone_triples!(context, x)
end

# constrains `(x[i], 1, t[i]) ∈ ExponentialCone` for each element of `x` and returns `t`. Optionally pass a permutation to reorder the triple before it is constrained.
function vectorized_exp_cone_triples!(
    context::Context{T},
    x,
    permutation = (1, 2, 3),
) where {T}
    @assert issetequal(permutation, (1, 2, 3))
    @assert length(permutation) == 3
    m, n = size(x)
    t = Variable(m, n)
    # Naive implementation:
    # t = Variable(m, n)
    # for i in 1:m, j in 1:n
    #     f = vcat(x[i, j], 1, t[i, j]) # up to permutation
    #     add_constraint!(context, GenericConstraint{MOI.ExponentialCone}(f))
    # end
    # return conic_form!(context, t)
    # This is slow, since we are indexing on the Convex side, and convex is based around
    # vector/matrix operations. We don't want to produce n*m IndexAtoms!
    # Instead, we will drop to the MOI level to implement this in terms of scalar operations.
    t_tape = conic_form!(context, t)
    # We just created `t`, so we know its "operation" is trivial.
    # Therefore we can just take the variables to get a vector of `MOI.VariableIndex`
    ts = t_tape.variables

    x_tape = conic_form!(context, x)
    xs = MOI.Utilities.scalarize(to_vaf(x_tape))
    # Now we have a vector of `m*n` ScalarAffineFunctions in order.

    # now we simply add the constraints
    for (xi, ti) in zip(xs, ts)
        tuple = (xi, T(1), ti)
        ordered_tuple = ntuple(i -> tuple[permutation[i]], Val(3))
        MOI.add_constraint(
            context.model,
            MOI.Utilities.operate(vcat, T, ordered_tuple...),
            MOI.ExponentialCone(),
        )
    end
    return t_tape
end
