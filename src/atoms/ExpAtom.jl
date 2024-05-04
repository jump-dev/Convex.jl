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
    # exp(x) \leq z  <=>  (x,1,z) \in ExpCone
    x = e.children[1]
    m, n = size(x)
    z = Variable(m, n)
    # Naive implementation:
    # for i in 1:m, j in 1:n
    # f = vcat(x[i, j], 1, z[i, j])
        # add_constraint!(context, GenericConstraint{MOI.ExponentialCone}(f))
    # end
    # return conic_form!(context, z)
    # This is slow, since we are indexing on the Convex side, and convex is based around
    # vector/matrix operations. We don't want to produce n*m IndexAtoms!
    # Instead, we will drop to the MOI level to implement this in terms of scalar operations.
    # First, we will get `x` as an MOI.VectorAffineFunction
    x_tape = conic_form!(context, x)
    vaf = to_vaf(x_tape)
    # Next, we can extract the individual components of `x` via `MOI.Utilities.scalarize`
    xs = MOI.Utilities.scalarize(vaf)
    # Now we have a vector of `m*n` ScalarAffineFunctions in order.
    # We can likewise lower `z` to a vector of `MOI.VariableIndex`
    z_tape = conic_form!(context, z)
    zs = z_tape.variables
    for i in eachindex(xs, zs)
        # Now, we wish to add the constraint `(x[i], 1, z[i]) ∈ MOI.ExponentialCone()`
        # however, we have 3 different types: x[i] is a ScalarAffineFunction, 1 is a constant,
        # and `z` is a VariableIndex.
        # So we can't use `MOI.Utilities.vectorize`. Instead, we will construct a VectorAffineFunction manually.
        # First, we construct the VectorAffineTerm's for the first and third components.
        terms = [
            [MOI.VectorAffineTerm(Int64(1), sat) for sat in xs[i].terms];
            MOI.VectorAffineTerm(Int64(3), MOI.ScalarAffineTerm(T(1), zs[i]))
        ]
        # Then we can add in the constants, and we are good to go.
        vaf_i = MOI.VectorAffineFunction(terms, [xs[i].constant, T(1), T(0)])
        MOI.add_constraint(context.model, vaf_i, MOI.ExponentialCone())
    end
    return z_tape
end
