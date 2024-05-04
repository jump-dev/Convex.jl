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
    # ```julia
    # for i in 1:m, j in 1:n
    #     f = vcat(x[i, j], 1, z[i, j])
    #     add_constraint!(context, GenericConstraint{MOI.ExponentialCone}(f))
    # end
    # return conic_form!(context, z)
    # ```
    # This is slow because we are indexing on the Convex side, and Convex is
    # based around vector/matrix operations. We don't want to produce n*m
    # IndexAtoms!
    #
    # Instead, we will drop to the MOI level to implement this in terms of
    # scalar operations.
    x_tape = conic_form!(context, x)
    # Since `ExpAtom` is restricted to `sign(x)` being real, `x_tape` is either
    # `SPARSE_VECTOR{T}` or `SparseTape{T}`. The `SPARSE_VECTOR{T}` case happens
    # when `x` is a constant (or a function of a constant). In this case, we can
    # just take `exp` directly.
    if x_tape isa SPARSE_VECTOR
        return exp.(x_tape)
    end
    # Otherwise, convert x_tape into a vector of `MOI.ScalarAffineFunction`s.
    xs = MOI.Utilities.scalarize(to_vaf(x_tape))
    z_tape = conic_form!(context, z)
    for (xi, zi) in zip(xs, z_tape.variables)
        MOI.add_constraint(
            context.model,
            MOI.Utilities.operate(vcat, T, xi, T(1), zi),
            MOI.ExponentialCone(),
        )
    end
    return z_tape
end
