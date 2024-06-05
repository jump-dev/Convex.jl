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

function new_conic_form!(context::Context{T}, e::ExpAtom) where {T}
    # exp(x) <= t  <=>  (x, 1, t) in ExponentialCone()
    x = e.children[1]
    x_tape = conic_form!(context, x)
    if x_tape isa SPARSE_VECTOR
        return exp.(x_tape)
    end
    return _add_vectorized_exp_cone(context, x, [1, 2, 3])
end

"""
    _add_vectorized_exp_cone(
        context::Context{T},
        x,
        permutation::Vector{Int},
    ) where {T}

Constrains `(x[i], 1, t[i]) ∈ ExponentialCone()` for each element of `x` and
returns `t`.

Permutation is a permuted vector of `[1, 2, 3]` to reorder the triple before it
is constrained. This is helpful for `LogAtom` and `EntropyAtom`.

## Motivation

A naive implementation of this method is:
```julia
    m, n = size(x)
    t = Variable(m, n)
    for i in 1:m, j in 1:n
        f = vcat(x[i, j], 1, t[i, j])[collect(permutation)]
        add_constraint!(context, Constraint{MOI.ExponentialCone}(f))
    end
    return conic_form!(context, t)
end
```
This is slow because we are indexing on the Convex side, and Convex is based
around vector/matrix operations. We don't want to produce n*m IndexAtoms!

Instead, we will drop to the MOI level to implement this in terms of scalar
operations.
"""
function _add_vectorized_exp_cone(
    context::Context{T},
    x,
    permutation::Vector{Int},
) where {T}
    @assert issetequal(permutation, (1, 2, 3))
    @assert length(permutation) == 3
    m, n = size(x)
    t = Variable(m, n)
    t_tape = conic_form!(context, t)
    ts = t_tape.variables::Vector{MOI.VariableIndex}
    x_tape = conic_form!(context, x)
    xs = MOI.Utilities.scalarize(to_vaf(x_tape))
    for (xi, ti) in zip(xs, ts)
        args = (xi, T(1), ti)[permutation]
        MOI.add_constraint(
            context.model,
            MOI.Utilities.operate(vcat, T, args...),
            MOI.ExponentialCone(),
        )
    end
    return t_tape
end
