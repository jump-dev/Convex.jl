# Copyright (c) 2014: Madeleine Udell and contributors
# Copyright (c) 2021: Hamza Fawzi
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

"""
    GeometricMeanEpiConeSquare(t::Rational, side_dimension::Int)

The constraint `(T, A, B) in GeometricMeanEpiConeSquare(t, side_dimension)`
constrains T to

  A #_t B ⪯ T

where:

  * A #_t B is the t-weighted geometric mean of A and B:
    A^{1/2} (A^{-1/2} B A^{-1/2})^t A^{1/2}
  * Parameter t must be in [-1, 0] or [1, 2].
  * Constraints A ⪰ 0, B ⪰ 0 are added.

## Reference

Ported from CVXQUAD which is based on the paper: "Lieb's concavity theorem,
matrix geometric means and semidefinite optimization" by Hamza Fawzi and James
Saunderson (arXiv:1512.03401)
"""
struct GeometricMeanEpiConeSquare <: MOI.AbstractVectorSet
    t::Rational
    side_dimension::Int

    function GeometricMeanEpiConeSquare(t::Rational, side_dimension::Int)
        if !(-1 <= t <= 0 || 1 <= t <= 2)
            throw(DomainError(t, "t must be in the range [-1, 0] or [1, 2]"))
        end
        return new(t, side_dimension)
    end
end

MOI.dimension(set::GeometricMeanEpiConeSquare) = 3 * set.side_dimension^2

function head(io::IO, ::GeometricMeanEpiConeSquare)
    return print(io, "GeometricMeanEpiConeSquare")
end

function Constraint(func::Tuple, set::GeometricMeanEpiConeSquare)
    for f in func
        n = LinearAlgebra.checksquare(f)
        if n != set.side_dimension
            throw(
                DimensionMismatch(
                    "Matrix of side dimension `$n` does not match set of side dimension `$(set.side_dimension)`",
                ),
            )
        end
    end
    return Constraint(vcat(vec.(func)...), set)
end

function _get_matrices(c::Constraint{GeometricMeanEpiConeSquare})
    n = c.set.side_dimension
    d = n^2
    T = reshape(c.child[1:d], n, n)
    A = reshape(c.child[d.+(1:d)], n, n)
    B = reshape(c.child[2d.+(1:d)], n, n)
    return T, A, B
end

# For t ∈ [-1, 0] ∪ [1, 2] the t-weighted matrix geometric mean is matrix convex
# (arxiv:1512.03401).
# So if A and B are convex sets, then T ⪰ A #_t B will be a convex set.
function vexity(constraint::Constraint{GeometricMeanEpiConeSquare})
    T, A, B = _get_matrices(constraint)
    if vexity(A) in (ConcaveVexity(), NotDcp()) ||
       vexity(B) in (ConcaveVexity(), NotDcp())
        return NotDcp()
    end
    return ConvexVexity() + -vexity(T)
end

function _add_constraint!(
    context::Context,
    constraint::Constraint{GeometricMeanEpiConeSquare},
)
    T, A, B = _get_matrices(constraint)
    t = constraint.set.t
    n = size(A, 1)
    Z = if sign(constraint.child) == ComplexSign()
        HermitianSemidefinite(n)
    else
        Semidefinite(n)
    end
    if t <= 0
        add_constraint!(context, [T A; A Z] ⪰ 0)
        add_constraint!(
            context,
            Constraint(
                (Z, A, B),
                GeometricMeanHypoConeSquare(-t, n, false),
            ),
        )
    else
        # range of t checked in GeometricMeanEpiConeSquare constructor.
        @assert t >= 1
        add_constraint!(context, [T B; B Z] ⪰ 0)
        add_constraint!(
            context,
            Constraint(
                (Z, A, B),
                GeometricMeanHypoConeSquare(2 - t, n, false),
            ),
        )
    end
    return
end
