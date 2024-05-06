# Copyright (c) 2014: Madeleine Udell and contributors
# Copyright (c) 2021: Hamza Fawzi
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

"""
    GeomMeanEpiConeSquare(t::Rational, side_dimension::Int)

The constraint `(T, A, B) in GeomMeanEpiConeSquare(t, side_dimension)`
constrains T to

  A #_t B ⪯ T

where:

  * A #_t B is the t-weighted geometric mean of A and B:
    A^{1/2} (A^{-1/2} B A^{-1/2})^t A^{1/2}
  * Parameter t must be in [-1, 0] or [1, 2].
  * Constraints A ⪰ 0, B ⪰ 0 are added.

Note on fullhyp: GeomMeanEpiConeSquare will always return a full epigraph cone
(unlike GeomMeanHypoCone) and so this parameter is not really used. It is here
just for consistency with the GeomMeanHypoCone function.

## Reference

Ported from CVXQUAD which is based on the paper: "Lieb's concavity theorem,
matrix geometric means and semidefinite optimization" by Hamza Fawzi and James
Saunderson (arXiv:1512.03401)
"""
struct GeomMeanEpiConeSquare <: MOI.AbstractVectorSet
    t::Rational
    side_dimension::Int

    function GeomMeanEpiConeSquare(t::Rational, side_dimension::Int)
        if !(-1 <= t <= 0 || 1 <= t <= 2)
            throw(DomainError(t, "t must be in the range [-1, 0] or [1, 2]"))
        end
        return new(t, side_dimension)
    end
end

MOI.dimension(set::GeomMeanEpiConeSquare) = 3 * set.side_dimension

head(io::IO, ::GeomMeanEpiConeSquare) = print(io, "∈(GeomMeanEpiConeSquare)")

function Base.in(func::Tuple, set::GeomMeanEpiConeSquare)
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
    return GenericConstraint(vcat(vec.(func)...), set)
end

# For t ∈ [-1, 0] ∪ [1, 2] the t-weighted matrix geometric mean is matrix convex
# (arxiv:1512.03401).
# So if A and B are convex sets, then T ⪰ A #_t B will be a convex set.
function vexity(constraint::GenericConstraint{GeomMeanEpiConeSquare})
    n = constraint.set.side_dimension
    d = n^2
    I = 1:d
    T = vexity(constraint.child[I])
    A = vexity(constraint.child[d.+I])
    B = vexity(constraint.child[2d.+I])
    if A in (ConcaveVexity(), NotDcp()) || B in (ConcaveVexity(), NotDcp())
        return NotDcp()
    end
    return ConvexVexity() + (-T)
end

function _add_constraint!(
    context::Context,
    constraint::GenericConstraint{GeomMeanEpiConeSquare},
)
    n = constraint.set.side_dimension
    d = n^2
    I = 1:d
    T = reshape(constraint.child[I], n, n)
    A = reshape(constraint.child[d.+I], n, n)
    B = reshape(constraint.child[2d.+I], n, n)
    t = constraint.set.t
    is_complex =
        sign(A) == ComplexSign() ||
        sign(B) == ComplexSign() ||
        sign(T) == ComplexSign()
    Z = if is_complex
        HermitianSemidefinite(size(A)[1])
    else
        Semidefinite(size(A)[1])
    end
    if t <= 0
        add_constraint!(context, [T A; A Z] ⪰ 0)
        add_constraint!(context, Z in GeomMeanHypoCone(A, B, -t, false))
    else
        @assert t >= 1 # range of t checked in GeomMeanEpiConeSquare constructor
        add_constraint!(context, [T B; B Z] ⪰ 0)
        add_constraint!(context, Z in GeomMeanHypoCone(A, B, 2 - t, false))
    end
    return
end
