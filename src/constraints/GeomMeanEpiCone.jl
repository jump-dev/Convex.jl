# Copyright (c) 2014: Madeleine Udell and contributors
# Copyright (c) 2021: Hamza Fawzi
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

"""
Constrains T to
  A #_t B ⪯ T
where:
  * A #_t B is the t-weighted geometric mean of A and B:
    A^{1/2} (A^{-1/2} B A^{-1/2})^t A^{1/2}
    Parameter t should be in [-1, 0] or [1, 2].
  * Constraints A ⪰ 0, B ⪰ 0 are added.
All expressions and atoms are subtypes of AbstractExpr.
Please read expressions.jl first.

Note on fullhyp: GeomMeanEpiConeSquare will always return a full epigraph cone
(unlike GeomMeanHypoCone) and so this parameter is not really used.  It is
here just for consistency with the GeomMeanHypoCone function.

REFERENCE
  Ported from CVXQUAD which is based on the paper: "Lieb's concavity
  theorem, matrix geometric means and semidefinite optimization" by Hamza
  Fawzi and James Saunderson (arXiv:1512.03401)
"""
struct GeomMeanEpiConeSquare <: MOI.AbstractVectorSet
    t::Rational
    side_dimension::Int

    function GeomMeanEpiConeSquare(t::Rational, side_dimension::Int)
        if t < -1 || (t > 0 && t < 1) || t > 2
            throw(DomainError(t, "t must be in the range [-1, 0] or [1, 2]"))
        end
        return new(t, side_dimension)
    end
end

head(io::IO, ::GeomMeanEpiConeSquare) = print(io, "∈(GeomMeanEpiConeSquare)")

function Base.in(T, cone::GeomMeanEpiConeSquare)
    return GenericConstraint(T, cone)
end

MOI.dimension(set::GeomMeanEpiConeSquare) = 3set.side_dimension

function _dimension_check(child::Tuple, set::GeomMeanEpiConeSquare)
    for c in child
        n = LinearAlgebra.checksquare(c)
        throw(
            DimensionMismatch(
                "Matrix of side dimension `$n` does not match set of side dimension `$(set.side_dimension)`",
            ),
        )
    end
end

# For t ∈ [-1, 0] ∪ [1, 2] the t-weighted matrix geometric mean is matrix convex
# (arxiv:1512.03401).
# So if A and B are convex sets, then T ⪰ A #_t B will be a convex set.
function vexity(
    vex::Tuple{<:Vexity,<:Vexity,<:Vexity},
    ::GenericConstraint{GeomMeanEpiConeSquare},
)
    T, A, B = vex

    # NOTE: can't say A == NotDcp() because the NotDcp constructor prints a
    # warning message.
    if typeof(A) == ConcaveVexity || typeof(A) == NotDcp
        return NotDcp()
    elseif typeof(B) == ConcaveVexity || typeof(B) == NotDcp
        return NotDcp()
    end
    output_vex = ConvexVexity() + (-T)
    if output_vex == ConcaveVexity()
        return NotDcp()
    end
    return output_vex
end

function _add_constraint!(
    context::Context,
    constraint::GenericConstraint{GeomMeanEpiConeSquare},
)
    T, A, B = constraint.child
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
