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

Note on fullhyp: GeomMeanEpiCone will always return a full epigraph cone
(unlike GeomMeanHypoCone) and so this parameter is not really used.  It is
here just for consistency with the GeomMeanHypoCone function.

REFERENCE
  Ported from CVXQUAD which is based on the paper: "Lieb's concavity
  theorem, matrix geometric means and semidefinite optimization" by Hamza
  Fawzi and James Saunderson (arXiv:1512.03401)
"""
struct GeomMeanEpiCone <: MOI.AbstractVectorSet
    side_dimension::Int
    t::Rational

    function GeomMeanEpiCone(
        side_dimension::Int,
        t::Rational,
    )
        if t < -1 || (t > 0 && t < 1) || t > 2
            throw(DomainError(t, "t must be in the range [-1, 0] or [1, 2]"))
        end
        return new(side_dimension, t)
    end
end

head(io::IO, ::GeomMeanEpiCone) = print(io, "∈(GeomMeanEpiCone)")

Base.in(T, cone::GeomMeanEpiCone) = GenericConstraint{GeomMeanEpiCone}(T, cone)

# For t ∈ [-1, 0] ∪ [1, 2] the t-weighted matrix geometric mean is matrix convex
# (arxiv:1512.03401).
# So if A and B are convex sets, then T ⪰ A #_t B will be a convex set.
function vexity(constraint::GenericConstraint{GeomMeanEpiCone})
    n = constraint.set.side_dimension
    T = constraint.child[1:n^2]
    A = constraint.child[n^2 .+ (1:n^2)]
    B = constraint.child[2n^2 .+ (1:n^2)]
    vT = vexity(T)
    vA = vexity(A)
    vB = vexity(B)

    # NOTE: can't say A == NotDcp() because the NotDcp constructor prints a
    # warning message.
    if typeof(vA) == ConcaveVexity || typeof(vA) == NotDcp
        return NotDcp()
    elseif typeof(vB) == ConcaveVexity || typeof(vB) == NotDcp
        return NotDcp()
    end
    vex = ConvexVexity() + (-vT)
    if vex == ConcaveVexity()
        return NotDcp()
    end
    return vex
end

function _add_constraint!(
    context::Context,
    constraint::GenericConstraint{GeomMeanEpiCone},
)
    n = constraint.set.side_dimension
    T = constraint.child[1:n^2]
    A = constraint.child[n^2 .+ (1:n^2)]
    B = constraint.child[2n^2 .+ (1:n^2)]
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
        @assert t >= 1 # range of t checked in GeomMeanEpiCone constructor
        add_constraint!(context, [T B; B Z] ⪰ 0)
        add_constraint!(context, Z in GeomMeanHypoCone(A, B, 2 - t, false))
    end
    return
end
