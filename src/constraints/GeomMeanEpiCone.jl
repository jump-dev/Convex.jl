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
mutable struct GeomMeanEpiCone
    A::AbstractExpr
    B::AbstractExpr
    t::Rational
    size::Tuple{Int,Int}

    function GeomMeanEpiCone(
        A::AbstractExpr,
        B::AbstractExpr,
        t::Rational,
        fullhyp::Bool = true,
    )
        if size(A) != size(B)
            throw(DimensionMismatch("A and B must be the same size"))
        end
        n = size(A)[1]
        if size(A) != (n, n)
            throw(DimensionMismatch("A and B must be square"))
        end
        if t < -1 || (t > 0 && t < 1) || t > 2
            throw(DomainError(t, "t must be in the range [-1, 0] or [1, 2]"))
        end
        return new(A, B, t, (n, n))
    end

    function GeomMeanEpiCone(
        A::Value,
        B::AbstractExpr,
        t::Rational,
        fullhyp::Bool = true,
    )
        return GeomMeanEpiCone(constant(A), B, t, fullhyp)
    end

    function GeomMeanEpiCone(
        A::AbstractExpr,
        B::Value,
        t::Rational,
        fullhyp::Bool = true,
    )
        return GeomMeanEpiCone(A, constant(B), t, fullhyp)
    end

    function GeomMeanEpiCone(
        A::Value,
        B::Value,
        t::Rational,
        fullhyp::Bool = true,
    )
        return GeomMeanEpiCone(constant(A), constant(B), t, fullhyp)
    end

    function GeomMeanEpiCone(
        A::Union{AbstractExpr,Value},
        B::Union{AbstractExpr,Value},
        t::Integer,
        fullhyp::Bool = true,
    )
        return GeomMeanEpiCone(A, B, t // 1, fullhyp)
    end
end

mutable struct GeomMeanEpiConeConstraint <: Constraint
    T::AbstractExpr
    cone::GeomMeanEpiCone

    function GeomMeanEpiConeConstraint(T::AbstractExpr, cone::GeomMeanEpiCone)
        if size(T) != cone.size
            throw(DimensionMismatch("T must be size $(cone.size)"))
        end
        return new(T, cone)
    end

    function GeomMeanEpiConeConstraint(T::Value, cone::GeomMeanEpiCone)
        return GeomMeanEpiConeConstraint(constant(T), cone)
    end
end

head(io::IO, ::GeomMeanEpiConeConstraint) = print(io, "∈(GeomMeanEpiCone)")

Base.in(T, cone::GeomMeanEpiCone) = GeomMeanEpiConeConstraint(T, cone)

function AbstractTrees.children(constraint::GeomMeanEpiConeConstraint)
    return (
        constraint.T,
        constraint.cone.A,
        constraint.cone.B,
        "t=$(constraint.cone.t)",
    )
end

# For t ∈ [-1, 0] ∪ [1, 2] the t-weighted matrix geometric mean is matrix convex
# (arxiv:1512.03401).
# So if A and B are convex sets, then T ⪰ A #_t B will be a convex set.
function vexity(constraint::GeomMeanEpiConeConstraint)
    A = vexity(constraint.cone.A)
    B = vexity(constraint.cone.B)
    T = vexity(constraint.T)

    # NOTE: can't say A == NotDcp() because the NotDcp constructor prints a
    # warning message.
    if typeof(A) == ConcaveVexity || typeof(A) == NotDcp
        return NotDcp()
    elseif typeof(B) == ConcaveVexity || typeof(B) == NotDcp
        return NotDcp()
    end
    vex = ConvexVexity() + (-T)
    if vex == ConcaveVexity()
        return NotDcp()
    end
    return vex
end

function _add_constraint!(
    context::Context,
    constraint::GeomMeanEpiConeConstraint,
)
    A = constraint.cone.A
    B = constraint.cone.B
    t = constraint.cone.t
    T = constraint.T
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
