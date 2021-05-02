#############################################################################
# geom_mean_epicone.jl
# Constrains T to
#   A #_t B ⪯ T
# where:
#   * A #_t B is the t-weighted geometric mean of A and B:
#     A^{1/2} (A^{-1/2} B A^{-1/2})^t A^{1/2}
#     Parameter t should be in [-1, 0] or [1, 2].
#   * Constraints A ⪰ 0, B ⪰ 0 are added.
# All expressions and atoms are subtypes of AbstractExpr.
# Please read expressions.jl first.
#
# Note on fullhyp: GeomMeanEpiCone will always return a full epigraph cone
# (unlike GeomMeanHypoCone) and so this parameter is not really used.  It is
# here just for consistency with the GeomMeanHypoCone function.
#
#REFERENCE
#   Ported from CVXQUAD which is based on the paper: "Lieb's concavity
#   theorem, matrix geometric means and semidefinite optimization" by Hamza
#   Fawzi and James Saunderson (arXiv:1512.03401)
#############################################################################

struct GeomMeanEpiCone
    A::AbstractExpr
    B::AbstractExpr
    t::Rational
    size::Tuple{Int, Int}

    function GeomMeanEpiCone(A::AbstractExpr, B::AbstractExpr, t::Rational, fullhyp::Bool=true)
        if size(A) != size(B)
            throw(DimensionMismatch("A and B must be the same size"))
        end
        n = size(A)[1]
        if size(A) != (n, n)
            throw(DimensionMismatch("A and B must be square"))
        end
        if t < -1 || (t > 0 && t < 1) || t > 2
            error("t must be in the range [-1, 0] or [1, 2]")
        end
        return new(A, B, t, (n, n))
    end

    GeomMeanEpiCone(A::Value,        B::AbstractExpr, t::Rational, fullhyp::Bool=true) = GeomMeanEpiCone(Constant(A), B, t, fullhyp)
    GeomMeanEpiCone(A::AbstractExpr, B::Value,        t::Rational, fullhyp::Bool=true) = GeomMeanEpiCone(A, Constant(B), t, fullhyp)
    GeomMeanEpiCone(A::Value,        B::Value,        t::Rational, fullhyp::Bool=true) = GeomMeanEpiCone(Constant(A), Constant(B), t, fullhyp)
end

struct GeomMeanEpiConeConstraint <: Constraint
    head::Symbol
    id_hash::UInt64
    T::AbstractExpr
    cone::GeomMeanEpiCone

    function GeomMeanEpiConeConstraint(T::AbstractExpr, cone::GeomMeanEpiCone)
        if size(T) != cone.size
            throw(DimensionMismatch("T must be size $(cone.size)"))
        end
        id_hash = hash((cone.A, cone.B, cone.t, :GeomMeanEpiCone))
        return new(:GeomMeanEpiCone, id_hash, T, cone)
    end

    GeomMeanEpiConeConstraint(T::Value, cone::GeomMeanEpiCone) = GeomMeanEpiConeConstraint(Constant(T), cone)
end

in(T, cone::GeomMeanEpiCone) = GeomMeanEpiConeConstraint(T, cone)

function AbstractTrees.children(constraint::GeomMeanEpiConeConstraint)
    return (constraint.T, constraint.cone.A, constraint.cone.B, "t=$(constraint.cone.t)")
end

# For t ∈ [-1, 0] ∪ [1, 2] the t-weighted matrix geometric mean is matrix convex (arxiv:1512.03401).
# So if A and B are convex sets, then T ⪰ A #_t B will be a convex set.
function vexity(constraint::GeomMeanEpiConeConstraint)
    A = vexity(constraint.cone.A)
    B = vexity(constraint.cone.B)
    T = vexity(constraint.T)

    # NOTE: can't say A == NotDcp() because the NotDcp constructor prints a warning message.
    if typeof(A) == ConcaveVexity || typeof(A) == NotDcp
        return NotDcp()
    end
    if typeof(B) == ConcaveVexity || typeof(B) == NotDcp
        return NotDcp()
    end
    # Copied from vexity(c::GtConstraint)
    vex = ConvexVexity() + (-T)
    if vex == ConcaveVexity()
        vex = NotDcp()
    end
    return vex
end

function conic_form!(constraint::GeomMeanEpiConeConstraint, unique_conic_forms::UniqueConicForms)
    if !has_conic_form(unique_conic_forms, constraint)
        A = constraint.cone.A
        B = constraint.cone.B
        t = constraint.cone.t
        T = constraint.T
        is_complex = sign(A) == ComplexSign() || sign(B) == ComplexSign() || sign(T) == ComplexSign()
        if is_complex
            make_temporary = () -> HermitianSemidefinite(size(A)[1])
        else
            make_temporary = () -> Semidefinite(size(A)[1])
        end

        Z = make_temporary()

        if t <= 0
            conic_form!([T A; A Z] ⪰ 0, unique_conic_forms)
            conic_form!(Z in GeomMeanHypoCone(A, B, -t, false), unique_conic_forms)
        elseif t >= 1
            conic_form!([T B; B Z] ⪰ 0, unique_conic_forms)
            conic_form!(Z in GeomMeanHypoCone(A, B, 2-t, false), unique_conic_forms)
        else
            error("t must be in the range [-1, 0] or [1, 2]")
        end

        cache_conic_form!(unique_conic_forms, constraint, Array{Convex.ConicConstr,1}())
    end
    return get_conic_form(unique_conic_forms, constraint)
end
