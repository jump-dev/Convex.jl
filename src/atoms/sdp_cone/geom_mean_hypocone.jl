#############################################################################
# geom_mean_hypocone.jl
# Constrains T to
#   A #_t B ⪰ T
# where:
#   * A #_t B is the t-weighted geometric mean of A and B:
#     A^{1/2} (A^{-1/2} B A^{-1/2})^t A^{1/2}
#     Parameter t should be in [0,1].
#   * Constraints A ⪰ 0, B ⪰ 0 are added.
#
# Note on parameter fullhyp:
#   In many applications one doesn't need the full hypograph
#     hyp_t = {(A,B,T) : A #_t B ⪰ T}
#   but rather it is enough to work with a convex set C_t that satisfies
#     (A,B,A #_t B) \in C_t
#     (A,B,T) \in C_t  =>  A #_t B ⪰ T
#   In this case one should set fullhyp = false. The SDP description will be
#   (slightly) smaller. (By default fullhyp is set to true).
#
# All expressions and atoms are subtypes of AbstractExpr.
# Please read expressions.jl first.
#
#REFERENCE
#   Ported from CVXQUAD which is based on the paper: "Lieb's concavity
#   theorem, matrix geometric means and semidefinite optimization" by Hamza
#   Fawzi and James Saunderson (arXiv:1512.03401)
#############################################################################

mutable struct GeomMeanHypoCone
    A::AbstractExpr
    B::AbstractExpr
    t::Rational
    size::Tuple{Int,Int}
    fullhyp::Bool

    function GeomMeanHypoCone(
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
        if t < 0 || t > 1
            throw(DomainError(t, "t must be in the range [0, 1]"))
        end
        return new(A, B, t, (n, n), fullhyp)
    end

    function GeomMeanHypoCone(
        A::Value,
        B::AbstractExpr,
        t::Rational,
        fullhyp::Bool = true,
    )
        return GeomMeanHypoCone(constant(A), B, t, fullhyp)
    end
    function GeomMeanHypoCone(
        A::AbstractExpr,
        B::Value,
        t::Rational,
        fullhyp::Bool = true,
    )
        return GeomMeanHypoCone(A, constant(B), t, fullhyp)
    end
    function GeomMeanHypoCone(
        A::Value,
        B::Value,
        t::Rational,
        fullhyp::Bool = true,
    )
        return GeomMeanHypoCone(constant(A), constant(B), t, fullhyp)
    end

    function GeomMeanHypoCone(
        A::AbstractExprOrValue,
        B::AbstractExprOrValue,
        t::Integer,
        fullhyp::Bool = true,
    )
        return GeomMeanHypoCone(A, B, t // 1, fullhyp)
    end
end

mutable struct GeomMeanHypoConeConstraint <: Constraint
    T::AbstractExpr
    cone::GeomMeanHypoCone

    function GeomMeanHypoConeConstraint(T::AbstractExpr, cone::GeomMeanHypoCone)
        if size(T) != cone.size
            throw(DimensionMismatch("T must be size $(cone.size)"))
        end
        return new(T, cone)
    end

    function GeomMeanHypoConeConstraint(T::Value, cone::GeomMeanHypoCone)
        return GeomMeanHypoConeConstraint(constant(T), cone)
    end
end

head(io::IO, ::GeomMeanHypoConeConstraint) = print(io, "∈(GeomMeanHypoCone)")

in(T, cone::GeomMeanHypoCone) = GeomMeanHypoConeConstraint(T, cone)

function AbstractTrees.children(constraint::GeomMeanHypoConeConstraint)
    return (
        constraint.T,
        constraint.cone.A,
        constraint.cone.B,
        "t=$(constraint.cone.t)",
    )
end

# For t ∈ [0,1] the t-weighted matrix geometric mean is matrix concave (arxiv:1512.03401).
# So if A and B are convex sets, then T ⪯ A #_t B will be a convex set.
function vexity(constraint::GeomMeanHypoConeConstraint)
    A = vexity(constraint.cone.A)
    B = vexity(constraint.cone.B)
    T = vexity(constraint.T)

    # NOTE: can't say A == NotDcp() because the NotDcp constructor prints a warning message.
    if typeof(A) == ConvexVexity || typeof(A) == NotDcp
        return NotDcp()
    end
    if typeof(B) == ConvexVexity || typeof(B) == NotDcp
        return NotDcp()
    end
    # Copied from vexity(c::GtConstraint)
    vex = -ConcaveVexity() + T
    if vex == ConcaveVexity()
        vex = NotDcp()
    end
    return vex
end

function _add_constraint!(
    context::Context,
    constraint::GeomMeanHypoConeConstraint,
)
    A = constraint.cone.A
    B = constraint.cone.B
    t = constraint.cone.t
    T = constraint.T
    fullhyp = constraint.cone.fullhyp

    is_complex =
        sign(A) == ComplexSign() ||
        sign(B) == ComplexSign() ||
        sign(T) == ComplexSign()
    if is_complex
        make_temporary = () -> HermitianSemidefinite(size(A)[1])
    else
        make_temporary = () -> Semidefinite(size(A)[1])
    end

    if fullhyp && t != 0 && t != 1
        W = make_temporary()
        add_constraint!(context, W in GeomMeanHypoCone(A, B, t, false))
        add_constraint!(context, W ⪰ T)
    else
        p = t.num
        q = t.den

        if t == 0
            #println("geom_mean_hypocone p=$p q=$q t==0")
            add_constraint!(context, A ⪰ 0)
            add_constraint!(context, B ⪰ 0)
            add_constraint!(context, A ⪰ T)
        elseif t == 1
            #println("geom_mean_hypocone p=$p q=$q t==1")
            add_constraint!(context, A ⪰ 0)
            add_constraint!(context, B ⪰ 0)
            add_constraint!(context, B ⪰ T)
        elseif t == 1 // 2
            #println("geom_mean_hypocone p=$p q=$q t==1/2")
            add_constraint!(context, [A T; T' B] ⪰ 0)
        elseif ispow2(q)
            #println("geom_mean_hypocone p=$p q=$q ispow2(q)")
            Z = make_temporary()
            if t < 1 / 2
                add_constraint!(
                    context,
                    Z in GeomMeanHypoCone(A, B, 2 * t, false),
                )
                add_constraint!(context, [A T; T' Z] ⪰ 0)
            else
                add_constraint!(
                    context,
                    Z in GeomMeanHypoCone(A, B, 2 * t - 1, false),
                )
                add_constraint!(context, [B T; T' Z] ⪰ 0)
            end
        elseif ispow2(p) && t > 1 // 2
            #println("geom_mean_hypocone p=$p q=$q ispow2(p) && t>1/2")
            Z = make_temporary()
            add_constraint!(
                context,
                Z in GeomMeanHypoCone(A, T, (2 * p - q) // p, false),
            )
            add_constraint!(context, [Z T; T B] ⪰ 0)
        elseif t < 1 / 2
            #println("geom_mean_hypocone p=$p q=$q t<1/2")
            X = make_temporary()
            # Decompose t = (p/2^l) * (2^l/q) where l=floor(log2(q))
            l = floor(Int, log2(q))
            add_constraint!(
                context,
                X in GeomMeanHypoCone(A, B, p // (2^l), false),
            )
            add_constraint!(
                context,
                T in GeomMeanHypoCone(A, X, (2^l) // q, false),
            )
        else
            #println("geom_mean_hypocone p=$p q=$q else")
            add_constraint!(context, T in GeomMeanHypoCone(B, A, 1 - t, false))
        end
    end
end
