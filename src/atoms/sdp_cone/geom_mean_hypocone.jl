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

struct GeomMeanHypoCone
    A::AbstractExpr
    B::AbstractExpr
    t::Rational
    size::Tuple{Int, Int}
    fullhyp::Bool

    function GeomMeanHypoCone(A::AbstractExpr, B::AbstractExpr, t::Rational, fullhyp::Bool=true)
        if size(A) != size(B)
            error("A and B must be the same size")
        end
        n = size(A)[1]
        if size(A) != (n, n)
            error("A and B must be square")
        end
        if t < 0 || t > 1
            error("t must be in the range [0, 1]")
        end
        return new(A, B, t, (n, n), fullhyp)
    end

    GeomMeanHypoCone(A::Value,        B::AbstractExpr, t::Rational, fullhyp::Bool=true) = GeomMeanHypoCone(Constant(A), B, t, fullhyp)
    GeomMeanHypoCone(A::AbstractExpr, B::Value,        t::Rational, fullhyp::Bool=true) = GeomMeanHypoCone(A, Constant(B), t, fullhyp)
    GeomMeanHypoCone(A::Value,        B::Value,        t::Rational, fullhyp::Bool=true) = GeomMeanHypoCone(Constant(A), Constant(B), t, fullhyp)
end

struct GeomMeanHypoConeConstraint <: Constraint
    head::Symbol
    id_hash::UInt64
    T::AbstractExpr
    cone::GeomMeanHypoCone

    function GeomMeanHypoConeConstraint(T::AbstractExpr, cone::GeomMeanHypoCone)
        if size(T) != cone.size
            error("T must be size $(cone.size)")
        end
        id_hash = hash((cone.A, cone.B, cone.t, :GeomMeanHypoCone))
        return new(:GeomMeanHypoCone, id_hash, T, cone)
    end

    GeomMeanHypoConeConstraint(T::Value, cone::GeomMeanHypoCone) = GeomMeanHypoConeConstraint(Constant(T), cone)
end

in(T, cone::GeomMeanHypoCone) = GeomMeanHypoConeConstraint(T, cone)

function AbstractTrees.children(constraint::GeomMeanHypoConeConstraint)
    return (constraint.T, constraint.cone.A, constraint.cone.B, "t=$(constraint.cone.t)")
end

# FIXME what is the meaning of the vexity of a variable? of a constraint?  Is this correct?
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

function conic_form!(constraint::GeomMeanHypoConeConstraint, unique_conic_forms::UniqueConicForms)
    if !has_conic_form(unique_conic_forms, constraint)
        A = constraint.cone.A
        B = constraint.cone.B
        t = constraint.cone.t
        T = constraint.T
        fullhyp = constraint.cone.fullhyp

        is_complex = sign(A) == ComplexSign() || sign(B) == ComplexSign() || sign(T) == ComplexSign()
        if is_complex
            make_temporary = () -> HermitianSemidefinite(size(A)[1])
        else
            make_temporary = () -> Semidefinite(size(A)[1])
        end

        if fullhyp && t != 0 && t != 1
            W = make_temporary()
            conic_form!(W in GeomMeanHypoCone(A, B, t, false), unique_conic_forms)
            conic_form!(W ⪰ T, unique_conic_forms)
            cache_conic_form!(unique_conic_forms, constraint, Array{Convex.ConicConstr,1}())
        else
            p = t.num
            q = t.den

            if t == 0
                #println("geom_mean_hypocone p=$p q=$q t==0")
                conic_form!(A ⪰ 0, unique_conic_forms)
                conic_form!(B ⪰ 0, unique_conic_forms)
                conic_form!(A ⪰ T, unique_conic_forms)
            elseif t == 1
                #println("geom_mean_hypocone p=$p q=$q t==1")
                conic_form!(A ⪰ 0, unique_conic_forms)
                conic_form!(B ⪰ 0, unique_conic_forms)
                conic_form!(B ⪰ T, unique_conic_forms)
            elseif t == 1//2
                #println("geom_mean_hypocone p=$p q=$q t==1/2")
                conic_form!([A T; T' B] ⪰ 0, unique_conic_forms)
            elseif ispow2(q)
                #println("geom_mean_hypocone p=$p q=$q ispow2(q)")
                Z = make_temporary()
                if t < 1/2
                    conic_form!(Z in GeomMeanHypoCone(A, B, 2*t, false), unique_conic_forms)
                    conic_form!([A T; T' Z] ⪰ 0, unique_conic_forms)
                else
                    conic_form!(Z in GeomMeanHypoCone(A, B, 2*t-1, false), unique_conic_forms)
                    conic_form!([B T; T' Z] ⪰ 0, unique_conic_forms)
                end
            elseif ispow2(p) && t > 1//2
                #println("geom_mean_hypocone p=$p q=$q ispow2(p) && t>1/2")
                Z = make_temporary()
                conic_form!(Z in GeomMeanHypoCone(A, T, (2*p-q)//p, false), unique_conic_forms)
                conic_form!([Z T; T B] ⪰ 0, unique_conic_forms)
            elseif t < 1/2
                #println("geom_mean_hypocone p=$p q=$q t<1/2")
                X = make_temporary()
                # Decompose t = (p/2^l) * (2^l/q) where l=floor(log2(q))
                l = floor(Int, log2(q));
                conic_form!(X in GeomMeanHypoCone(A, B, p//(2^l), false), unique_conic_forms)
                conic_form!(T in GeomMeanHypoCone(A, X, (2^l)//q, false), unique_conic_forms)
            else
                #println("geom_mean_hypocone p=$p q=$q else")
                conic_form!(T in GeomMeanHypoCone(B, A, 1-t, false), unique_conic_forms)
            end

            cache_conic_form!(unique_conic_forms, constraint, Array{Convex.ConicConstr,1}())
            #cache_conic_form!(unique_conic_forms, constraint, ConicConstr([A,B,T], :geom_mean_hypocone, [size(A), size(B), size(T)]))
            #cache_conic_form!(unique_conic_forms, conic_form!(T, unique_conic_forms))
        end
    end
    return get_conic_form(unique_conic_forms, constraint)
end
