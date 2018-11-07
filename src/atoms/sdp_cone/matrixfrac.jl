#############################################################################
# matrixfrac.jl
# implements the atom for x^T*P^{-1}*x, where P is a positive semidefinite
# matrix.
# All expressions and atoms are subtypes of AbstractExpr.
# Please read expressions.jl first.
#############################################################################
export matrixfrac

struct MatrixFracAtom <: AbstractExpr
    head::Symbol
    id_hash::UInt64
    children::Tuple{AbstractExpr, AbstractExpr}
    size::Tuple{Int, Int}

    function MatrixFracAtom(x::AbstractExpr, P::AbstractExpr)
        if x.size[2] != 1
            error("first argument of matrix frac must be a vector")
        elseif P.size[1] != P.size[2]
            error("second argument of matrix frac must be square")
        elseif x.size[1] != P.size[1]
            error("sizes must agree for arguments of matrix frac")
        end
        children = (x, P)
        return new(:matrixfrac, hash(children), children, (1,1))
    end
end

function sign(m::MatrixFracAtom)
    return Positive()
end

function monotonicity(m::MatrixFracAtom)
    return (NoMonotonicity(), NoMonotonicity())
end

function curvature(m::MatrixFracAtom)
    return ConvexVexity()
end

function evaluate(m::MatrixFracAtom)
    x = evaluate(m.children[1])
    return x'*inv(evaluate(m.children[2]))*x
end

matrixfrac(x::AbstractExpr, P::AbstractExpr) = MatrixFracAtom(x, P)
matrixfrac(x::Value, P::AbstractExpr) = MatrixFracAtom(Constant(x), P)
matrixfrac(x::AbstractExpr, P::Value) = MatrixFracAtom(x, Constant(P))

function conic_form!(m::MatrixFracAtom, unique_conic_forms::UniqueConicForms=UniqueConicForms())
    if !has_conic_form(unique_conic_forms, m)
        x = m.children[1]
        P = m.children[2]
        t = Variable()
        # the matrix [t x'; x P] has Schur complement t - x'*P^{-1}*x
        # this matrix is PSD <=> t >= x'*P^{-1}*x
        p = minimize(t, [t x'; x P] ⪰ 0)
        cache_conic_form!(unique_conic_forms, m, p)
    end
    return get_conic_form(unique_conic_forms, m)
end
