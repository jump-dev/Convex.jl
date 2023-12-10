#############################################################################
# matrixfrac.jl
# implements the atom for x^T*P^{-1}*x, where P is a positive semidefinite
# matrix.
# All expressions and atoms are subtypes of AbstractExpr.
# Please read expressions.jl first.
#############################################################################

mutable struct MatrixFracAtom <: AbstractExpr
    children::Tuple{AbstractExpr,AbstractExpr}
    size::Tuple{Int,Int}

    function MatrixFracAtom(x::AbstractExpr, P::AbstractExpr)
        if x.size[2] != 1
            error("first argument of matrix frac must be a vector")
        elseif P.size[1] != P.size[2]
            error("second argument of matrix frac must be square")
        elseif x.size[1] != P.size[1]
            error("sizes must agree for arguments of matrix frac")
        end
        children = (x, P)
        return new(children, (1, 1))
    end
end

head(io::IO, ::MatrixFracAtom) = print(io, "matrixfrac")

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
    return x' * inv(evaluate(m.children[2])) * x
end

matrixfrac(x::AbstractExpr, P::AbstractExpr) = MatrixFracAtom(x, P)
matrixfrac(x::Value, P::AbstractExpr) = MatrixFracAtom(constant(x), P)
matrixfrac(x::AbstractExpr, P::Value) = MatrixFracAtom(x, constant(P))

function new_conic_form!(context::Context, m::MatrixFracAtom)
    x = m.children[1]
    P = m.children[2]
    t = Variable()
    # the matrix [t x'; x P] has Schur complement t - x'*P^{-1}*x
    # this matrix is PSD <=> t >= x'*P^{-1}*x
    p = minimize(t, [t x'; x P] ⪰ 0)
    return conic_form!(context, p)
end
