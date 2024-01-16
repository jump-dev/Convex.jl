mutable struct MatrixFracAtom <: AbstractExpr
    children::Tuple{AbstractExpr,AbstractExpr}
    size::Tuple{Int,Int}

    function MatrixFracAtom(x::AbstractExpr, P::AbstractExpr)
        if x.size[2] != 1
            error(
                "[MatrixFracAtom] first argument of matrixfrac must be a vector",
            )
        elseif P.size[1] != P.size[2]
            error(
                "[MatrixFracAtom] second argument of matrixfrac must be square",
            )
        elseif x.size[1] != P.size[1]
            error(
                "[MatrixFracAtom] sizes must agree for arguments of matrixfrac",
            )
        end
        return new((x, P), (1, 1))
    end
end

head(io::IO, ::MatrixFracAtom) = print(io, "matrixfrac")

Base.sign(::MatrixFracAtom) = Positive()

monotonicity(::MatrixFracAtom) = (NoMonotonicity(), NoMonotonicity())

curvature(::MatrixFracAtom) = ConvexVexity()

function evaluate(m::MatrixFracAtom)
    x = evaluate(m.children[1])
    return x' * (evaluate(m.children[2]) \ x)
end

matrixfrac(x::AbstractExpr, P::AbstractExpr) = MatrixFracAtom(x, P)

matrixfrac(x::Value, P::AbstractExpr) = MatrixFracAtom(constant(x), P)

matrixfrac(x::AbstractExpr, P::Value) = MatrixFracAtom(x, constant(P))

function new_conic_form!(context::Context, m::MatrixFracAtom)
    x, P = m.children
    t = Variable()
    # the matrix [t x'; x P] has Schur complement t - x'*P^{-1}*x
    # this matrix is PSD <=> t >= x'*P^{-1}*x
    add_constraint!(context, [t x'; x P] ⪰ 0)
    return conic_form!(context, t)
end
