mutable struct OperatorNormAtom <: AbstractExpr
    children::Tuple{AbstractExpr}
    size::Tuple{Int,Int}

    OperatorNormAtom(x::AbstractExpr) = new((x,), (1, 1))
end

head(io::IO, ::OperatorNormAtom) = print(io, "opnorm")

Base.sign(::OperatorNormAtom) = Positive()

monotonicity(::OperatorNormAtom) = (NoMonotonicity(),)

curvature(::OperatorNormAtom) = ConvexVexity()

# In julia, `norm` on matrices is the operator norm.
function evaluate(x::OperatorNormAtom)
    return LinearAlgebra.opnorm(evaluate(x.children[1]), 2)
end

sigmamax(x::AbstractExpr) = OperatorNormAtom(x)

function LinearAlgebra.opnorm(x::AbstractExpr, p::Real = 2)
    if length(size(x)) <= 1 || minimum(size(x)) == 1
        throw(ArgumentError("argument to `opnorm` must be a matrix"))
    end
    if p == 1
        return maximum(sum(abs(x), dims = 1))
    elseif p == 2
        return OperatorNormAtom(x)
    elseif p == Inf
        return maximum(sum(abs(x), dims = 2))
    else
        throw(
            ArgumentError("matrix p-norms only defined for p = 1, 2, and Inf"),
        )
    end
end

Base.@deprecate operatornorm(x::AbstractExpr) LinearAlgebra.opnorm(x)

function new_conic_form!(context::Context{T}, x::OperatorNormAtom) where {T}
    A = x.children[1]
    if sign(A) == ComplexSign()
        # Create the equivalent conic problem:
        #   minimize t
        #   subject to
        #            [tI_m A; A' tI_n] is positive semidefinite
        # see eg Recht, Fazel, Parillo 2008 "Guaranteed Minimum-Rank Solutions
        # of Linear Matrix Equations via Nuclear Norm Minimization"
        # http://arxiv.org/pdf/0706.4138v1.pdf
        m, n = size(A)
        t = Variable()
        p = minimize(t, [t*spidentity(T, m) A; A' t*spidentity(T, n)] ⪰ 0)
        return conic_form!(context, p)
    end
    t = conic_form!(context, Variable())
    f = operate(vcat, T, sign(x), t, conic_form!(context, A))
    m, n = size(A)
    MOI_add_constraint(context.model, f, MOI.NormSpectralCone(m, n))
    return t
end
