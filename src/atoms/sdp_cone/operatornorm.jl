#############################################################################
# operatornorm.jl
# Handles matrix operator norm (the maximum singular value of a matrix)
# and creates the alias sigmamax
# All expressions and atoms are subtypes of AbstractExpr.
# Please read expressions.jl first.
#############################################################################
import LinearAlgebra: opnorm

### Operator norm

struct OperatorNormAtom <: AbstractExpr
    children::Tuple{AbstractExpr}
    size::Tuple{Int,Int}

    function OperatorNormAtom(x::AbstractExpr)
        children = (x,)
        return new(children, (1, 1))
    end
end

head(io::IO, ::OperatorNormAtom) = print(io, "opnorm")

function sign(x::OperatorNormAtom)
    return Positive()
end

# The monotonicity
function monotonicity(x::OperatorNormAtom)
    return (NoMonotonicity(),)
end

function curvature(x::OperatorNormAtom)
    return ConvexVexity()
end

# in julia, `norm` on matrices is the operator norm
function evaluate(x::OperatorNormAtom)
    return opnorm(evaluate(x.children[1]), 2)
end

sigmamax(x::AbstractExpr) = OperatorNormAtom(x)

function opnorm(x::AbstractExpr, p::Real = 2)
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

Base.@deprecate operatornorm(x::AbstractExpr) opnorm(x)

function conic_form!(context::Context{T}, x::OperatorNormAtom) where {T}
    A = x.children[1]
    if sign(A) == ComplexSign()
        # Create the equivalent conic problem:
        #   minimize t
        #   subject to
        #            [tI_m A; A' tI_n] is positive semidefinite
        # see eg Recht, Fazel, Parillo 2008 "Guaranteed Minimum-Rank Solutions of Linear Matrix Equations via Nuclear Norm Minimization"
        # http://arxiv.org/pdf/0706.4138v1.pdf
        m, n = size(A)
        t = Variable()
        p = minimize(
            t,
            [t*sparse(one(T) * I, m, m) A; A' t*sparse(one(T) * I, n, n)] ⪰ 0,
        )
        return conic_form!(context, p)
    else
        t = conic_form!(context, Variable())
        f = operate(vcat, T, sign(x), t, conic_form!(context, A))
        m, n = size(A)
        MOI_add_constraint(context.model, f, MOI.NormSpectralCone(m, n))
        return t
    end
end
