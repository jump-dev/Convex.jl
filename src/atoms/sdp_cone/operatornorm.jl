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
    head::Symbol
    id_hash::UInt64
    children::Tuple{AbstractExpr}
    size::Tuple{Int, Int}

    function OperatorNormAtom(x::AbstractExpr)
        children = (x,)
        return new(:opnorm, hash(children), children, (1,1))
    end
end

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
    opnorm(evaluate(x.children[1]), 2)
end

sigmamax(x::AbstractExpr) = OperatorNormAtom(x)

function opnorm(x::AbstractExpr, p::Real=2)
    if length(size(x)) <= 1 || minimum(size(x)) == 1
        throw(ArgumentError("argument to `opnorm` must be a matrix"))
    end
    if p == 1
        return maximum(sum(abs(x), dims=1))
    elseif p == 2
        return OperatorNormAtom(x)
    elseif p == Inf
        return maximum(sum(abs(x), dims=2))
    else
        throw(ArgumentError("matrix p-norms only defined for p = 1, 2, and Inf"))
    end
end

Base.@deprecate operatornorm(x::AbstractExpr) opnorm(x)


function template(x::OperatorNormAtom, context::Context{T}) where T
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
        p = minimize(t, [t*sparse(1.0I, m, m) A; A' t*sparse(1.0I, n, n)] ⪰ 0)
        return template(p, context)
    else
        t = template(Variable(), context)
        f = operate(vcat, T, t, template(A, context))
        m, n = size(A)
        MOI_add_constraint(context.model, f, MOI.NormSpectralCone(m,n))
        return t
    end
end
