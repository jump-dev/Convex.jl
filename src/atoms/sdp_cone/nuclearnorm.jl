#############################################################################
# nuclearnorm.jl
# Handles nuclear norm (the sum of the singular values of a matrix),
# All expressions and atoms are subtypes of AbstractExpr.
# Please read expressions.jl first.
#############################################################################

### Nuclear norm

struct NuclearNormAtom <: AbstractExpr
    head::Symbol
    id_hash::UInt64
    children::Tuple{AbstractExpr}
    size::Tuple{Int, Int}

    function NuclearNormAtom(x::AbstractExpr)
        children = (x,)
        return new(:nuclearnorm, hash(children), children, (1,1))
    end
end

function sign(x::NuclearNormAtom)
    return Positive()
end

# The monotonicity
function monotonicity(x::NuclearNormAtom)
    return (NoMonotonicity(),)
end

function curvature(x::NuclearNormAtom)
    return ConvexVexity()
end

function evaluate(x::NuclearNormAtom)
    return sum(svdvals(evaluate(x.children[1])))
end

nuclearnorm(x::AbstractExpr) = NuclearNormAtom(x)


function template(x::NuclearNormAtom, context::Context{T}) where T
    A = only(children(x))
    if sign(A) == ComplexSign()
        # I'm not sure how to use MOI's `NormNuclearCone` in this case, so we'll just do the extended formulation as an SDP ourselves:
        #   minimize (tr(U) + tr(V))/2
        #   subject to
        #            [U A; A' V] ⪰ 0
        # see eg Recht, Fazel, Parillo 2008 "Guaranteed Minimum-Rank Solutions of Linear Matrix Equations via Nuclear Norm Minimization"
        # http://arxiv.org/pdf/0706.4138v1.pdf
        m, n = size(A)
        U = Variable(m,m)
        V = Variable(n,n)
        p = minimize((tr(U) + tr(V))/2, [U A; A' V] ⪰ 0)
        return template(p, context)
    else
        t = template(Variable(), context)
        f = operate(vcat, T, t, template(A, context))
        m, n = size(A)
        MOI_add_constraint(context.model, f, MOI.NormNuclearCone(m,n))
        return t
    end
end
