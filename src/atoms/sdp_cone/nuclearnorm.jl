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
    size::Tuple{Int,Int}

    function NuclearNormAtom(x::AbstractExpr)
        children = (x,)
        return new(:nuclearnorm, hash(children), children, (1, 1))
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

# Create the equivalent conic problem:
#   minimize (tr(U) + tr(V))/2
#   subject to
#            [U A; A' V] ⪰ 0
# see eg Recht, Fazel, Parillo 2008 "Guaranteed Minimum-Rank Solutions of Linear Matrix Equations via Nuclear Norm Minimization"
# http://arxiv.org/pdf/0706.4138v1.pdf
#
# The complex case is example 1.20 of Watrous' "The Theory of Quantum Information"
# (the operator A is negated but this doesn't affect the norm)
# https://cs.uwaterloo.ca/~watrous/TQI/TQI.pdf
function conic_form!(x::NuclearNormAtom, unique_conic_forms)
    if !has_conic_form(unique_conic_forms, x)
        A = x.children[1]
        m, n = size(A)
        if sign(A) == ComplexSign()
            U = ComplexVariable(m, m)
            V = ComplexVariable(n, n)
        else
            U = Variable(m, m)
            V = Variable(n, n)
        end
        p = minimize(0.5 * real(tr(U) + tr(V)), [U A; A' V] ⪰ 0)
        cache_conic_form!(unique_conic_forms, x, p)
    end
    return get_conic_form(unique_conic_forms, x)
end
