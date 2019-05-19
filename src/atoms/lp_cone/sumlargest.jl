#############################################################################
# sumlargest.jl
# The sum of the k largest/smallest elements of an expression
# All expressions and atoms are subtpyes of AbstractExpr.
# Please read expressions.jl first.
#############################################################################

export sumlargest, sumsmallest
export sign, curvature, monotonicity, evaluate

struct SumLargestAtom <: AbstractExpr
    head::Symbol
    id_hash::UInt64
    children::Tuple{AbstractExpr}
    size::Tuple{Int, Int}
    k::Int

    function SumLargestAtom(x::AbstractExpr, k::Int)
        if sign(x) == ComplexSign()
            error("argument should be real instead it is $(sign(x))")
        else
            if k <= 0
                error("sumlargest and sumsmallest only support positive values of k")
            end
            if k > length(x)
                error("k cannot be larger than the number of entries in x")
            end
            children = (x,)
            return new(:sumlargest, hash((children, k)), children, (1,1), k)
        end
    end
end

function sign(x::SumLargestAtom)
    return sign(x.children[1])
end

function monotonicity(x::SumLargestAtom)
    return (Nondecreasing(), )
end

function curvature(x::SumLargestAtom)
    return ConvexVexity()
end

function evaluate(x::SumLargestAtom)
    return sum(sort(vec(evaluate(x.children[1])), rev=true)[1:x.k])
end

function conic_form!(x::SumLargestAtom, unique_conic_forms::UniqueConicForms=UniqueConicForms())
    if !has_conic_form(unique_conic_forms, x)
        c = x.children[1]
        t = Variable(size(c))
        q = Variable()
        # sum k largest given by the solution to
        # minimize sum(t) + k*q
        # subject to c <= t + q, t >= 0
        objective = conic_form!(sum(t) + x.k*q, unique_conic_forms)
        conic_form!(c <= t + q, unique_conic_forms)
        conic_form!(t >= 0, unique_conic_forms)
        cache_conic_form!(unique_conic_forms, x, objective)
    end
    return get_conic_form(unique_conic_forms, x)
end

sumlargest(x::AbstractExpr, k::Int) = SumLargestAtom(x, k)
sumsmallest(x::AbstractExpr, k::Int) = -SumLargestAtom(-x, k)
