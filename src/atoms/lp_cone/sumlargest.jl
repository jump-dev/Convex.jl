#############################################################################
# sumlargest.jl
# The sum of the k largest/smallest elements of an expression
# All expressions and atoms are subtpyes of AbstractExpr.
# Please read expressions.jl first.
#############################################################################

struct SumLargestAtom <: AbstractExpr
    head::Symbol
    id_hash::UInt64
    children::Tuple{AbstractExpr}
    size::Tuple{Int, Int}
    k::Int

    function SumLargestAtom(x::AbstractExpr, k::Int)
        if sign(x) == ComplexSign()
            error("Argument should be real instead it is $(sign(x))")
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

function template(x::SumLargestAtom, context::Context)
    c = x.children[1]
    t = Variable(size(c))
    q = Variable()
    # sum k largest given by the solution to
    # minimize sum(t) + k*q
    # subject to c <= t + q, t >= 0
    objective = template(sum(t) + x.k*q, context)
    add_constraints_to_context(c <= t + q, context)
    add_constraints_to_context(t >= 0, context)
    return objective
end

sumlargest(x::AbstractExpr, k::Int) = SumLargestAtom(x, k)
sumsmallest(x::AbstractExpr, k::Int) = -SumLargestAtom(-x, k)
