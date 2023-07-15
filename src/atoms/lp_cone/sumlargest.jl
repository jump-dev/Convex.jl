#############################################################################
# sumlargest.jl
# The sum of the k largest/smallest elements of an expression
# All expressions and atoms are subtpyes of AbstractExpr.
# Please read expressions.jl first.
#############################################################################

mutable struct SumLargestAtom <: AbstractExpr
    children::Tuple{AbstractExpr}
    size::Tuple{Int,Int}
    k::Int

    function SumLargestAtom(x::AbstractExpr, k::Int)
        if sign(x) == ComplexSign()
            error("Argument should be real instead it is $(sign(x))")
        else
            if k <= 0
                error(
                    "sumlargest and sumsmallest only support positive values of k",
                )
            end
            if k > length(x)
                error("k cannot be larger than the number of entries in x")
            end
            children = (x,)
            return new(children, (1, 1), k)
        end
    end
end

head(io::IO, ::SumLargestAtom) = print(io, "sumlargest")

function sign(x::SumLargestAtom)
    return sign(x.children[1])
end

function monotonicity(x::SumLargestAtom)
    return (Nondecreasing(),)
end

function curvature(x::SumLargestAtom)
    return ConvexVexity()
end

function evaluate(x::SumLargestAtom)
    return sum(sort(vec(evaluate(x.children[1])), rev = true)[1:x.k])
end

function _conic_form!(context::Context, x::SumLargestAtom)
    c = x.children[1]
    t = Variable(size(c))
    q = Variable()
    # sum k largest given by the solution to
    # minimize sum(t) + k*q
    # subject to c <= t + q, t >= 0
    objective = conic_form!(context, sum(t) + x.k * q)
    add_constraint!(context, c <= t + q)
    add_constraint!(context, t >= 0)
    return objective
end

function sumlargest(x::AbstractExpr, k::Int)
    return k == 0 ? Constant(0) : SumLargestAtom(x, k)
end
function sumsmallest(x::AbstractExpr, k::Int)
    return k == 0 ? Constant(0) : -SumLargestAtom(-x, k)
end
