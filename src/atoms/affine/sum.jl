#############################################################################
# sum.jl
# Handles sums of expressions.
# All expressions and atoms are subtpyes of AbstractExpr.
# Please read expressions.jl first.
#############################################################################
import Base.sum

### Sum Atom
struct SumAtom <: AbstractExpr
    head::Symbol
    id_hash::UInt64
    children::Tuple{AbstractExpr}
    size::Tuple{Int, Int}

    function SumAtom(x::AbstractExpr)
        children = (x,)
        return new(:sum, hash(children), children, (1, 1))
    end
end

function sign(x::SumAtom)
    return sign(x.children[1])
end

# The monotonicity
function monotonicity(x::SumAtom)
    return (Nondecreasing(),)
end

# If we have h(x) = f o g(x), the chain rule says h''(x) = g'(x)^T f''(g(x))g'(x) + f'(g(x))g''(x);
# this represents the first term
function curvature(x::SumAtom)
    return ConstVexity()
end

function evaluate(x::SumAtom)
    return sum(evaluate(x.children[1]))
end


function template(A::SumAtom, context::Context{T}) where T
    subobj = template(only(children(A)), context)
    obj = operate(sum, T, subobj)
    return obj
end

# Dispatch to an internal helper function that handles the dimension argument in
# the same manner as Base, with dims=: denoting a regular sum
sum(x::AbstractExpr; dims=:) = _sum(x, dims)

_sum(x::AbstractExpr, ::Colon) = SumAtom(x)

function _sum(x::AbstractExpr, dimension::Integer)
    if dimension == 1
        return Constant(ones(1, x.size[1]), Positive()) * x
    elseif dimension == 2
        return x * Constant(ones(x.size[2], 1), Positive())
    else
        error("Sum not implemented for dimension $dimension")
    end
end

Base.@deprecate sum(x::AbstractExpr, dim::Int) sum(x, dims=dim)
