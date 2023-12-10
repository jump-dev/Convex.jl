#############################################################################
# sum.jl
# Handles sums of expressions.
# All expressions and atoms are subtpyes of AbstractExpr.
# Please read expressions.jl first.
#############################################################################
import Base.sum

### Sum Atom
mutable struct SumAtom <: AbstractExpr
    children::Tuple{AbstractExpr}
    size::Tuple{Int,Int}

    function SumAtom(x::AbstractExpr)
        children = (x,)
        return new(children, (1, 1))
    end
end

head(io::IO, ::SumAtom) = print(io, "sum")

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

function new_conic_form!(context::Context{T}, A::SumAtom) where {T}
    subobj = conic_form!(context, only(children(A)))
    obj = operate(sum, T, sign(A), subobj)
    return obj
end

# Dispatch to an internal helper function that handles the dimension argument in
# the same manner as Base, with dims=: denoting a regular sum
sum(x::AbstractExpr; dims = :) = _sum(x, dims)

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

Base.@deprecate sum(x::AbstractExpr, dim::Int) sum(x, dims = dim)
