#############################################################################
# minimum.jl
# Compute the minimum value of an array.
# All expressions and atoms are subtpyes of AbstractExpr.
# Please read expressions.jl first.
#############################################################################
import Base.minimum

### Minimum Atom
mutable struct MinimumAtom <: AbstractExpr
    children::Tuple{AbstractExpr}
    size::Tuple{Int,Int}

    function MinimumAtom(x::AbstractExpr)
        if sign(x) == ComplexSign()
            error("Argument should be real instead it is $(sign(x))")
        else
            children = (x,)
            return new(children, (1, 1))
        end
    end
end

head(io::IO, ::MinimumAtom) = print(io, "minimum")

function sign(x::MinimumAtom)
    return sign(x.children[1])
end

# The monotonicity
function monotonicity(x::MinimumAtom)
    return (Nondecreasing(),)
end

# If we have h(x) = f o g(x), the chain rule says h''(x) = g'(x)^T f''(g(x))g'(x) + f'(g(x))g''(x);
# this represents the first term
function curvature(x::MinimumAtom)
    return ConcaveVexity()
end

function evaluate(x::MinimumAtom)
    return Base.minimum(evaluate(x.children[1]))
end

function _conic_form!(context::Context, x::MinimumAtom)
    t = Variable()
    add_constraint!(context, t <= x.children[1])
    return conic_form!(context, t)
end

minimum(x::AbstractExpr) = MinimumAtom(x)
