#############################################################################
# maximum.jl
# Compute the maximum value of an array.
# All expressions and atoms are subtpyes of AbstractExpr.
# Please read expressions.jl first.
#############################################################################
import Base.maximum

### Maximum Atom
mutable struct MaximumAtom <: AbstractExpr
    children::Tuple{AbstractExpr}
    size::Tuple{Int,Int}

    function MaximumAtom(x::AbstractExpr)
        if sign(x) == ComplexSign()
            error("Argument should be real instead it is $(sign(x))")
        else
            children = (x,)
            return new(children, (1, 1))
        end
    end
end

head(io::IO, ::MaximumAtom) = print(io, "maximum")

function sign(x::MaximumAtom)
    return sign(x.children[1])
end

# The monotonicity
function monotonicity(x::MaximumAtom)
    return (Nondecreasing(),)
end

# If we have h(x) = f o g(x), the chain rule says h''(x) = g'(x)^T f''(g(x))g'(x) + f'(g(x))g''(x);
# this represents the first term
function curvature(x::MaximumAtom)
    return ConvexVexity()
end

function evaluate(x::MaximumAtom)
    return Base.maximum(evaluate(x.children[1]))
end

function _conic_form!(context::Context, x::MaximumAtom)
    t = Variable()
    add_constraint!(context, t >= x.children[1])
    return conic_form!(context, t)
end

maximum(x::AbstractExpr) = MaximumAtom(x)
