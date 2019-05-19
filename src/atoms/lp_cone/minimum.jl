#############################################################################
# minimum.jl
# Compute the minimum value of an array.
# All expressions and atoms are subtpyes of AbstractExpr.
# Please read expressions.jl first.
#############################################################################
import Base.minimum
export minimum

### Minimum Atom
struct MinimumAtom <: AbstractExpr
    head::Symbol
    id_hash::UInt64
    children::Tuple{AbstractExpr}
    size::Tuple{Int, Int}

    function MinimumAtom(x::AbstractExpr)
        if sign(x) == ComplexSign()
            error("argument should be real instead it is $(sign(x))")
        else
            children = (x,)
            return new(:minimum, hash(children), children, (1, 1))
        end
    end
end

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

# x >= this if minimum(x) = this
# so, x - this will be in the :NonNeg cone
function conic_form!(x::MinimumAtom, unique_conic_forms::UniqueConicForms=UniqueConicForms())
    if !has_conic_form(unique_conic_forms, x)
        this = Variable()
        objective = conic_form!(this, unique_conic_forms)
        conic_form!(this <= x.children[1], unique_conic_forms)
        cache_conic_form!(unique_conic_forms, x, objective)
    end
    return get_conic_form(unique_conic_forms, x)
end

minimum(x::AbstractExpr) = MinimumAtom(x)
