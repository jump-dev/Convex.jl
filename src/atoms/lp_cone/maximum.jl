#############################################################################
# maximum.jl
# Compute the maximum value of an array.
# All expressions and atoms are subtpyes of AbstractExpr.
# Please read expressions.jl first.
#############################################################################
import Base.maximum
export maximum

### Maximum Atom
struct MaximumAtom <: AbstractExpr
    head::Symbol
    id_hash::UInt64
    children::Tuple{AbstractExpr}
    size::Tuple{Int, Int}

    function MaximumAtom(x::AbstractExpr)
        if sign(x)==ComplexSign()
            error("argument should be real instead it is $(sign(x))")
        else
            children = (x,)
            return new(:maximum, hash(children), children, (1, 1))
        end
    end
end

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

# x <= this if maximum(x) = this
# so, this - x will be in the :NonNeg cone
function conic_form!(x::MaximumAtom, unique_conic_forms::UniqueConicForms=UniqueConicForms())
    if !has_conic_form(unique_conic_forms, x)
        this = Variable()
        objective = conic_form!(this, unique_conic_forms)
        conic_form!(this >= x.children[1], unique_conic_forms)
        cache_conic_form!(unique_conic_forms, x, objective)
    end
    return get_conic_form(unique_conic_forms, x)
end

maximum(x::AbstractExpr) = MaximumAtom(x)
