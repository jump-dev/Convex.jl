#############################################################################
# exp.jl
# e raised to the power of an expression
# All expressions and atoms are subtpyes of AbstractExpr.
# Please read expressions.jl first.
#############################################################################

import Base.exp
export exp
export sign, curvature, monotonicity, evaluate

### Exponential

struct ExpAtom <: AbstractExpr
    head::Symbol
    id_hash::UInt64
    children::Tuple{AbstractExpr}
    size::Tuple{Int, Int}

    function ExpAtom(x::AbstractExpr)
        if sign(x) == ComplexSign()
            error("argument should be real but it's instead complex")
        else
            children = (x,)
            return new(:exp, hash(children), children, x.size)
        end
    end
end

function sign(x::ExpAtom)
    return Positive()
end

function monotonicity(x::ExpAtom)
    return (Nondecreasing(),)
end

function curvature(x::ExpAtom)
    return ConvexVexity()
end

function evaluate(x::ExpAtom)
    return exp.(evaluate(x.children[1]))
end

exp(x::AbstractExpr) = ExpAtom(x)

function conic_form!(e::ExpAtom, unique_conic_forms::UniqueConicForms=UniqueConicForms())
    if !has_conic_form(unique_conic_forms, e)
        # exp(x) \leq z  <=>  (x,ones(),z) \in ExpCone
        x = e.children[1]
        y = Constant(ones(size(x)))
        z = Variable(size(x))
        objective = conic_form!(z, unique_conic_forms)
        conic_form!(ExpConstraint(x, y, z), unique_conic_forms)
        cache_conic_form!(unique_conic_forms, e, objective)
    end
    return get_conic_form(unique_conic_forms, e)
end
