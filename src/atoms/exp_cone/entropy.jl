#############################################################################
# entropy.jl
# entropy (ie, sum_i( -x_i log (x_i) ) of an expression x
# All expressions and atoms are subtypes of AbstractExpr.
# Please read expressions.jl first.
#############################################################################

export entropy
export sign, curvature, monotonicity, evaluate

### Entropy: sum_i -x_i log (x_i)

# TODO: make this work for a *list* of inputs, rather than just for scalar/vector/matrix inputs

# Entropy atom: -xlogx entrywise
mutable struct EntropyAtom <: AbstractExpr
    head::Symbol
    id_hash::UInt64
    children::Tuple{AbstractExpr}
    size::Tuple{Int, Int}

    function EntropyAtom(x::AbstractExpr)
        if sign(x) == ComplexSign()
            error("argument should be real but it's instead complex")
        else
            children = (x,)
            # TODO check positivity or enforce it
            return new(:entropy, hash(children), children, size(x))
        end
    end
end

function sign(x::EntropyAtom)
    return NoSign()
end

function monotonicity(x::EntropyAtom)
    return (NoMonotonicity(),)
end

function curvature(x::EntropyAtom)
    return ConcaveVexity()
end

function evaluate(x::EntropyAtom)
    c = evaluate(x.children[1])
    return -c .* log.(c)
end

function conic_form!(e::EntropyAtom, unique_conic_forms::UniqueConicForms=UniqueConicForms())
    if !has_conic_form(unique_conic_forms, e)
        # -x log x >= t  <=>  x exp(t/x) <= 1  <==>  (t,x,1) in exp cone
        t = Variable(e.size)
        x = e.children[1]
        objective = conic_form!(t, unique_conic_forms)
        for i = 1:size(x, 1)
            for j = 1:size(x, 2)
                conic_form!(ExpConstraint(t[i, j], x[i, j], 1), unique_conic_forms)
            end
        end
        cache_conic_form!(unique_conic_forms, e, objective)
    end
    return get_conic_form(unique_conic_forms, e)
end

entropy(x::AbstractExpr) = sum(EntropyAtom(x))
