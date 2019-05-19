#############################################################################
# logsumexp.jl
# log of sum of exponentials of each entry in an expression
# All expressions and atoms are subtypes of AbstractExpr.
# Please read expressions.jl first.
#############################################################################

export logsumexp, logisticloss
export sign, curvature, monotonicity, evaluate

### LogSumExp

# TODO: make this work for a *list* of inputs, rather than just for vector/matrix inputs

struct LogSumExpAtom <: AbstractExpr
    head::Symbol
    id_hash::UInt64
    children::Tuple{AbstractExpr}
    size::Tuple{Int, Int}

    function LogSumExpAtom(x::AbstractExpr)
        if sign(x) == ComplexSign()
            error("argument should be real but it's instead complex")
        else
            children = (x,)
            return new(:logsumexp, hash(children), children, (1,1))
        end
    end
end

function sign(x::LogSumExpAtom)
    return NoSign()
end

function monotonicity(x::LogSumExpAtom)
    return (Nondecreasing(),)
end

function curvature(x::LogSumExpAtom)
    return ConvexVexity()
end

function evaluate(x::LogSumExpAtom)
    return log(sum(exp.(evaluate(x.children[1]))))
end

logsumexp(x::AbstractExpr) = LogSumExpAtom(x)

function conic_form!(e::LogSumExpAtom, unique_conic_forms::UniqueConicForms=UniqueConicForms())
    if !has_conic_form(unique_conic_forms, e)
        # log(sum(exp(x))) <= t  <=>  sum(exp(x)) <= exp(t) <=> sum(exp(x - t)) <= 1
        t = Variable(e.size)
        z = sum(exp(e.children[1] - t))
        objective = conic_form!(t, unique_conic_forms)
        conic_form!(z, unique_conic_forms)
        conic_form!(z <= 1, unique_conic_forms)

        cache_conic_form!(unique_conic_forms, e, objective)
    end
    return get_conic_form(unique_conic_forms, e)
end

function logisticloss(e::AbstractExpr)
    s = 0
    length(e) == 1 && return logsumexp([e, 0])
    for i = 1:length(e)
        s += logsumexp([e[i]; 0])
    end
    return s
end
