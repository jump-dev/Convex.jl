#############################################################################
# logsumexp.jl
# log of sum of exponentials of each entry in an expression
# All expressions and atoms are subtypes of AbstractExpr.
# Please read expressions.jl first.
#############################################################################

### LogSumExp

# TODO: make this work for a *list* of inputs, rather than just for vector/matrix inputs

mutable struct LogSumExpAtom <: AbstractExpr
    children::Tuple{AbstractExpr}
    size::Tuple{Int,Int}

    function LogSumExpAtom(x::AbstractExpr)
        if sign(x) == ComplexSign()
            error("The argument should be real but it's instead complex")
        else
            children = (x,)
            return new(children, (1, 1))
        end
    end
end

head(io::IO, ::LogSumExpAtom) = print(io, "logsumexp")

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
    _x = evaluate(x.children[1])
    max_x = maximum(_x)
    return max_x + log(sum(exp.(_x .- max_x)))
end

logsumexp(x::AbstractExpr) = LogSumExpAtom(x)

function _conic_form!(context::Context, e::LogSumExpAtom)
    # log(sum(exp(x))) <= t  <=>  sum(exp(x)) <= exp(t) <=> sum(exp(x - t)) <= 1
    t = Variable()
    z = sum(exp(e.children[1] - t * ones(size(e.children[1]))))
    objective = conic_form!(context, t)
    add_constraint!(context, z <= 1)
    return objective
end

function logisticloss(e::AbstractExpr)
    s = 0
    length(e) == 1 && return logsumexp([e; 0])
    for i in 1:length(e)
        s += logsumexp([e[i]; 0])
    end
    return s
end
