#############################################################################
# entropy.jl
# entropy (ie, sum_i( -x_i log (x_i) ) of an expression x
# All expressions and atoms are subtypes of AbstractExpr.
# Please read expressions.jl first.
#############################################################################

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
            error("The argument should be real but it's instead complex")
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

entropy(x::AbstractExpr) = sum(EntropyAtom(x))

function template(e::EntropyAtom, context::Context)
    # -x log x >= t  <=>  x exp(t/x) <= 1  <==>  (t,x,1) in exp cone
    t = Variable(e.size)
    x = e.children[1]
    add_constraints_to_context(ExpConstraint(t, x, ones(e.size...)), context)

    return template(t, context)
end
