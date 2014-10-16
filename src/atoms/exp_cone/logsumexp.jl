#############################################################################
# logsumexp.jl
# log of sum of exponentials of each entry in an expression
# All expressions and atoms are subtypes of AbstractExpr.
# Please read expressions.jl first.
#############################################################################

export logsumexp, logistic_loss
export sign, curvature, monotonicity, evaluate

### LogSumExp

# TODO: make this work for a *list* of inputs, rather than just for vector/matrix inputs

type LogSumExpAtom <: AbstractExpr
  head::Symbol
  id_hash::Uint64
  children::(AbstractExpr,)
  size::(Int64, Int64)

  function LogSumExpAtom(x::AbstractExpr)
    children = (x,)
    return new(:logsumexp, hash(children), children, x.size)
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
  return logsumexp(evaluate(x.children[1]))
end

logsumexp(x::AbstractExpr) = LogSumExpAtom(x)

function conic_form(e::LogSumExpAtom, unique_conic_forms::UniqueConicForms)
  if !has_conic_form(unique_conic_forms, e)
    # log(sum(exp(z))) <= t  <=>  sum(exp(z)) <= exp(t)
    t = Variable()
    objective = conic_form(t, unique_conic_forms)
    conic_form(sum(exp(e.children[1])) <= exp(t), unique_conic_forms)

    cache_conic_form!(unique_conic_forms, e, objective)
  end
  return get_conic_form(unique_conic_forms, e)
end

logistic_loss(e::AbstractExpr) = logsumexp([e, 0])
