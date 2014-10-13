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
  children_hash::Uint64
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

function conic_form(e::LogSumExpAtom, unique_constr)
  if !((e.head, e.children_hash) in keys(unique_constr))
    # log(sum(exp(z))) <= t  <=>  sum(exp(z)) <= exp(t)
    t = Variable()
    objective, constraints = conic_form(e.children[1], unique_constr)
    _, new_constraints = conic_form(sum(exp(e.children[1])) <= exp(t), unique_constr)
    append!(constraints, new_constraints)

    unique_constr[(e.head, e.children_hash)] = (objective, constraints)
  end
  return safe_copy(unique_constr[(e.head, e.children_hash)])
end

logistic_loss(e::AbstractExpr) = logsumexp([e, 0])
