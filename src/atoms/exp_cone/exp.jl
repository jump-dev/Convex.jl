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

type ExpAtom <: AbstractExpr
  head::Symbol
  children_hash::Uint64
  children::(AbstractExpr,)
  size::(Int64, Int64)

  function ExpAtom(x::AbstractExpr)
    children = (x,)
    return new(:exp, hash(children), children, x.size)
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
  return exp(evaluate(x.children[1]))
end

exp(x::AbstractExpr) = ExpAtom(x)

function conic_form(e::ExpAtom, unique_constr)
  if !((e.head, e.children_hash) in keys(unique_constr))
    # exp(x) \leq t  <=>  (x,ones(),t) \in ExpCone
    x = e.children[1]
    t = Variable(size(x))
    objective, constraints = conic_form(x, unique_constr)
    append!(constraints, conic_form(ExpConstraint(x,ones(size(x)),t), unique_constr)[2])
    unique_constr[(e.head, e.children_hash)] = (objective, constraints)
  end
  return safe_copy(unique_constr[(e.head, e.children_hash)])
end
