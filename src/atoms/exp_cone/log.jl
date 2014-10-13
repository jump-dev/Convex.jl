#############################################################################
# log.jl
# natural logarithm of an logression
# All expressions and atoms are subtpyes of AbstractExpr.
# Please read expressions.jl first.
#############################################################################

import Base.log
export log
export sign, curvature, monotonicity, evaluate

### Logarithm

type LogAtom <: AbstractExpr
  head::Symbol
  children_hash::Uint64
  children::(AbstractExpr,)
  size::(Int64, Int64)

  function LogAtom(x::AbstractExpr)
    children = (x,)
    return new(:log, hash(children), children, x.size)
  end
end

function sign(x::LogAtom)
  return NoSign()
end

function monotonicity(x::LogAtom)
  return (Nondecreasing(),)
end

function curvature(x::LogAtom)
  return ConcaveVexity()
end

function evaluate(x::LogAtom)
  return log(evaluate(x.children[1]))
end

log(x::AbstractExpr) = LogAtom(x)

function conic_form(e::LogAtom, unique_constr)
  if !((e.head, e.children_hash) in keys(unique_constr))
    # log(z) \geq x  <=>  (x,ones(),z) \in ExpCone
    z = e.children[1]
    y = Constant(ones(size(z)))
    x = Variable(size(z))

    # x is the objective, which is the first child of ExpConstraint
    objective, constraints = conic_form(ExpConstraint(1, x, y, z), unique_constr)

    unique_constr[(e.head, e.children_hash)] = (objective, constraints)
  end
  return safe_copy(unique_constr[(e.head, e.children_hash)])
end
