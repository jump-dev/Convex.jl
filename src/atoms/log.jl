#############################################################################
# log.jl
# natural logarithm of an logression
# All expressions and atoms are subtpyes of AbstractExpr.
# Please read logressions.jl first.
#############################################################################

import Base.log
export log
export sign, curvature, monotonicity, evaluate

### Logonential

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
    # log(x) \geq t  <=>  (t,ones(),x) \in ExpCone
    x = e.children[1]
    t = Variable(size(x))
    objective, constraints = conic_form(x, unique_constr)
    append!(constraints, conic_form(ExpConstraint(t,ones(size(e.x)),x), unique_constr)[2])
    unique_constr[(x.head, x.children_hash)] = (objective, constraints)
  end
  return safe_copy(unique_constr[(x.head, x.children_hash)])
end
